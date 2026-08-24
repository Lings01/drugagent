"""Agent 'brain': robust wrapper around an OpenAI-compatible LLM endpoint.

Handles llama.cpp quirks (reasoning_content field, empty content) and
records every decision to the project decision log.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from loguru import logger

from .config import LLMConfig
from .utils import jload, jsave


@dataclass
class Decision:
    """A recorded agent judgment."""
    node: str
    question: str
    answer: Any
    rationale: str = ""
    model: str = ""
    raw: str = ""

    def to_dict(self) -> dict:
        return {
            "node": self.node, "question": self.question, "answer": self.answer,
            "rationale": self.rationale, "model": self.model, "raw": self.raw[:2000],
        }


class AgentBrain:
    """Thin, fault-tolerant LLM client with decision logging."""

    def __init__(self, cfg: LLMConfig | None = None, project_dir: str | Path | None = None):
        self.cfg = cfg or LLMConfig.from_env()
        from openai import OpenAI  # lazy import (langchain-openai dep)
        self.client = OpenAI(base_url=self.cfg.base_url, api_key=self.cfg.api_key,
                             timeout=self.cfg.timeout)
        self.project_dir = Path(project_dir) if project_dir else None
        if self.project_dir:
            self.project_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    def _chat(self, system: str, user: str, temperature: float | None = None,
              max_tokens: int | None = None) -> str:
        resp = self.client.chat.completions.create(
            model=self.cfg.model,
            messages=[{"role": "system", "content": system},
                      {"role": "user", "content": user}],
            temperature=self.cfg.temperature if temperature is None else temperature,
            max_tokens=max_tokens or self.cfg.max_tokens,
        )
        msg = resp.choices[0].message
        content = (msg.content or "").strip()
        reasoning = getattr(msg, "reasoning_content", "") or ""
        if not content and reasoning:
            # thinking model exhausted the token budget inside reasoning;
            # take the tail of the reasoning as a best-effort answer
            content = reasoning[-1500:]
        return content

    # ------------------------------------------------------------------ #
    def chat_tools(self, messages: list[dict], tools: list[dict] | None = None,
                   temperature: float | None = None,
                   max_tokens: int | None = None) -> dict:
        """Native function-calling chat (ReAct loop).

        Returns {"content": str, "tool_calls": [{"id","name","args"}] | None}.
        """
        resp = self.client.chat.completions.create(
            model=self.cfg.model,
            messages=messages,
            tools=tools or None,
            temperature=self.cfg.temperature if temperature is None else temperature,
            max_tokens=max_tokens or self.cfg.max_tokens,
        )
        m = resp.choices[0].message
        content = (m.content or "").strip()
        reasoning = getattr(m, "reasoning_content", "") or ""
        tool_calls = []
        for tc in (getattr(m, "tool_calls", None) or []):
            fn = tc.function
            try:
                args = json.loads(fn.arguments or "{}")
                if not isinstance(args, dict):
                    args = {"value": args}
            except json.JSONDecodeError:
                args = {"_raw_arguments": fn.arguments}
            tool_calls.append({"id": getattr(tc, "id", None) or
                               f"call_{len(tool_calls)}",
                               "name": fn.name, "args": args})
        if not content and not tool_calls and reasoning:
            content = reasoning[-800:]
        return {"content": content,
                "tool_calls": tool_calls or None}

    # ------------------------------------------------------------------ #
    def ask(self, node: str, question: str, context: str = "", *,
            system: str = "You are a computational drug-discovery expert. "
                          "Be concise. Answer in Chinese unless asked otherwise.",
            temperature: float | None = None) -> str:
        user = question if not context else f"Context:\n{context}\n\nQuestion: {question}"
        try:
            answer = self._chat(system, user, temperature)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"LLM call failed in {node}: {e}; falling back to empty")
            answer = ""
        self._log(node, question, answer, answer[:200])
        return answer

    # ------------------------------------------------------------------ #
    def decide(self, node: str, question: str, context: str = "", *,
               choices: list[str] | None = None,
               expect: str = "text") -> Decision:
        """Structured judgment.

        expect: 'text' | 'json' | 'choice'
        Falls back to `choices[0]` (or empty text) when parsing fails.
        """
        sys_prompt = (
            "You are a senior computational drug-discovery scientist making a "
            "specific judgment call. Be decisive and concise. "
            + ('Respond ONLY with a JSON object like {"value": ..., "rationale": "..."}. '
               if expect in ("json", "choice") else
               "Respond in plain text, at most 3 sentences. ")
        )
        if expect == "choice" and choices:
            sys_prompt += f"\"value\" MUST be one of: {json.dumps(choices, ensure_ascii=False)}."
        user = question
        if context:
            user = f"Context:\n{context}\n\nDecision needed: {question}"
        try:
            raw = self._chat(sys_prompt, user)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"LLM decide() failed in {node}: {e}")
            raw = ""

        value, rationale = self._parse(raw, expect, choices)
        dec = Decision(node=node, question=question, answer=value,
                       rationale=rationale, model=self.cfg.model, raw=raw)
        self._log(node, question, value, rationale, raw)
        return dec

    @staticmethod
    def _parse(raw: str, expect: str, choices: list[str] | None):
        if not raw:
            return (choices[0] if choices else ""), "LLM unavailable; default used"
        # try to extract JSON
        m = re.search(r"\{.*\}", raw, re.S)
        if m and expect in ("json", "choice"):
            try:
                d = json.loads(m.group(0))
                value = d.get("value", d)
                rationale = str(d.get("rationale", d.get("reason", "")))
                if expect == "choice" and choices:
                    value = str(value)
                    if value not in choices:
                        low = value.lower()
                        match = next((c for c in choices if c.lower() == low), None)
                        value = match if match else choices[0]
                return value, rationale
            except json.JSONDecodeError:
                pass
        text = raw.strip()
        if expect == "choice" and choices:
            low = text.lower()
            match = next((c for c in choices if c.lower() in low), None)
            if match:
                return match, text
            return choices[0], text
        return text, ""

    # ------------------------------------------------------------------ #
    def _log(self, node: str, question: str, answer: Any, rationale: str,
             raw: str = "") -> None:
        dec = Decision(node=node, question=question, answer=answer,
                       rationale=rationale, model=self.cfg.model, raw=raw)
        logger.info(f"[decision@{node}] {question[:120]} -> {str(answer)[:120]}")
        if self.project_dir:
            log_path = self.project_dir / "decisions.json"
            log = jload(log_path) if log_path.exists() else []
            log.append(dec.to_dict())
            jsave(log_path, log)

    # ------------------------------------------------------------------ #
    @staticmethod
    def available(cfg: LLMConfig | None = None) -> bool:
        """Ping the endpoint cheaply."""
        cfg = cfg or LLMConfig.from_env()
        try:
            import requests
            r = requests.get(f"{cfg.base_url.rstrip('/')}/models",
                             headers={"Authorization": f"Bearer {cfg.api_key}"}, timeout=10)
            return r.status_code == 200
        except Exception:  # noqa: BLE001
            return False
