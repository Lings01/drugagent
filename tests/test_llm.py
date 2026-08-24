from drugagent.llm import AgentBrain


def test_parse_json():
    val, rat = AgentBrain._parse(
        '好的。{"value": 12.5, "rationale": "分布合理"} 就这样。', "json", None)
    assert val == 12.5
    assert "分布合理" in rat


def test_parse_choice_match():
    val, _ = AgentBrain._parse(
        '我认为应该选 keep_ligand, 因为配体明确。', "choice",
        ["clean_and_use", "keep_ligand", "remodel"])
    assert val == "keep_ligand"


def test_parse_choice_fallback():
    val, _ = AgentBrain._parse("完全看不懂", "choice", ["a", "b"])
    assert val == "a"


def test_parse_empty():
    val, rat = AgentBrain._parse("", "choice", ["x", "y"])
    assert val == "x"
    assert "default" in rat


def test_parse_text():
    val, _ = AgentBrain._parse("直接对接即可, 库不大。", "text", None)
    assert "对接" in val


def test_available(has_net):
    if not has_net:
        return
    assert AgentBrain.available() is True
