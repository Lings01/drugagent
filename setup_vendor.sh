#!/bin/bash
# vendor/openfold: full aqlaboratory/openfold + the 2 CPU patches
set -euo pipefail
ROOT=/home/data/lrs/drug/drugagent
VDIR=$ROOT/drugagent/vendor
mkdir -p "$VDIR"
cd /tmp
curl -sL --max-time 600 "https://codeload.github.com/aqlaboratory/openfold/tar.gz/refs/heads/main" -o of.tar.gz
tar xzf of.tar.gz
rm -rf "$VDIR/openfold"
cp -r /tmp/openfold-main/openfold "$VDIR/openfold"
python3.12 - <<'EOF'
import pathlib
root = pathlib.Path("/home/data/lrs/drug/drugagent/drugagent/vendor/openfold")
# patch 1: utils/kernel/attention_core.py
p = root / "utils/kernel/attention_core.py"
t = p.read_text()
old = 'attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")'
new = 'try:\n    attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")\nexcept ImportError:\n    attn_core_inplace_cuda = None'
assert old in t
t = t.replace(old, new)
old2 = """        attn_core_inplace_cuda.forward_(
            attention_logits, 
            reduce(mul, attention_logits.shape[:-1]),
            attention_logits.shape[-1],
        )"""
new2 = """        if attn_core_inplace_cuda is not None:
            attn_core_inplace_cuda.forward_(
                attention_logits,
                reduce(mul, attention_logits.shape[:-1]),
                attention_logits.shape[-1],
            )
        else:
            torch.nn.functional.softmax(attention_logits, dim=-1,
                                        out=attention_logits)"""
if old2 in t:
    t = t.replace(old2, new2)
    print("patch1 fwd block applied")
p.write_text(t)
# patch 2: model/structure_module.py
p = root / "model/structure_module.py"
t = p.read_text()
old = 'attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")'
if old in t:
    new = 'try:\n    attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")\nexcept ImportError:\n    attn_core_inplace_cuda = None'
    t = t.replace(old, new)
old2 = "        if (inplace_safe):"
new2 = "        if (inplace_safe and attn_core_inplace_cuda is not None):"
if old2 in t:
    t = t.replace(old2, new2, 1)
    print("patch2 inplace guard applied")
p.write_text(t)
print("openfold patched")
EOF
echo VENDOR_DONE
