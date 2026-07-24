#!/usr/bin/env python3
"""Generate tda_context_init() from non-zero file-scope initialisers."""
import re
from pathlib import Path

SKIP = {"t_xwin.c","t_xlib.c"}
VAR_DEF = re.compile(
    r"^(?:(?:unsigned|signed|short|long|const)\s+)*"
    r"(?:int|double|float|char|short|long|FILE|size_t|void)\s+"
    r"\*{0,3}([A-Za-z_][A-Za-z0-9_]*)"
    r"(?:\[[^\]]*\])*"
    r"\s*=\s*([^;{]+);"
)
ZERO = {"0","0.0","0.0f","NULL","0L","'\\0'"}

src = Path("/tmp/tda")
inits = []
for c in sorted(src.glob("*.c.orig")):
    if c.name.replace(".orig","") in SKIP: continue
    depth = 0
    for raw in c.read_text(errors="replace").splitlines():
        s = raw.strip()
        if depth == 0 and raw and raw[0] not in (' ','\t','\n','#','/','*'):
            m = VAR_DEF.match(s)
            if m:
                val = m.group(2).strip()
                if val not in ZERO and not val.startswith("{"):
                    inits.append((m.group(1), val))
        depth = max(0, depth + s.count('{') - s.count('}'))

# load global_names to filter only actual context members
ctx_h = (src / "tda_context.h").read_text()
DEFINED_CONSTS = set(re.findall(r"^#define\s+([A-Za-z_][A-Za-z0-9_]*)", ctx_h, re.MULTILINE))
_TW = {"int","double","float","char","short","long","void","FILE",
       "size_t","unsigned","signed","const","TDAContext","struct"}
global_names = set()
for mm in re.finditer(r"^\s+(\S[^{};]*);", ctx_h, re.MULTILINE):
    for nm in re.findall(r"\b([A-Za-z_][A-Za-z0-9_]+)\b", mm.group(1)):
        if nm not in _TW and nm not in DEFINED_CONSTS:
            global_names.add(nm)

lines = [
    '#include "tda_context.h"',
    '#include <math.h>',
    '',
    'void tda_context_init(TDAContext *ctx) {',
]
skipped = []
for name, val in inits:
    if name in global_names:
        lines.append(f'    ctx->{name} = {val};')
    else:
        skipped.append((name, val))

lines += ['}', '']

Path(src / "tda_context_init.c").write_text('\n'.join(lines))
print(f"wrote tda_context_init.c: {len(inits)-len(skipped)} assignments")
if skipped:
    print(f"skipped {len(skipped)} (not in struct):", [n for n,v in skipped[:5]])
