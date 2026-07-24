#!/bin/sh
# Reproducible coverage audit for one engine file.
#   sh tools/coverage_audit.sh t_stata.c
# Builds a fresh instrumented binary (no stale-object line shifts:
# the build dir is removed first), runs the full cf suite, and prints
# the gcov percentage plus a categorized count of every missed line.
# The categories sum exactly to the missed total -- if they don't,
# the script says so and exits nonzero.
set -e
F=${1:?usage: sh tools/coverage_audit.sh t_stata.c}
cd "$(dirname "$0")/.."
rm -rf /tmp/tda-cov
sh tools/quickbuild.sh cov "$F" > /dev/null 2>&1
rm -f /tmp/tda-cov/*.gcda
python3 tests/check.py /tmp/tda-cov/tda \
    examples/exam examples/ehhnew examples/tests examples/coverage \
    > /dev/null 2>&1 || true
cd /tmp/tda-cov
gcov "$F" > /dev/null 2>&1
python3 - "$F" <<'PYEOF'
import re, sys
f = sys.argv[1]
total = hit = 0
missed = []
for ln in open(f + ".gcov"):
    m = re.match(r"\s*(#####|\d+):\s*(\d+):(.*)", ln)
    if not m:
        continue
    total += 1
    if m.group(1) == "#####":
        missed.append((int(m.group(2)), m.group(3).strip()))
    else:
        hit += 1
cats = {"error/guard/cleanup": [], "write path": [], "read path": []}
# function ranges from the source itself, not guessed line numbers
src = open(f).read().splitlines()
fnstart = {}
for i, l in enumerate(src, 1):
    m = re.match(r"(?:int|void|double)\s+([a-z0-9_]+)\(TDAContext", l)
    if m and not l.rstrip().endswith(";"):
        fnstart[i] = m.group(1)
starts = sorted(fnstart)
def owner(lno):
    o = "(file scope)"
    for s in starts:
        if s <= lno:
            o = fnstart[s]
        else:
            break
    return o
for lno, c in missed:
    if re.search(r"p_err|Error|Warning|Probably|exceeded|Can.t|goto \w+Fin"
                 r"|return\(-1\)|free\(|p_clean|!= *\d+ *\)|error", c):
        cats["error/guard/cleanup"].append(lno)
    elif owner(lno).startswith("wr_"):
        cats["write path"].append(lno)
    else:
        cats["read path"].append(lno)
print(f"{f}: {hit}/{total} lines executed = {100*hit/total:.2f}%  "
      f"({len(missed)} missed)")
ssum = 0
for k, v in sorted(cats.items(), key=lambda kv: -len(kv[1])):
    print(f"  {k:22s} {len(v):4d}")
    ssum += len(v)
if ssum != len(missed):
    print(f"CATEGORY SUM {ssum} != MISSED {len(missed)} -- audit is broken")
    sys.exit(1)
print(f"  categories sum to the missed total: {ssum} == {len(missed)}")
PYEOF
