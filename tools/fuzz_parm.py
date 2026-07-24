#!/usr/bin/env python3
"""Option-space fuzzer for TDA's command parser (t_parm) and command
surface.  Design per the project's agreed philosophy: structure from
the documentation, not random bytes -- seeds are PASSING cases from
the suite, and mutations iterate each present option through a pool
of documented-plausible values.  The oracle is (a) no crash (signal),
(b) no hang (timeout), (c) after the mutated command file, a known
-good file still runs (engine state integrity).  Ordinary command
errors are expected and only counted.  Usage:
    python3 tools/fuzz_parm.py [seed] [max_runs]
"""
import random, re, shutil, subprocess, sys, tempfile
from pathlib import Path

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 42
MAXR = int(sys.argv[2]) if len(sys.argv) > 2 else 300
random.seed(SEED)
ROOT = Path(__file__).resolve().parent.parent
import os
TDA = Path(os.environ.get("TDA_FUZZ_BIN", ROOT / "src" / "tda"))
CASES = sorted((ROOT / "examples" / "tests").glob("*.cf"))
POOL = ["0", "1", "2", "3", "5", "-1", "100", "0.5", "1.e-8", "1.e8",
        "10.4", "12.8", "1,2", "0,1"]
CANARY = "nvar(A[4.0]=1); dstat;\n"
opt_pat = re.compile(r"\b([a-z][a-z0-9]{1,7})=([^,()\s;]+)")

def run(workdir, cf):
    try:
        p = subprocess.run([str(TDA), f"cf={cf}"], cwd=workdir,
                           capture_output=True, timeout=20)
        if p.returncode < 0:
            return ("crash", -p.returncode)
        # TDA exits 0 even on command errors; classify from the output
        out = p.stdout.decode(errors="replace")
        return ("ok", 1 if ("Error" in out or "error" in out) else 0)
    except subprocess.TimeoutExpired:
        return ("hang", None)

n = crashes = hangs = errors = 0
findings = []
while n < MAXR:
    case = random.choice(CASES)
    text = case.read_text()
    opts = list(opt_pat.finditer(text))
    if not opts:
        continue
    m = random.choice(opts)
    val = random.choice(POOL)
    mutated = text[:m.start(2)] + val + text[m.end(2):]
    wd = tempfile.mkdtemp(prefix="fzp")
    for f in case.parent.glob("*"):
        if f.suffix != ".ref" and f.is_file():
            shutil.copy(f, wd)
    Path(wd, "mut.cf").write_text(mutated)
    Path(wd, "canary.cf").write_text(CANARY)
    st, rc = run(wd, "mut.cf")
    n += 1
    desc = f"{case.name}: {m.group(1)}={m.group(2)} -> {val}"
    if st == "crash":
        crashes += 1; findings.append(("CRASH sig" + str(rc), desc, wd)); continue
    if st == "hang":
        hangs += 1; findings.append(("HANG", desc, wd)); continue
    if rc != 0:
        errors += 1
    st2, rc2 = run(wd, "canary.cf")
    if st2 != "ok" or rc2 != 0:
        findings.append(("STATE after", desc, wd)); continue
    shutil.rmtree(wd, ignore_errors=True)
for kind, desc, wd in findings:
    print(f"{kind:12s} {desc}   [{wd}]")
print(f"fuzz_parm seed {SEED}: {n} runs, {errors} clean command errors, "
      f"{crashes} crashes, {hangs} hangs, {len(findings)} findings")
sys.exit(1 if findings else 0)
