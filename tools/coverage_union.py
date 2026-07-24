#!/usr/bin/env python3
"""Union coverage across the standalone and package profiles.

Standalone profile: /home/claude/covbuild/standalone -- built there as
`cc --coverage -O0 -g -o tda *.c -lm`, which names the objects
tda-<stem>.o; the notes/data files follow (tda-<stem>.gcno/.gcda), so
gcov must be addressed with `gcov -t tda-<stem>` from that directory.
Addressing it as `gcov -o . <stem>.c` fails silently per file
("cannot open notes file") and yields a union table that is really
just the package side -- that cost half a session to notice.

Package profile: tdaR/src -- R CMD INSTALL with --coverage Makevars
compiles in place with plain object names, so `gcov -t -o . <stem>.c`
works there.  Run tools/coverage.sh with COV_DIR pointed AWAY from
/home/claude/covbuild or it will rm -rf the standalone profile.

R-only sources (tda_zoo.c and friends) exist only in the package
profile and are counted from that side alone.
"""
import subprocess, pathlib, re

SA = pathlib.Path("/home/claude/covbuild/standalone")
PK = pathlib.Path(__file__).resolve().parents[1] / "tdaR" / "src"

def parse(txt):
    execd, execable = set(), set()
    for ln in txt.splitlines():
        m = re.match(r"\s*([0-9#=-]+[*]?):\s*(\d+):", ln)
        if not m: continue
        cnt, no = m.group(1).rstrip("*"), int(m.group(2))
        if no == 0 or cnt == "-": continue
        execable.add(no)
        if cnt not in ("#####", "====="):
            execd.add(no)
    return execd, execable

def lines_sa(stem):
    if not (SA / f"tda-{stem}.gcda").exists(): return set(), set()
    r = subprocess.run(["gcov", "-t", f"tda-{stem}"], cwd=SA,
                       capture_output=True, text=True)
    return parse(r.stdout)

def lines_pk(stem):
    if not (PK / f"{stem}.gcda").exists(): return set(), set()
    r = subprocess.run(["gcov", "-t", "-o", ".", f"{stem}.c"], cwd=PK,
                       capture_output=True, text=True)
    return parse(r.stdout)

def main():
    stems = sorted(p.stem for p in SA.glob("*.c"))
    stems += sorted(set(p.stem for p in PK.glob("*.c")) - set(stems))
    rows, tot_e, tot_a = [], 0, 0
    for st in stems:
        e1, a1 = lines_sa(st)
        e2, a2 = lines_pk(st)
        a = a1 | a2; e = e1 | e2
        if not a: continue
        rows.append((100.0*len(e)/len(a), len(a), st + ".c"))
        tot_e += len(e); tot_a += len(a)
    print(f"UNION aggregate: {100.0*tot_e/tot_a:.2f}% of {tot_a} lines")
    print(f"    pct  lines  file   (all {len(rows)} files, ascending)")
    for pct, na, n in sorted(rows):
        print(f"{pct:7.2f}% {na:6d}  {n}")

if __name__ == "__main__":
    main()
