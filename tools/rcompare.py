import sys, re
from pathlib import Path
src = open('../../src/check.py').read().split('def run_suite')[0]
src = src.replace('binary = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./tda")', 'binary=None')
src = src.replace('if not binary.exists():\n    sys.exit("ERROR: binary %s not found" % binary)', '')
exec(src)

tot_ok = tot = 0
for suite in ("exam", "ehhnew"):
    refdir = Path('../../examples')/suite
    outdir = Path('rout')/suite
    ok, fails = 0, []
    for ref in sorted(refdir.glob('*.ref')):
        rout = outdir/(ref.stem + '.rout')
        if not rout.exists():
            continue
        got = lines_of(strip_residue(rout.read_text(errors='replace')))
        want = lines_of(strip_residue(ref.read_text(errors='replace')))
        rtol = RTOL_OVERRIDE.get(ref.stem + '.cf', RTOL)
        sf = (ref.stem + '.cf') in SIGN_FREE
        if lines_match(want, got, rtol, sf):
            ok += 1
        else:
            fails.append((ref.stem, first_diff(want, got, rtol, sf)))
    print("%s: %d/%d" % (suite, ok, ok+len(fails)))
    for nm, d in fails[:12]:
        print("  FAIL %-12s line %d\n    want %r\n    got  %r" % (nm, d[0], d[1], d[2]))
    tot_ok += ok; tot += ok+len(fails)
print("Total: %d/%d" % (tot_ok, tot))
