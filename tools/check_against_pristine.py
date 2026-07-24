#!/usr/bin/env python3
"""Run every shipped command file through both binaries and diff.

    python3 tools/check_against_pristine.py <pristine_tda> <our_tda> <examples_dir>

Answers "is this difference ours?" for the whole example suite at once,
instead of one command file at a time. Any difference is either a fault
we introduced or a fix we made deliberately; both want knowing about.

Console output is compared after stripping what legitimately varies:
timestamps, memory figures, and the version banner.
"""
import os, re, subprocess, sys, tempfile, shutil

# absolute: each run happens in a scratch directory, so a relative path
# to either binary silently fails there and every file looks different
PRIS, OURS, EX = (os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]),
                  os.path.abspath(sys.argv[3]))

VARY = [
    (re.compile(r'^.*Analysis of Transition Data.*$', re.M), ''),
    (re.compile(r'^.*(Current|Max) memory.*$', re.M), ''),
    (re.compile(r'^.*End of program.*$', re.M), ''),
    (re.compile(r'\b\d{2}:\d{2}:\d{2}\b'), ''),
]


def norm(t):
    for pat, rep in VARY:
        t = pat.sub(rep, t)
    return [l.rstrip() for l in t.split('\n') if l.strip()]


def run(binary, cf, srcdir):
    d = tempfile.mkdtemp()
    try:
        for f in os.listdir(srcdir):
            p = os.path.join(srcdir, f)
            if os.path.isfile(p):
                shutil.copy(p, d)
        r = subprocess.run([binary, 'cf=' + cf], cwd=d, capture_output=True,
                           text=True, timeout=120)
        return r.stdout
    except Exception as e:
        return 'RUNFAIL %s' % e   # never silently equal
    finally:
        shutil.rmtree(d, ignore_errors=True)


cfs = sorted(f for f in os.listdir(EX) if f.endswith('.cf'))
diff = []
for cf in cfs:
    a, b = norm(run(PRIS, cf, EX)), norm(run(OURS, cf, EX))
    if a != b:
        n = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
        diff.append((cf, n))
print('command files run: %d' % len(cfs))
print('differing: %d' % len(diff))
for cf, n in sorted(diff, key=lambda t: -t[1]):
    print('  %-14s %d line(s)' % (cf, n))
