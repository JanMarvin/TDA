#!/usr/bin/env python3
"""Control-flow keywords, a pristine TDA source tree against the port.

    python3 tools/check_control_flow.py <pristine_dir> src

A mechanical rewrite can drop a `break`, an `else` or a `goto` without
the compiler noticing. This counts them per file and reports only
FEWER-in-the-port: more is routine, since the port adds guards and
error paths throughout.

Strip character literals BEFORE string literals. A '"' inside a char
literal otherwise opens a string that swallows the rest of the file --
the first run of this reported t_gen.c as having lost eleven breaks
when the port in fact has twenty.
"""
import os, re, sys

KW = ("break", "continue", "goto", "return", "else", "switch", "case")


def strip(s):
    s = re.sub(r'/\*.*?\*/', ' ', s, flags=re.S)
    s = re.sub(r'//[^\n]*', ' ', s)
    s = re.sub(r"'(?:[^'\\]|\\.)'", "'c'", s)
    return re.sub(r'"(?:[^"\\]|\\.)*"', '""', s)


def tally(p):
    s = strip(open(p, errors="replace").read())
    return {k: len(re.findall(r'\b%s\b' % k, s)) for k in KW}


odir, udir = sys.argv[1], sys.argv[2]
bad = 0
for f in sorted(os.listdir(odir)):
    if not f.endswith(".c"):
        continue
    u = os.path.join(udir, f)
    if not os.path.exists(u):
        continue
    o, n = tally(os.path.join(odir, f)), tally(u)
    d = {k: (o[k], n[k]) for k in KW if n[k] < o[k]}
    if d:
        bad += 1
        print("  %-12s %s" % (f, "  ".join("%s %d->%d" % (k, a, b)
                                           for k, (a, b) in d.items())))
print("\n%d file(s) with a control-flow keyword fewer in the port" % bad)
