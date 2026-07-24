#!/usr/bin/env python3
"""Every command in t_cmd.c must be reachable.

The dispatcher is a chain of strncmp(p, "name", n) tests, so a command
whose name starts with an earlier command's name is swallowed by it:
"mdsn1(...)" matched the four-character test for "mdsn", ran the wrong
command, and then reported a syntax error on the leftover "1(...)".
That made mdsn1 look like 238 lines of dead code for thirty years.

lsreg1 is ordered before lsreg and is fine; mdsn1 had no entry at all.
This check fails if any command is shadowed that way again.
"""
import re
import sys

src = open("tdaR/src/t_cmd.c").read()
cmds = [(m.group(1), int(m.group(2)))
        for m in re.finditer(r'strncmp\(p,"([a-z0-9_]+)",(\d+)\)', src)]

seen, bad = [], []
for name, n in cmds:
    for pname, pn in seen:
        if len(name) > pn and name[:pn] == pname[:pn]:
            bad.append((pname, name))
            break
    seen.append((name, n))

for a, b in bad:
    print('  "%s" is tested before "%s" and swallows it' % (a, b))
print("%d commands dispatched, %d shadowed" % (len(cmds), len(bad)))
sys.exit(1 if bad else 0)
