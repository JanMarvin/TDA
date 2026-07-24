#!/usr/bin/env python3
"""Invoke every command in TDA's dispatch table and check it survives.

The 350 worked examples exercise 93 of the 274 commands.  This does not test
what the other 181 compute -- it only checks that each one is still reachable
from the parser and does not crash when called without arguments.  That is a
low bar, but it is mechanical, needs no knowledge of what any command does, and
it is exactly the bar that two real bugs failed:

  conj          the command string was renamed along with the C function when
                conj() became tda_conj() to avoid the C99 complex builtin, so
                strncmp(p,"tda_conj",4) sat inside case 'c' comparing against
                "tda_" and could never match.

  rfit, rfit1   `nn` bounds a cleanup loop at the function's exit label, which
                is reached by goto before nn is assigned when parm() fails.
                Uninitialised, the loop walked AcI/AcIPtr out of bounds and
                called free() on whatever it found.  Present in the 6.4
                release, so it predates this port.

Usage: python3 tests/smoke.py src/tda [-v]
"""
import re, subprocess, sys, tempfile, os
from pathlib import Path

binary = Path(sys.argv[1] if len(sys.argv) > 1 else "./tda").resolve()
if not binary.exists():
    sys.exit("ERROR: binary %s not found" % binary)
VERBOSE = "-v" in sys.argv

# The dispatch table is a switch on the first character, so the command names
# are the string literals compared against p inside t_cmd.c.
# Skip commented-out entries (rdxf) and anything behind #if S_XWIN, which the
# portable build compiles out (xshow).
src = (Path(__file__).parent.parent / "src" / "t_cmd.c").read_text(errors="replace")
lines, keep, xwin = src.split("\n"), [], 0
for l in lines:
    if re.match(r"\s*#if\s+S_XWIN", l):
        xwin = 1
    elif re.match(r"\s*#endif", l) and xwin:
        xwin = 0
    elif not xwin and not re.match(r"\s*/\*", l):
        keep.append(l)
cmds = sorted(set(re.findall(r'strn?cmp\(p,"([a-z][a-z0-9_]*)"', "\n".join(keep))))
if not cmds:
    sys.exit("ERROR: no commands found in t_cmd.c")

work = tempfile.mkdtemp(prefix="tdasmoke")
cf = os.path.join(work, "s.cf")

unknown, crashed, hung = [], [], []
for c in cmds:
    seen_known = False
    for form in ("%s();\n", "%s;\n"):
        Path(cf).write_text(form % c)
        try:
            r = subprocess.run([str(binary), "cf=" + cf], cwd=work, timeout=20,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except subprocess.TimeoutExpired:
            hung.append(c)
            break
        # 128+n means killed by signal n: segfault, abort, bus error
        if r.returncode >= 128:
            crashed.append("%s (rc=%d, as %r)" % (c, r.returncode, form.strip()))
            break
        if b"Unknown command" not in r.stdout:
            seen_known = True
            break
    else:
        if not seen_known:
            unknown.append(c)

print("commands in dispatch table: %d" % len(cmds))
print("  reachable                : %d" % (len(cmds) - len(unknown) - len(crashed) - len(hung)))
print("  not reachable            : %d" % len(unknown))
print("  crashed                  : %d" % len(crashed))
print("  timed out                : %d" % len(hung))

if unknown and VERBOSE:
    print("\nnot reachable (may need arguments to parse, or be loop-only):")
    print("  " + " ".join(unknown))
for c in crashed:
    print("CRASH: " + c)
for c in hung:
    print("HANG:  " + c)

sys.exit(1 if crashed or hung else 0)
