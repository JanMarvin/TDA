#!/usr/bin/env python3
"""Coverage census: which dispatched commands have no .cf case at all.

The corpus is examples/exam, examples/ehhnew, examples/tests and
examples/coverage.  A command counts as exercised if it appears at the
start of a statement in any .cf file there ("name(" or "name;" or
"name ="), which is how TDA's own parser sees it.

For every unexercised command the script also digs out the syntax box
from its own source file's header comment, so a batch of cases can be
written from one listing instead of one grep per command.

    python3 tools/coverage_census.py            the listing
    python3 tools/coverage_census.py --count    just the numbers
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# every dispatched command, in dispatch order
src = (ROOT / "src" / "t_cmd.c").read_text(errors="replace")
# strip comments first: a dispatch line commented out (rdxf) is not
# a dispatched command
src = re.sub(r"/\*.*?\*/", "", src, flags=re.DOTALL)
cmds = []
for m in re.finditer(r'strncmp\(p,"([a-z0-9_#]+)",(\d+)\)', src):
    if m.group(1) not in [c for c, _ in cmds]:
        cmds.append((m.group(1), int(m.group(2))))

# which appear in the corpus
corpus = []
for suite in ("exam", "ehhnew", "tests", "coverage"):
    corpus += list((ROOT / "examples" / suite).glob("*.cf"))
text = "\n".join(p.read_text(errors="replace") for p in corpus)
# strip comments
text = re.sub(r"#[^\n]*", "", text)

used = set()
for name, _ in cmds:
    if re.search(r"(?m)^\s*%s\s*[(;=]" % re.escape(name), text):
        used.add(name)

# commands that cannot be exercised in this build: the X11 window
# commands need a screen and live window state (xshow's stub refusal IS
# pinned, in xstub.cf); rdxf is commented out of the dispatcher.
UNTESTABLE = {"xconh", "xreg", "xplotf"}

missing = [name for name, _ in cmds
           if name not in used and name not in UNTESTABLE]

if "--count" in sys.argv:
    print("%d dispatched, %d exercised, %d without any case"
          % (len(cmds), len(used), len(missing)))
    sys.exit(0)

# locate each missing command's syntax box: the comment block above its
# own "int <name>(TDAContext *ctx)" definition, in whichever t_*.c holds it
def syntax_box(name):
    for f in sorted((ROOT / "src").glob("t_*.c")):
        s = f.read_text(errors="replace")
        m = re.search(r"int %s\(TDAContext \*ctx\)\s*\n\{" % re.escape(name), s)
        if not m:
            continue
        # TDA writes its boxes one /* ... */ comment per LINE, so walk
        # upward collecting consecutive comment lines
        lines = s[:m.start()].splitlines()
        box = []
        for line in reversed(lines):
            t = line.strip()
            if not t or t.startswith("/*") or t.startswith("*") \
               or t.endswith("*/"):
                box.append(line)
                if re.match(r"/\* -{10,}", t):
                    break
            else:
                break
        return f.name, "\n".join(reversed(box)).strip("\n")
    return None, ""

print("%d dispatched, %d exercised, %d without any case\n"
      % (len(cmds), len(used), len(missing)))
byfile = {}
for name in missing:
    f, box = syntax_box(name)
    byfile.setdefault(f or "?", []).append((name, box))

for f in sorted(byfile, key=lambda k: -len(byfile[k])):
    print("== %s: %d commands" % (f, len(byfile[f])))
    for name, box in byfile[f]:
        print("-- %s" % name)
        for line in box.splitlines():
            line = line.strip().strip("/*").rstrip("*/").rstrip()
            if line.strip():
                print("   %s" % line)
        print()
