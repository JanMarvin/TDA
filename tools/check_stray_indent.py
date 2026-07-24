#!/usr/bin/env python3
"""Statements whose indentation disagrees with their block.

    python3 tools/check_stray_indent.py src/*.c

The port's TDAContext refactor rewrote expressions throughout. Where a
compound assignment carried a side effect in a subscript --

    AcYF[k++] /= (double)AcS[i];

-- rewriting it for explicit casts could leave the increment behind as a
separate statement, and a brace-less loop body then took only the first
statement:

    for (j = 0; j < ndim; ++j)
        ctx->AcYF[k] = (float)((double)ctx->AcYF[k] / (double)ctx->AcS[i]);
    k++;                     <- runs once, not ndim times

The compiler cannot see this and neither can a reference test whose
output happens not to depend on it. What gives it away is that the
stray statement sits at an indentation matching nothing around it.

Reports a brace-less `for`/`while`/`if` whose next-next line is indented
LESS than the body but is not a recognised continuation -- and, more
narrowly, any such case where that line is a bare increment.
"""
import re, sys

HEAD = re.compile(r'^(\s*)(for|while|if)\s*\(')
BARE_INC = re.compile(r'^\s*[A-Za-z_][\w.\->\[\]]*\s*(\+\+|--)\s*;\s*$')
CLOSERS = ("}", "else", "#", "/*", "*", "return", "break", "continue")


def scan(path):
    lines = open(path, errors="replace").read().split("\n")
    out = []
    for i, l in enumerate(lines):
        m = HEAD.match(l)
        if not m or l.rstrip().endswith("{") or l.rstrip().endswith(";"):
            continue
        # a head whose condition wraps over several lines is not
        # brace-less: find where the condition actually closes
        depth = l.count("(") - l.count(")")
        i2 = i
        while depth > 0 and i2 + 1 < len(lines):
            i2 += 1
            depth += lines[i2].count("(") - lines[i2].count(")")
        if i2 != i:
            if lines[i2].rstrip().endswith(("{", ";")):
                continue
            i = i2

        # brace-less head: find the single statement that is its body
        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines):
            continue
        body = lines[j]
        if body.strip().startswith("{") or body.rstrip().endswith("{"):
            continue
        bi = len(body) - len(body.lstrip())
        hi = len(m.group(1))
        if bi <= hi:
            continue
        # the statement AFTER the body
        k = j + 1
        while k < len(lines) and not lines[k].strip():
            k += 1
        if k >= len(lines):
            continue
        nxt = lines[k]
        s = nxt.strip()
        if not s or s.startswith(CLOSERS):
            continue
        ni = len(nxt) - len(nxt.lstrip())
        # Misindented but still inside the enclosing braces is cosmetic:
        # what matters is whether the statement falls OUTSIDE the block
        # the head opened.  Count braces from the head to here -- if any
        # are still open, the statement is inside something and runs with
        # it, however it is indented.
        seg = "\n".join(lines[i:k])
        if seg.count("{") > seg.count("}"):
            continue
        if ni < bi and ni != hi:
            # Rank it.  The fault this exists for is a LOOP whose body
            # subscripts an array with the very index the stray statement
            # increments: the index then advances once instead of once per
            # iteration.  A stray increment after a brace-less `if` inside
            # a loop is the common, harmless misindentation.
            inc = BARE_INC.match(nxt)
            level = "note"
            if inc:
                var = re.match(r'\s*([A-Za-z_][\w.\->]*)', s).group(1)
                loop = m.group(2) in ("for", "while")
                subscripted = re.search(r'\[\s*%s\s*[\]+]' % re.escape(var),
                                        body) is not None
                level = "SUSPECT" if (loop and subscripted) else "increment"
            out.append((k + 1, level, lines[i].strip(), body.strip(), s))
    return out


VERBOSE = "-v" in sys.argv
args = [a for a in sys.argv[1:] if a != "-v"]
total = 0
for p in args:
    for ln, level, head, body, stray in scan(p):
        if level == "note" and not VERBOSE:
            continue
        total += 1
        print("%s:%d  %s" % (p, ln, level))
        print("        %s" % head[:66])
        print("            %s" % body[:62])
        print("      >>> %s" % stray[:66])
print("\n%d flagged (run with -v for every misindented statement)" % total)
