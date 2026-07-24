#!/usr/bin/env python3
"""Compare a pristine TDA source with the ported one, counting the side
effects that live inside subscripts.

    python3 tools/check_subscript_increments.py <original.c> <ours.c>

The TDAContext refactor rewrote expressions throughout.  Where the
original carried an increment in a subscript --

    AcYF[k++] /= (double)AcS[i];

-- splitting the compound assignment for the explicit casts can leave
the increment behind as a separate statement.  If the enclosing body is
brace-less it then runs once instead of once per iteration, and nothing
warns: not the compiler, not -Wall -Wextra, and not a reference test
whose output does not happen to depend on that value.  One such slip in
t_gdf.c survived 494 of them.

Comments are stripped first: a comment quoting the original expression
would otherwise be counted as code, which is how the first version of
this check reported a file as clean when it was not.

A count that is lower in the port is not proof of a fault -- the
increment may have been split out correctly, inside braces -- but every
one of them is worth reading.
"""
import os, re, sys

# The refactor moved globals onto the context, so both the array and the
# index may have picked up a "ctx->" that the original does not have --
# ESTyp[ESCnt++] became ctx->ESTyp[ctx->ESCnt++].  Matching only a bare
# \w+ on each side misses every one of those and reports the file clean.
NAME = r'(?:ctx->)?(\w+)'
SUB = re.compile(NAME + r'\s*\[\s*' + NAME + r'\s*\+\+\s*\]')
PRE = re.compile(NAME + r'\s*\[\s*\+\+\s*' + NAME + r'\s*\]')


def strip_comments(s):
    s = re.sub(r'/\*.*?\*/', ' ', s, flags=re.S)
    return re.sub(r'//[^\n]*', ' ', s)


def tally(path, pat):
    s = strip_comments(open(path, errors="replace").read())
    out = {}
    for m in pat.finditer(s):
        k = (m.group(1), m.group(2))
        out[k] = out.get(k, 0) + 1
    return out


RISK = re.compile(r'([A-Za-z_][\w.\->]*)\s*\[\s*(\+\+|--)?\s*([\w>.\-]+?)'
                  r'\s*(\+\+|--)?\s*\]\s*(\+=|-=|\*=|/=|%=)')


def risky(path):
    """Compound assignments whose subscript carries a side effect.

    These are the ones that cannot be split naively: x[i++] += y reads
    and writes one element and advances i once, while x[i] = x[i] + y
    with a separate i++ evaluates the subscript twice and, if the body
    is brace-less, advances i once per loop rather than once per pass.
    Both faults found in this port were of exactly this shape.
    """
    out = []
    for m in RISK.finditer(strip_comments(open(path, errors="replace").read())):
        if m.group(2) or m.group(4):
            out.append((m.group(1).replace("ctx->", ""),
                        (m.group(2) or "") +
                        m.group(3).replace("ctx->", "") + (m.group(4) or "")))
    return out


def sweep(odir, udir):
    """Every file in odir that has a counterpart in udir."""
    import collections
    tot, bad = 0, []
    for f in sorted(os.listdir(odir)):
        if not f.endswith(".c"):
            continue
        u = os.path.join(udir, f)
        if not os.path.exists(u):
            print("no counterpart in the port: %s" % f)
            continue
        o = risky(os.path.join(odir, f))
        tot += len(o)
        co, cu = collections.Counter(o), collections.Counter(risky(u))
        for k in co:
            if cu[k] < co[k]:
                bad.append((f, "%s[%s]" % k, co[k], cu[k]))
    print("compound assignments with a side effect in the subscript: %d" % tot)
    print("not matched one-for-one in the port: %d" % len(bad))
    for r in bad:
        print("   %-12s %-16s original %d  ported %d" % r)
    return len(bad)


def main():
    if os.path.isdir(sys.argv[1]):
        sys.exit(1 if sweep(sys.argv[1], sys.argv[2]) else 0)
    a, b = sys.argv[1], sys.argv[2]
    bad = 0
    for name, pat in (("array[idx++]", SUB), ("array[++idx]", PRE)):
        o, u = tally(a, pat), tally(b, pat)
        print("%-14s original %3d   ported %3d"
              % (name, sum(o.values()), sum(u.values())))
        for k in sorted(set(o) | set(u)):
            if u.get(k, 0) != o.get(k, 0):
                bad += 1
                print("    %-18s original %d, ported %d  <<< read both"
                      % ("%s[%s++]" % k, o.get(k, 0), u.get(k, 0)))
    print("\n%d subscript side effect(s) differ in count" % bad)


main()
