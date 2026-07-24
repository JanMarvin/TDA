#!/usr/bin/env python3
"""Enumerate every printer in the TDA sources that formats a number, and
whether the function containing it has an export producer.

The NUMBER RULE: every printer that formats a number gets a producer.
Only zero-number printers (rulers, padding, pure text, name echoes) are
excluded, and that property is checkable rather than asserted -- this
script is the check.  Run from the repository root:

    python3 tools/export_census.py [--detail]

Counts functions, not call sites: one function's table is one producer.
"""
import os, re, sys

SRC = "tdaR/src"
# a print call that formats at least one number: a %-conversion other
# than %s/%c, or one of TDA's own format-string variables
NUMFMT = re.compile(r'%[-#0-9.*+ ]*(?:l|h|ll)?[diouxXeEfgGa]')
FMTVAR = re.compile(r'(PM\w*FmtS|VPFmtS|PMNFmtS|PMTFmtS|PMFmtS)')
PRINTC = re.compile(r'\b(printf1|printfe|fprintf|prn_sfmt|prn_f1mat)\s*\(')
# printf1/printfe are tapped in t_gen.c; fprintf to a context-held
# stream is redirected to tda_fprintf() by tda.h.  Both export every
# number they format without a producer at the call site.
# prn_sfmt() and prn_f1mat() are not extra entry points: they format
# through printf1/fprintf themselves, so their numbers reach the taps
# too -- verified in t_pgen.c, not assumed from the names.
TAPPED = re.compile(r'\b(printf1|printfe|fprintf|prn_sfmt|prn_f1mat)\s*\(')
WINDOW = 12   # lines between a print and its producer
FUNCDEF = re.compile(r'^[A-Za-z_][A-Za-z0-9_ *]*\b(\w+)\s*\(TDAContext\b[^;]*$')

def functions(path):
    """Yield (name, start, end, lines) for each function definition."""
    lines = open(path, errors="replace").read().split("\n")
    starts = []
    for i, l in enumerate(lines):
        m = FUNCDEF.match(l)
        if m and not l.rstrip().endswith(";"):
            # a definition is followed by '{' on this or the next lines
            for k in range(i, min(i + 4, len(lines))):
                if lines[k].strip() == "{":
                    starts.append((m.group(1), i))
                    break
    for n, (name, s) in enumerate(starts):
        e = starts[n + 1][1] if n + 1 < len(starts) else len(lines)
        yield name, s + 1, e, lines[s:e]

def main():
    detail = "--detail" in sys.argv
    rows = []
    for fn in sorted(os.listdir(SRC)):
        if not fn.endswith(".c") or fn.startswith("tda_"):
            continue
        path = os.path.join(SRC, fn)
        for name, start, end, body in functions(path):
            text = "\n".join(body)
            sites = [i for i, l in enumerate(body)
                     if PRINTC.search(l) and
                     (NUMFMT.search(l) or FMTVAR.search(l))]
            if not sites:
                continue
            # A site counts as covered when a producer call sits within
            # WINDOW lines of it -- the producer rule puts the export
            # BESIDE the print, so proximity is the checkable form of
            # that.  Marking a whole function covered because it
            # exports something somewhere flattered every partially
            # switched printer.
            exp = [i for i, l in enumerate(body)
                   if re.search(r'tda_export_(mat|row|cell|strings|str_row)\b', l)]
            ncov = sum(1 for i in sites
                       if any(abs(i - e) <= WINDOW for e in exp))
            # printf1()/printfe() are tapped centrally in t_gen.c, so
            # every number they print is exported already -- generically
            # (line, value), not as a named table with a shape.  Counted
            # apart from the named producers so neither figure flatters
            # the other.
            ntap = sum(1 for i in sites if TAPPED.search(body[i]))
            rows.append((fn, name, start, len(sites), ncov, ntap))

    nsite = sum(r[3] for r in rows)
    ncov = sum(r[4] for r in rows)
    ntap = sum(r[5] for r in rows)
    either = sum(max(r[4], r[5]) for r in rows)
    full = [r for r in rows if r[4] == r[3]]
    part = [r for r in rows if 0 < r[4] < r[3]]
    none = [r for r in rows if r[4] == 0]
    pct = lambda k: 100.0 * k / nsite if nsite else 0.0
    print("number-printing functions: %d" % len(rows))
    print("  fully covered (named):   %d" % len(full))
    print("  partly covered (named):  %d" % len(part))
    print("  no named producer:       %d" % len(none))
    print("  sites, named producer:   %d of %d  (%.1f%%)" % (ncov, nsite, pct(ncov)))
    print("  sites, generic tap:      %d of %d  (%.1f%%)" % (ntap, nsite, pct(ntap)))
    print("  sites, either:           %d of %d  (%.1f%%)"
          % (either, nsite, pct(either)))
    print("  sites, neither:          %d" % (nsite - either))
    if detail:
        print("\n-- uncovered sites, by count --")
        print("(uncovered = neither a named producer nor the console tap)")
        for fn, name, start, n, c, t in sorted(
                rows, key=lambda r: -(r[3] - max(r[4], r[5]))):
            u = n - max(c, t)
            if u <= 0:
                continue
            print("  %-14s %-20s %s:%d  (%d of %d uncovered)"
                  % (fn, name, fn, start, u, n))
    return 0

if __name__ == "__main__":
    sys.exit(main())
