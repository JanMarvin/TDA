#!/usr/bin/env python3
"""No internal helper may sit directly under a roxygen block.

roxygen attaches a #' block to whatever definition FOLLOWS it.  Insert a
`.helper <- function(...)` between the block and the exported function
it documents, and the docs -- examples included -- silently move onto
the helper: the exported function loses its page, and R CMD check then
tries to RUN those examples against an internal function that takes
different arguments.  It fails as an ERROR in "checking examples", a
long way from the edit that caused it.

This happened twice while adding export-first constructors
(.ineq_table, .dstat_from_exports, .seq_info_from_exports), so it is
worth a gate rather than vigilance.

    python3 tools/check_roxygen_helpers.py     # rc 1 if any found
"""
import re, sys, pathlib

bad = []
for p in sorted(pathlib.Path("tdaR/R").glob("*.R")):
    lines = p.read_text().split("\n")
    for i, l in enumerate(lines):
        m = re.match(r"^(\.[A-Za-z0-9_.]+)\s*<-\s*function", l)
        if not m:
            continue
        k = i - 1
        while k >= 0 and (lines[k].strip() == "" or
                          (lines[k].lstrip().startswith("#") and
                           not lines[k].lstrip().startswith("#'"))):
            k -= 1
        if k >= 0 and lines[k].lstrip().startswith("#'"):
            bad.append((str(p), i + 1, m.group(1)))

if bad:
    print("internal helpers sitting under a roxygen block:")
    for f, ln, nm in bad:
        print("    %s:%d  %s" % (f, ln, nm))
    print("move each BELOW the function its block documents")
    sys.exit(1)
print("no internal helper is adopting an exported function's docs")
