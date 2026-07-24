#!/usr/bin/env python3
"""Code that escaped its chunk, and chunk fences that do not balance.

    python3 tools/check_vignette_fences.py tdaR/vignettes/*.Rmd

Editing prose around a chunk can eat the opening ```{r} and leave the
code sitting in the markdown as plain text -- it then renders as a wall
of raw source with no output and no box, which is easy to miss in a
document this long. It can also eat the first line of the chunk, so the
code that survives refers to an object that was never built.

Reports an unbalanced fence, and any line outside a chunk that looks
like R (an assignment or a call). The second check has false positives:
a sentence ending "... structure (Box 4)." looks like a call. Read them,
do not just count them.
"""
import re, sys

bad = 0
for path in sys.argv[1:]:
    lines = open(path, errors="replace").read().split("\n")
    depth = 0
    for i, l in enumerate(lines, 1):
        if l.startswith("```{r"):
            if depth:
                print("%s:%d  a chunk opens inside another" % (path, i))
                bad += 1
            depth = 1
        elif l.rstrip() == "```":
            if not depth:
                print("%s:%d  a fence closes with none open" % (path, i))
                bad += 1
            depth = 0
    if depth:
        print("%s: the last chunk never closes" % path)
        bad += 1

    depth = 0
    for i, l in enumerate(lines, 1):
        if l.startswith("```{r"):
            depth = 1
            continue
        if l.rstrip() == "```":
            depth = 0
            continue
        if depth or l.startswith(("|", "#", ">", "    ")):
            continue
        if re.match(r'^\s*[a-zA-Z_.][\w.]*\s*(<-|\()', l):
            print("%s:%d  looks like code outside a chunk: %s"
                  % (path, i, l.strip()[:56]))
            bad += 1
print("\n%d line(s) to look at" % bad)
