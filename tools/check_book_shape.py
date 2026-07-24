#!/usr/bin/env python3
"""Does each ehhnew box print what its command file's .ref file holds?

    python3 tools/check_book_shape.py <html> [ref_dir] [qa_csv]

check_vignette_refs.py pulls the ESTIMATE ROWS out of a .ref file --
"Idx SN Org Des MT Variable Coeff Error C/Error" and the lsreg/glm
shape -- and asks whether the vignette box reproduces those numbers.
That is a real check of the coefficients and blind to everything else:
a survivor table, a residual listing, a frequency count that TDA prints
and the box does not show is invisible to it, exactly as the manual's
membership test was blind to a half-missing box.

This counts instead, over the whole .ref file, with the same exclusions
the manual shape check uses: the banner and settings TDA echoes before
it estimates, prose, and figures' axis runs. Each occurrence in the .ref
must be consumed by a distinct occurrence in the box, exact matches
first and the leftovers retried at four significant figures.

Exact equality is not the bar. TDA's console carries per-iteration
traces and memory figures a vignette would never print, so read the
list; the useful signal is a box short by a whole table's worth.
"""
import collections, csv, html as H, os, re, sys

HTML = sys.argv[1]
REF = sys.argv[2] if len(sys.argv) > 2 else "examples/ehhnew"
QA = sys.argv[3] if len(sys.argv) > 3 else "tests/vignette-qa/vignette-qa-book.csv"

NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
# TDA's console furniture: what it says about reading files, building
# its data matrix and its own settings, as opposed to what it reports as
# a result. None of it belongs in a vignette box.
SKIP = re.compile(
    r"TDA\.|Current memory|Max memory|End of program|Reading |^\s*$|"
    r"Algorithm|Number of model parameters|Type of covariance|"
    r"Maximum number of|Convergence criterion|Tolerance|Mue of |"
    r"Minimum of step|Scaling factor|Method:|Model:|Idx |^-+$|^=+$|"
    r"Read records|Number of cases|Number of variables|Free format|"
    r"Using data file|Missing values|Creating|Created|single episode|"
    r"Number of episodes|Sum of|^\s*nvar\(|^\s*edef\(|Maximum of log|"
    r"Norm of final|Last absolute|Last relative|Numerical problem|"
    # TDA prints a Sum row under each spell's transitions and an episode
    # total under the table. Both are sums of the rows themselves --
    # tapply(episodes, sn, sum) gives 201, 162, 107, 62 and 532 for
    # ehd7, exactly as it prints them -- so a box that shows the table
    # has not lost them.
    r"^Sum\s|Number of episodes",
    re.I)


def decimals(t):
    m = re.match(r'^-?\d*(?:\.(\d*))?(?:[eE]([-+]?\d+))?$', t)
    if not m:
        return 0
    return len(m.group(1) or "") - int(m.group(2) or 0)


def ref_values(cf):
    path = os.path.join(REF, cf + ".ref")
    if not os.path.exists(path):
        return None
    text = open(path, errors="replace").read()
    # A command file that opens a PostScript file draws a figure, and
    # its box shows the figure. What its .ref holds besides is setup --
    # the case selection, the file it opened -- not a table anyone would
    # print.
    if "psfile=" in text or "Opened new PostScript file" in text:
        return None
    out = []
    # Two blocks in a .ref are TDA talking about its own progress rather
    # than reporting a result, and no vignette box prints either: the
    # per-iteration trace ("Iter / Function Value / Norm of Gradient"),
    # and the echo of the data definition ("Idx Variable T S PFmt").
    skipping = False
    for line in open(path, errors="replace"):
        if re.search(r"^\s*Iter\s+Function Value|^Idx Variable", line):
            skipping = True
            continue
        if skipping:
            if re.match(r"^=+\s*$", line) or re.match(r"^\s*$", line):
                skipping = False
            continue
        if SKIP.search(line):
            continue
        if len(re.findall(r"[A-Za-z]{2,}", line)) >= 8:
            continue
        # An estimate row is "Idx SN Org Des MT Variable Coeff Error
        # C/Error Signif". The first five are bookkeeping -- the row's
        # number, and which transition it belongs to, constant for a
        # single-transition model -- and summary() carries them as row
        # names instead. Signif is TDA's 1 - p where summary() prints p.
        # Count the three that are the estimate.
        m = re.match(r"\s*\d+\s+\d+\s+\d+\s+\d+\s+\w\s+\S+\s+"
                     r"(-?[\d.]+(?:[eE][-+]?\d+)?)\s+"
                     r"(-?[\d.]+(?:[eE][-+]?\d+)?)\s+"
                     r"(-?[\d.]+(?:[eE][-+]?\d+)?)", line)
        if m:
            out += [(float(g), decimals(g)) for g in m.groups()]
            continue
        out += [(float(t), decimals(t)) for t in NUM.findall(line)]
    return out


def boxes(path):
    doc = open(path, errors="replace").read()
    out = {}
    for m in re.finditer(r'<span class="lab">Box (\d+)</span>(.*?)</p>', doc, re.S):
        st = doc.find('<div class="box">', m.end())
        if st < 0:
            continue
        d, i = 0, st
        for t in re.finditer(r"<(/?)div\b", doc[st:]):
            d += -1 if t.group(1) else 1
            if d == 0:
                i = st + t.end()
                break
        out[int(m.group(1))] = H.unescape(re.sub("<[^>]+>", " ", doc[st:i]))
    return out


# Differences looked into and explained, listed apart rather than left
# at the top of a queue of work.
SETTLED = {
    "ehf4": "ehf4.cf splits inside edef(split=) and so reports the "
            "episodes BEFORE the split (600: 142/458); the box splits "
            "first with tda_split() and reports the split ones. The "
            "coefficients match the .ref exactly, all nine rows.",
    "ehf3": "same route difference as ehf4",
}

bx = boxes(HTML)
rows, settled = [], []
# the book table's columns: book box, page, caption, command file,
# vignette box, verdict, note
for r in list(csv.reader(open(QA)))[1:]:
    if len(r) < 5 or not r[4].strip():
        continue
    cfs = re.findall(r'\b([a-z]{1,6}\d*[a-z]?)(?:\.cf)?\s*$', r[3].strip())
    if not cfs:
        cfs = re.findall(r'\b([a-z]{1,6}\d*[a-z]?)\.cf', r[2])
    nums = [int(x) for x in re.findall(r"\d+", r[4])]
    if not cfs or not nums:
        continue
    pool = []
    for n in nums:
        pool += [float(t) for t in NUM.findall(bx.get(n, ""))]
    want = []
    for cf in cfs:
        v = ref_values(cf)
        if v:
            want += v
    if not want:
        continue
    left = sorted(pool)
    leftover = []
    for x, d in want:
        tol = 0.5 * 10 ** (-d) + 1e-12
        hit, best = None, None
        for k, y in enumerate(left):
            dy = abs(y - x)
            if dy <= tol and (best is None or dy < best):
                hit, best = k, dy
        if hit is None:
            leftover.append(x)
        else:
            del left[hit]
    short = 0
    for x in leftover:
        hit, best = None, None
        for k, y in enumerate(left):
            dy = abs(y - x)
            if dy <= max(abs(x), abs(y)) * 1e-3 + 1e-9 and (best is None or dy < best):
                hit, best = k, dy
        if hit is None:
            short += 1
        else:
            del left[hit]
    if short:
        key = cfs[0] if cfs else ""
        (settled if key in SETTLED else rows).append(
            (short, len(want), ",".join(cfs), r[2][:40]))

rows.sort(reverse=True)
print("ehhnew boxes printing fewer numbers than their .ref file: %d\n" % len(rows))
print("%6s %6s  %-14s %s" % ("short", "total", "command file", "box"))
for short, tot, cf, cap in rows[:25]:
    print("%6d %6d  %-14s %s" % (short, tot, cf, cap))
if settled:
    seen = set()
    print("\nand explained, in doc/changes-from-tda.md:")
    for short, tot, cf, cap in settled:
        k = cf.split(",")[0]
        if k in seen:
            continue
        seen.add(k)
        print("  %-8s %s" % (k, SETTLED[k]))
