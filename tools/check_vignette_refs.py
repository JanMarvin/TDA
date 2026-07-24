#!/usr/bin/env python3
"""Check the ehhnew vignette's numbers against TDA's own reference output.

For every row of tests/vignette-qa/vignette-qa-book.csv, take the command files it names,
read the coefficient table and the two log likelihoods out of the matching
examples/ehhnew/*.ref, and look for each of those values in the vignette box
the row points at.  Matching is numeric, at whichever of the two printed
precisions is coarser, since the vignette carries full doubles where TDA's
text is at tfmt.

    Rscript -e 'rmarkdown::render("tdaR/vignettes/ehhnew.Rmd", output_dir=...)'
    python3 tools/check_vignette_refs.py <rendered ehhnew.html>

Rows whose entry has no .ref of its own -- the gamma mixtures of chapter 10,
which ship no command file, and the pseudoresidual file of Box 8.2.3 -- are
listed as skipped; those are checked against the book by hand and the result
is recorded in the csv's note.
"""
import os, sys
import re

REF = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "examples",
    os.environ.get("TDA_EXAMPLE_DIR", "ehhnew"))

def decimals(tok):
    """Digits after the point, in the value's own scale.

    -2.5140e+03 carries four decimals in its mantissa but states the
    value to 0.05: the exponent shifts the precision. Counting after the
    "." gives 8 and demands 1e-8, which no 5-figure printout can meet.
    """
    m = re.match(r'^-?\d*(?:\.(\d*))?(?:[eE]([-+]?\d+))?$', tok)
    if not m:
        return len(tok.split(".")[1]) if "." in tok else 0
    return len(m.group(1) or "") - int(m.group(2) or 0)

def ref_numbers(cf):
    path = os.path.join(REF, cf + ".ref")
    if not os.path.exists(path):
        return None
    out = []
    gof = False
    for line in open(path, errors="replace"):
        # lsreg and glm print "Idx Wave Variable Coeff Error Coeff/E Signif",
        # a different shape from the rate family's "Idx SN Org Des MT ...";
        # both carry values the vignette has to reproduce
        m2 = re.match(r"\s*\d+\s+(?:-|\d+)\s+(\S+)\s+"
                      r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if m2:
            out.append(("coef " + m2.group(1), m2.group(2)))
            out.append(("se   " + m2.group(1), m2.group(3)))
            out.append(("t    " + m2.group(1), m2.group(4)))
            continue
        m = re.match(r"\s*\d+\s+\d+\s+\d+\s+\d+\s+\w\s+(\S+)\s+"
                     r"(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if m:
            out.append(("coef " + m.group(1), m.group(2)))
            out.append(("se   " + m.group(1), m.group(3)))
            continue
        # the global goodness-of-fit table: SN Org Des TStat DF Signif,
        # a third shape again -- neither the rate family's nor lsreg's
        if gof:
            m3 = re.match(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+"
                          r"(\d+)\s+(-?\d+\.\d+)", line)
            if m3:
                tag = "gof %s-%s" % (m3.group(2), m3.group(3))
                out.append((tag + " stat", m3.group(4)))
                out.append((tag + " sig", m3.group(6)))
                continue
            if not line.strip() or line.strip().startswith("-"):
                pass
            else:
                gof = False
        if re.match(r"\s*SN\s+Org\s+Des\s+TStat", line):
            gof = True
            continue
        m = re.match(r"\s*Log likelihood \(([^)]+)\):\s*(-?\d+\.\d+)", line)
        if m:
            out.append(("logLik " + m.group(1)[:5], m.group(2)))
    return out

def tokens(text):
    return [(t, decimals(t)) for t in re.findall(r"-?\d+(?:\.\d+)?", text)]

def found(target, toks):
    x, dx = float(target), decimals(target)
    for t, dt in toks:
        d = min(dx, dt)
        if abs(float(t) - x) <= 0.5 * 10 ** (-d) + 1e-12:
            return t
    return None


import html, csv

HTML = sys.argv[1] if len(sys.argv) > 1 else "ehhnew.html"
h = open(HTML, errors="replace").read()


def boxes_of(doc):
    """Caption and body of each numbered box.

    Splitting on the label alone runs each box into the prose that
    follows it, which then counts as that box's content -- so the body is
    taken by matching the <div class="box"> that follows the caption,
    depth-counted, and stops where the box does.
    """
    out = {}
    for m in re.finditer(r'<span class="lab">Box (\d+)</span>(.*?)</p>', doc,
                         re.S):
        n = int(m.group(1))
        start = doc.find('<div class="box">', m.end())
        if start < 0:
            continue
        depth, i = 0, start
        for t in re.finditer(r"<(/?)div\b", doc[start:]):
            depth += -1 if t.group(1) else 1
            if depth == 0:
                i = start + t.end()
                break
        out[n] = html.unescape(re.sub("<[^>]+>", " ", m.group(2) + doc[start:i]))
    return out


box = boxes_of(h)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# the book vignette's table is keyed on the book's boxes and names the
# command file in its own column; the manual's is keyed on our box and
# names the command file inside the caption
CSV = os.environ.get("TDA_QA_CSV",
                     os.path.join(ROOT, "doc", "vignette-qa-book.csv"))
BOOK = os.path.basename(CSV) == "vignette-qa-book.csv"
rows = list(csv.reader(open(CSV)))[1:]

def boxnums(spec):
    out = []
    for part in spec.split(","):
        part = part.strip()
        m = re.fullmatch(r"(\d+)-(\d+)", part)
        if m:
            out += list(range(int(m.group(1)), int(m.group(2)) + 1))
        elif part.isdigit():
            out.append(int(part))
    return out

tot = miss_tot = 0
strong = [0, 0]
report = []
for r in rows:
    if BOOK:
        bk, cfs, spec = r[0], r[3], r[4]
    else:
        bk, cfs, spec = r[3] or r[0], r[2], r[1]
    if not cfs:
        print("%-20s %-10s box %-7s no command file of its own; see the csv "
              "note" % (bk, "--", spec))
        continue
    cands = (re.findall(r"\b([a-z]{1,4}\d+[a-z]?)\.cf", cfs) if not BOOK
             else re.findall(r"\beh[a-z]\d+[a-z]?\b", cfs))
    if not cands:
        continue
    if not boxnums(spec):
        # no box names this command file: unchecked, not mismatched, and
        # counting its values as misses would make the total read as a
        # failure rate when it is a coverage gap
        print("%-20s %-10s NOT MAPPED -- no box names this command file"
              % (bk, cfs))
        continue
    text = " ".join(box.get(n, "") for n in boxnums(spec))
    toks = tokens(text)
    for cf in cands:
        nums = ref_numbers(cf)
        if nums is None:
            report.append((bk, cf, spec, "no .ref", [])); continue
        if not nums:
            report.append((bk, cf, spec, "no estimates in .ref", [])); continue
        missing = [(lab, v) for lab, v in nums if not found(v, toks)]
        tot += len(nums); miss_tot += len(missing)
        # A value of one or two significant digits matches almost any box.
        # Counted apart, so the headline does not rest on coincidence.
        for lab, v in nums:
            sig = len(str(v).replace("-", "").replace(".", "").lstrip("0"))
            if sig >= 3:
                strong[1] += 1
                if found(v, toks):
                    strong[0] += 1
        report.append((bk, cf, spec, "%d/%d" % (len(nums) - len(missing), len(nums)),
                       missing))

for bk, cf, spec, stat, missing in report:
    flag = "" if stat.count("/") and stat.split("/")[0] == stat.split("/")[1] else "  <<<"
    print("%-20s %-10s box %-7s %s%s" % (bk, cf, spec, stat, flag))
    for lab, v in missing[:8]:
        print("        missing: %-18s %s" % (lab, v))
print("\nTOTAL %d/%d reference values found in the matching vignette box"
      % (tot - miss_tot, tot))
print("      of which %d/%d carry three or more significant digits -- "
      "the ones that could not match by coincidence" % (strong[0], strong[1]))
