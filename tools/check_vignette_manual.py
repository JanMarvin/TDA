#!/usr/bin/env python3
"""Check the tdaR vignette against the numbers the MANUAL prints.

    python3 tools/check_vignette_manual.py <manual_pages.tsv> <tdaR.html>

tman1.pdf puts two logical pages on each landscape sheet; qa/manual_extract.py
splits them at x = 420pt and normalises TeX's ligatures, writing one row per
logical page. This reads that, finds each box, and looks for the numbers it
prints in the vignette box that reproduces it.

Boxes are keyed by (section, number) because the manual numbers them within
a section. Only OUTPUT boxes are compared: a "Syntax for ..." box is
reference material, a "Command file x.cf" box is TDA's own command syntax,
and a "Data file" box is a listing of the shipped data -- the R reaches all
three by other means and reproduces none of their tokens.

Matching is numeric at whichever of the two printed precisions is coarser.
"""
import re, sys, html, collections

PAGES, HTML = sys.argv[1], sys.argv[2]

def kind(cap):
    c = cap.lower()
    if re.match(r'^(syntax|options|parameters?|types? of|storage|default|'
                r'general commands|commands for|list of|overview|format|'
                r'notation|summary of|variables? in|models available|'
                r'built-in|structure of|algorithm for)', c):
        return "syntax"
    if "command file" in c:
        return "command file"
    if re.search(r'\bdata files?\b|\.dat\b|\.df\b|example data', c):
        return "data file"
    return "output"

pages = collections.defaultdict(list)
for line in open(PAGES):
    s, p, t = line.rstrip("\n").split("\t", 2)
    if s:
        pages[s].append((int(p or 0), t.replace("\\n", "\n")))

# The running header of every logical page -- "5.1.4.11   eigenvalues and
# eigenvectors   6" -- sits inside any box that spans a page break, and
# its section number and page number are not content.  Likewise a print
# format ("mfmt = 13.10;") is a command, not a value.
HDR_LINE = re.compile(r'^\s*[\d.]+\s+[a-z][a-z ,\-]+\s+\d+\s*$')
FMT_LINE = re.compile(r'^\s*\w*fmt\s*=')


def strip_furniture(text):
    return "\n".join(l for l in text.split("\n")
                     if not HDR_LINE.match(l) and not FMT_LINE.match(l))


def boxes_of_section(sec):
    text = strip_furniture("\n".join(t for _, t in sorted(pages[sec])))
    out = []
    for m in re.finditer(r'^\s*Box\s+(\d+)\s+(\S.{3,90})$', text, re.M):
        rest = text[m.end():]
        stop = re.search(r'^\s*(Box\s+\d+\s+\S|\d+\.\d[\d.]*\s+[A-Z])', rest, re.M)
        out.append((m.group(1), " ".join(m.group(2).split()),
                    rest[:stop.start()] if stop else rest))
    return out

def vignette_boxes(doc):
    out = {}
    for m in re.finditer(r'<span class="lab">Box (\d+)</span>(.*?)</p>', doc, re.S):
        n = int(m.group(1))
        st = doc.find('<div class="box">', m.end())
        if st < 0:
            continue
        d, i = 0, st
        for t in re.finditer(r"<(/?)div\b", doc[st:]):
            d += -1 if t.group(1) else 1
            if d == 0:
                i = st + t.end()
                break
        out[n] = html.unescape(re.sub("<[^>]+>", " ", m.group(2) + doc[st:i]))
    return out


def decimals(t):
    """Digits after the point, in the value's own scale.

    A number written as -2.5140e+03 carries four decimals in its
    mantissa but states the value to 0.05: the exponent shifts the
    precision. Counting the characters after the "." gives 8 and demands
    agreement to 1e-8, which no 5-significant-figure printout can meet.
    Every scientific-notation value in the manual was being compared
    that way.
    """
    m = re.match(r'^-?\d*(?:\.(\d*))?(?:[eE]([-+]?\d+))?$', t)
    if not m:
        return len(t.split(".")[1]) if "." in t else 0
    frac = len(m.group(1) or "")
    exp = int(m.group(2) or 0)
    return frac - exp


def numbers(text):
    return [(t, decimals(t)) for t in
            re.findall(r"-?\d+(?:\.\d+)?(?:e[-+]?\d+)?", text)]


def found(tok, d, toks):
    x = float(tok)
    for t, dt in toks:
        k = min(d, dt)
        if abs(float(t) - x) <= 0.5 * 10 ** (-k) + 1e-9:
            return True
    return False


def section_boxes(doc):
    """Which vignette boxes sit under which manual section.

    Taken from the vignette's own headings rather than from the QA table:
    a box belongs to the last "### N.N.N" heading before it. The table's
    box numbers go stale every time a box is added, the headings do not,
    and a section split across several boxes needs no special case.
    """
    marks = []
    for m in re.finditer(r'<h[23][^>]*>\s*([\d.]+)\s', doc):
        marks.append((m.start(), m.group(1)))
    for m in re.finditer(r'<span class="lab">Box (\d+)</span>', doc):
        marks.append((m.start(), int(m.group(1))))
    marks.sort()
    out, cur = collections.defaultdict(list), None
    for _, v in marks:
        if isinstance(v, str):
            cur = v
        elif cur:
            out[cur].append(v)
    return out


def main():
    doc = open(HTML, errors="replace").read()
    bysec = section_boxes(doc)
    vb = vignette_boxes(doc)

    tot = hit = 0
    strong = [0, 0]
    report = []
    for sec in sorted(pages, key=lambda s: [int(x) for x in s.split(".")]):
        for num, cap, body in boxes_of_section(sec):
            k = kind(cap)
            if k != "output":
                report.append((sec, num, k, 0, 0, cap, []))
                continue
            boxes = sorted(bysec.get(sec, ()))
            if not boxes:
                report.append((sec, num, "NOT MAPPED", 0, 0, cap, []))
                continue
            have = []
            for n in boxes:
                have += numbers(vb.get(n, ""))
            want = numbers(body)
            miss = [t for t, d in want if not found(t, d, have)]
            tot += len(want)
            hit += len(want) - len(miss)
            # A token of one or two significant digits matches almost
            # anything: a box with a handful of integers in it contains a
            # "1" and a "2" whatever it is of.  Counted separately, so the
            # headline cannot rest on them.
            for t, d in want:
                sig = len(t.replace("-", "").replace(".", "").lstrip("0"))
                if sig >= 3:
                    strong[1] += 1
                    if found(t, d, have):
                        strong[0] += 1
            report.append((sec, num, ",".join(map(str, boxes)),
                           len(want) - len(miss), len(want), cap, miss))
    for sec, num, where, h, w, cap, miss in report:
        if where in ("syntax", "command file", "data file"):
            continue
        flag = "" if (w and h == w) else "   <<<"
        s = "%d/%d" % (h, w) if w else where
        print("%-10s Box %-3s %-10s %-9s %s%s"
              % (sec, num, where if w else "", s, cap[:40], flag))
        for t in miss[:6]:
            print("        missing: %s" % t)
    # The denominator counts only sections the vignette has a box under.
    # Report the rest too: a section with no box is not compared at all,
    # and a total that silently leaves it out looks better than it is.
    skipped = 0
    for sec in pages:
        if bysec.get(sec):
            continue
        for num, cap, body in boxes_of_section(sec):
            if kind(cap) == "output":
                skipped += len(numbers(body))
    print("\nTOTAL %d/%d numbers printed in the manual's output boxes found "
          "in the vignette" % (hit, tot))
    print("      of which %d/%d carry three or more significant digits -- "
          "the ones that could not match by coincidence"
          % (strong[0], strong[1]))
    print("      %d more in %d section(s) the vignette has no box for; "
          "%d in the manual's output boxes altogether"
          % (skipped, sum(1 for sec in pages if not bysec.get(sec)
                          and any(kind(c) == "output"
                                  for _, c, _ in boxes_of_section(sec))),
             tot + skipped))


main()
