#!/usr/bin/env python3
"""Check the ehhnew vignette against the numbers the BOOK prints.

    python3 tools/check_vignette_book.py book.txt ehhnew.html \
            tests/vignette-qa/vignette-qa-book.csv

where book.txt is `pdftotext -layout` over a scan of Blossfeld and Rohwer.
For every box and figure, take the numbers printed there and look for each
in the vignette box tests/vignette-qa/vignette-qa-book.csv maps it to.

READ THE OUTPUT WITH CARE.  This is a lead-finder, not a gate; the
authoritative comparison is tools/check_vignette_refs.py, which runs
against TDA's own reference output.  A miss here is usually one of:

  - an OCR error in the scan.  Box 4.3.2b prints the constant of ehd7 as
    -6.0179 with C/Error -12.0226; TDA prints -5.0179 / -12.0225, and
    -6.0179/0.4174 would be -14.4, so the book's own two columns
    disagree and the 5 was read as a 6.  The scan also turns 1 into I
    and 0 into O.
  - a number in the command-file listing with no counterpart in the R:
    column indices (c6, c12), the century-month constants 468/504/588/624
    behind COHO2 and COHO3, which tda_rrdat() derives for you.
  - a row of a long table that the book prints and the vignette shows
    only the head of.  These are the ones worth chasing.
  - a page folio or a year that leaked out of the surrounding prose.
"""
import re, csv, html, sys, collections

BOOK = sys.argv[1]
HTML = sys.argv[2]
CSV = sys.argv[3]

CAP = re.compile(r'\s*([Bb][Oo][Xx]|Figure)\s+(\d+\.\d+\.\d+[a-z]?)\s+(\S.*)')
SKIP = re.compile(r"^(shows|show|illustrates|displays|gives|lists|is |are |in |"
                  r"of |and |for |above|below|presents|contains)", re.I)


def book_boxes(path):
    txt = open(path, errors="replace").read().split("\n")
    starts = []
    for i, l in enumerate(txt):
        m = CAP.match(l)
        if m and not SKIP.match(m.group(3).strip()):
            starts.append((i, ("Box " if m.group(1).lower() == "box"
                               else "Fig. ") + m.group(2)))
    out = {}
    for k, (i, name) in enumerate(starts):
        j = starts[k + 1][0] if k + 1 < len(starts) else len(txt)
        out.setdefault(name, "\n".join(body_lines(txt[i + 1:j])))
    return out


def body_lines(lines):
    """The box's own content, without the prose that follows it.

    A box in the book is a framed listing or table; the scan gives no
    frame, only the text.  A table row or a command-file line carries
    numbers and few words; the prose around it is the other way round.
    Page headers (a folio and a running title) go the same way.
    """
    keep = []
    for l in lines:
        t = l.strip()
        if not t:
            continue
        nums = re.findall(r"-?\d+(?:\.\d+)?", t)
        words = re.findall(r"[A-Za-z]{2,}", t)
        if len(words) > 8 and len(nums) < 3:
            continue                      # a sentence
        if len(nums) == 1 and len(words) > 3:
            continue                      # prose with one figure in it
        if re.search(r"\b(section|chapter|see|Table|Figure|cf\.)\b", t, re.I):
            continue                      # a cross-reference, not content
        if "%" in t and len(words) > 2:
            continue                      # a sentence quoting a percentage
        keep.append(t)
    return keep


def boxes_of(doc):
    out = {}
    for m in re.finditer(r'<span class="lab">Box (\d+)</span>(.*?)</p>', doc, re.S):
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


def decimals(t):
    return len(t.split(".")[1]) if "." in t else 0


def numbers(text, drop_refs=False):
    if drop_refs:
        # "(Box 4.1.1)", "Box 7.2 3", section numbers in the caption line
        text = re.sub(r"\(?\b[Bb]ox[^)\n]{0,18}\)?", " ", text)
        text = re.sub(r"\bFigure\s+[\d. ]+", " ", text)
    return [(t, decimals(t)) for t in re.findall(r"-?\d+(?:\.\d+)?", text)]


def found(target, dt, toks):
    x = float(target)
    for t, d in toks:
        k = min(dt, d)
        if abs(float(t) - x) <= 0.5 * 10 ** (-k) + 1e-9:
            return True
    return False


bb = book_boxes(BOOK)
vb = boxes_of(open(HTML, errors="replace").read())
rows = list(csv.reader(open(CSV)))[1:]

tot = hit = 0
report = []
for r in rows:
    name = r[0] if r[0].startswith("Fig.") else "Box " + r[0]
    spec = r[4]
    if name not in bb:
        report.append((name, "not found in the scan", 0, 0, []))
        continue
    # A "Command file X.cf" box is TDA's own command syntax: column
    # references (c3, c7), the century-month constants behind COHO2 and
    # COHO3, model numbers.  The R does the same work through different
    # means and reproduces none of those tokens, so they are not content
    # to compare -- the result box that follows carries what has to match.
    if re.match(r"(Part of )?[Cc]ommand file", r[2].strip()):
        report.append((name, "command syntax, not output", 0, 0, []))
        continue
    if name.startswith("Fig."):
        # a figure prints no numbers of its own; only its caption and the
        # prose around it, which is not content to reproduce.  Figures are
        # checked by eye against the page.
        report.append((name, "graph, checked against the page by eye", 0, 0, []))
        continue
    want = numbers(bb[name], drop_refs=True)
    if not spec:
        report.append((name, "no vignette box", 0, len(want), []))
        continue
    have = []
    for n in re.findall(r"\d+", spec):
        have += numbers(vb.get(int(n), ""))
    miss = [t for t, d in want if not found(t, d, have)]
    tot += len(want)
    hit += len(want) - len(miss)
    report.append((name, "", len(want) - len(miss), len(want), miss))

for name, note, h, w, miss in report:
    flag = "" if (h == w) else "   <<<"
    print("%-12s %s%s%s" % (name, note or "%d/%d" % (h, w), flag,
                            ("  missing: " + " ".join(miss[:10])) if miss[:10] else ""))
print("\nTOTAL %d/%d numbers printed in the book found in the vignette box "
      "that reproduces it" % (hit, tot))
