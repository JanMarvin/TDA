#!/usr/bin/env python3
"""Every number of the manual's output boxes the vignette does not show,
with the line it came from, classified.

    python3 tools/manual_residual.py <manual_pages.tsv> <tdaR.html>

Writes /tmp/manual-misses.tsv -- OUTSIDE the tree, because each row
quotes a line of the manual, which is not ours to distribute.

The point is to separate what the vignette genuinely does not reproduce
from what the box extractor swept in: a box runs from its caption to the
next, so surrounding prose, a figure's axis labels and the running
header all land inside it.
"""
import re, sys, os, collections

sys.argv = sys.argv[:3]
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "check_vignette_manual.py")).read()
ns = {"__name__": "x"}
# the module ends with a bare main(); drop that call, not the word
exec(src.rsplit("\nmain()", 1)[0], ns)

pages, boxes_of_section = ns["pages"], ns["boxes_of_section"]
vignette_boxes, numbers, found = (ns["vignette_boxes"], ns["numbers"],
                                  ns["found"])
kind, section_boxes = ns["kind"], ns["section_boxes"]

doc = open(sys.argv[2], errors="replace").read()
bysec, vb = section_boxes(doc), vignette_boxes(doc)

rows = []
for sec in pages:
    for num, cap, body in boxes_of_section(sec):
        if kind(cap) != "output" or not bysec.get(sec):
            continue
        have = []
        for n in sorted(bysec[sec]):
            have += numbers(vb.get(n, ""))
        for line in body.split("\n"):
            for t, d in numbers(line):
                if not found(t, d, have):
                    rows.append((sec, num, t, " ".join(line.split())[:100]))

with open("/tmp/manual-misses.tsv", "w") as f:
    for r in rows:
        f.write("\t".join(r) + "\n")


SETTING = re.compile(r"^(Maximum number of|Tolerance|Convergence criterion|"
                     r"Scaling factor|Type of covariance|Algorithm|Number of "
                     r"model parameters|Level of significance)", re.I)


# Boxes investigated and written up in doc/changes-from-tda.md. Each difference is
# explained there; they are listed so the residual report shows what is
# outstanding rather than what has already been answered.
SETTLED = {
    "6.12.5 Box 2": "manual out of date (checked against Rohwer's binary)",
    "8.4.1 Box 2": "manual out of date (checked against Rohwer's binary)",
    "6.15.2.1 Box 4": "TDA reads its data file as 4-byte floats; we use doubles",
    "6.2.2 Box 12": "gdf4.cf round-trips its data through a [10.6] file",
    "5.1.4.11 Box 1": "defective eigenvalues: repeated roots, machine-dependent",
    "5.1.4.11 Box 2": "defective eigenvalues: repeated roots, machine-dependent",
    "5.1.4.11 Box 3": "defective eigenvalues: repeated roots, machine-dependent",
    "6.14.2 Box 2": "sig is not identified on these data",
}


def classify(sec, num, tok, line):
    k = "%s Box %s" % (sec, num)
    if k in SETTLED:
        return "settled: " + SETTLED[k]
    # TDA prints its own memory use as it goes
    if re.search(r"[Cc]urrent memory|Max memory", line):
        return "a memory figure"
    # the manual's own LaTeX source stamp at the foot of each page
    if re.search(r"\.tex\s+[A-Z][a-z]+ \d+, \d{4}", line):
        return "the manual's own page footer"
    # a figure's axis labels come through as a run of round numbers
    if re.fullmatch(r"[-\d. ]+", line) and len(line.split()) >= 6 and \
       all(float(x) == int(float(x)) and int(float(x)) % 5 == 0
           for x in line.split() if re.fullmatch(r"-?\d+\.?\d*", x)):
        return "a figure's axis labels"
    # TDA echoes the settings it was given before it estimates; the
    # package's own audit counts these as settings, not results
    if SETTING.match(line.strip()):
        return "a setting echoed back, not a result"
    if re.search(r"function evaluations", line):
        return "stored in $convergence, not printed by the box"
    # the last step of a converged fit is machine-epsilon noise: the
    # manual's run and ours both stop at the tolerance, from slightly
    # different directions (2.75435e-16 against 6.88141e-16). There is
    # nothing to match.
    if re.search(r"(scaled parameter|absolute) change", line):
        try:
            if abs(float(tok)) < 1e-12:
                return "convergence noise at the tolerance"
        except ValueError:
            pass
    # a bibliographic citation: Author [year, p. NNN]
    if re.search(r"\[\d{4}[^\]]*\]", line):
        return "a bibliographic citation"
    # a run of evenly spaced round numbers is a figure's axis
    parts = line.split()
    if len(parts) >= 5 and all(re.fullmatch(r"-?\d+", x) for x in parts):
        d = [int(parts[i + 1]) - int(parts[i]) for i in range(len(parts) - 1)]
        if len(set(d)) == 1 and d[0] > 0:
            return "a figure's axis labels"
    if re.fullmatch(r"(19|20)\d\d", tok):
        return "a year in the surrounding prose"
    if re.search(r"\b(section|chapter|see|Table|Figure|Box)\b", line, re.I):
        return "a cross-reference, not output"
    if re.match(r"^\s*%s\s+[a-z]" % re.escape(sec), line):
        return "the page's running header"
    if len(tok) > 8 and "." in tok and len(tok.split(".")[0]) >= 5:
        return "two columns run together in the scan"
    if re.search(r"\bc\d{1,2}\b", line):
        return "a command file's column reference"
    # a box that quotes TDA's own command syntax -- mdef(T,3,3) = ... --
    # is showing the input, not a result
    if re.search(r"\b[a-z][a-z0-9_]{2,}\s*\([^)]*\)\s*=", line):
        return "command syntax quoted in the box"
    if len(re.findall(r"[A-Za-z]{3,}", line)) > 5:
        return "a figure in a sentence"
    return "UNCLASSIFIED"


c = collections.Counter(classify(*r) for r in rows)
tot = sum(c.values())
print("numbers of the manual's output boxes not shown by the vignette: %d\n"
      % tot)
for k, v in c.most_common():
    print("  %5d  %4.1f%%  %s" % (v, 100 * v / tot, k))
un = [r for r in rows if classify(*r) == "UNCLASSIFIED"]
print("\nunclassified, by section:")
for s, n in collections.Counter("%s Box %s" % (r[0], r[1]) for r in un).most_common(12):
    print("  %-22s %d" % (s, n))
