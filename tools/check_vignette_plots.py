#!/usr/bin/env python3
"""Compare every plot in the ehhnew vignette with its command file.

    python3 tools/check_vignette_plots.py <rendered ehhnew.html>

Each plotting .cf in examples/ehhnew declares its coordinate system
(pxa/pya), its tick spacing (plxa/plya sc= and ic=), how many series it
draws (one plot= per series) and its labels (pltext, with the point it
is drawn at).  The vignette's own tda_ps() calls declare the same
things.  Anything that disagrees is reported: a picture on the wrong
axes, a missing series, a label the book has and the vignette does not.
"""
import re, sys, os, html, glob

# which example directory the vignette's command files live in: the book
# vignette uses examples/ehhnew, the manual's uses examples/exam
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "examples",
                  os.environ.get("TDA_EXAMPLE_DIR", "ehhnew"))


def cf_spec(path):
    t = open(path, errors="replace").read()
    t = re.sub(r"#[^\n]*", "", t)          # comments, incl. #plframe
    g = lambda p: (re.search(p, t) or [None, None])
    spec = {}
    m = re.search(r"pxa\s*=\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)", t)
    if m:
        spec["xlim"] = (float(m.group(1)), float(m.group(2)))
    m = re.search(r"pya\s*=\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)", t)
    if m:
        spec["ylim"] = (float(m.group(1)), float(m.group(2)))
    for ax in ("plxa", "plya"):
        m = re.search(ax + r"\s*\(([^)]*)\)", t)
        if m:
            sc = re.search(r"sc\s*=\s*([\d.]+)", m.group(1))
            ic = re.search(r"ic\s*=\s*(\d+)", m.group(1))
            spec[ax] = (float(sc.group(1)) if sc else None,
                        int(ic.group(1)) if ic else None)
    spec["series"] = len(re.findall(r"^\s*plot\s*[(=]", t, re.M))
    spec["text"] = [(s.strip().strip('"').rstrip(";").strip(), float(x), float(y))
                    for x, y, s in
                    re.findall(r'pltext\s*\(\s*xy\s*=\s*(-?[\d.]+)\s*,\s*'
                               r'(-?[\d.]+)\s*\)\s*=\s*([^;]+);', t)]
    spec["frame"] = bool(re.search(r"^\s*plframe", t, re.M))
    return spec


def vignette_plots(doc):
    """Each box's tda_ps()/tda_pl() source, keyed by box number."""
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
        # keep the caption apart from the body: the tags are stripped
        # below, so a later split on the box <div> would never match and
        # the whole box would be read as its own caption
        cap = html.unescape(re.sub("<[^>]+>", " ", m.group(2)))
        body = html.unescape(re.sub("<[^>]+>", " ", doc[start:i]))
        if "tda_ps(" in body:
            out[n] = (cap, body)
    return out


def vig_spec(src):
    spec = {}
    m = re.search(r"xlim\s*=\s*c\(\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)", src)
    if m:
        spec["xlim"] = (float(m.group(1)), float(m.group(2)))
    m = re.search(r"ylim\s*=\s*c\(\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)", src)
    if m:
        spec["ylim"] = (float(m.group(1)), float(m.group(2)))
    for ax, key in (("plxa", "plxa"), ("plya", "plya")):
        # ic= is optional in both the command files and the R
        m = re.search(r'"%s"\s*,\s*sc\s*=\s*([\d.]+)(?:\s*,\s*ic\s*=\s*(\d+))?'
                      % ax, src)
        if m:
            spec[key] = (float(m.group(1)),
                         int(m.group(2)) if m.group(2) else None)
    # by= draws one series per level in a single call, and band= adds two
    # more plot() commands per series (the grey upper and white lower TDA
    # shades between); neither is countable from the source, so a call
    # carrying either is not compared on series count
    spec["series"] = len(re.findall(r"tda_pl_lines\(|tda_pl_points\(", src))
    spec["countable"] = not re.search(r"by\s*=|band\s*=", src)
    spec["text"] = [(s.strip().strip('"').rstrip(";").strip(), float(x), float(y))
                    for s, x, y in
                    re.findall(r'tda_pl_text\([^,]+,\s*("(?:[^"\\]|\\.)*")\s*,\s*'
                               r'at\s*=\s*c\(\s*(-?[\d.e]+)\s*,\s*(-?[\d.e]+)', src)]
    # tda_pl_axes(p, sc=, ic=) sets both axes in one call, and
    # tda_pl_frame() is the frame; the gallery code uses these forms
    # sc= and ic= take either one value for both axes or c(x, y)
    def pair(arg):
        m = re.search(r'tda_pl_axes\([^\n]*?%s\s*=\s*c\(\s*([\d.]+)\s*,'
                      r'\s*([\d.]+)' % arg, src)
        if m:
            return float(m.group(1)), float(m.group(2))
        m = re.search(r'tda_pl_axes\([^\n]*?%s\s*=\s*([\d.]+)' % arg, src)
        return (float(m.group(1)),) * 2 if m else None
    sc, ic = pair("sc"), pair("ic")
    if sc:
        spec.setdefault("plxa", (sc[0], int(ic[0]) if ic else None))
        spec.setdefault("plya", (sc[1], int(ic[1]) if ic else None))
    spec["frame"] = '"plframe"' in src or "tda_pl_frame(" in src
    return spec


def exceptions():
    """Labels the vignette adds on purpose, and why (doc/...-exceptions.tsv)."""
    path = os.path.join(ROOT, "doc", "vignette-plot-exceptions.tsv")
    out = {}
    if os.path.exists(path):
        for line in open(path, errors="replace"):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out.setdefault(f[0], set()).add(f[1])
    return out


def main():
    allowed = exceptions()
    doc = open(sys.argv[1], errors="replace").read()
    plots = vignette_plots(doc)
    bad = 0
    for n in sorted(plots):
        src = plots[n]
        # in the order the caption names them, which is the order the box
        # draws them -- sorting puts ehc10 before ehc6 and pairs the plots
        # with the wrong command files
        seen, cfs = set(), []
        # only what the CAPTION names.  Reading the code body instead
        # pairs a data frame called dh1 with dh1.cf, which is how the
        # kernel-density box came to be compared against a histogram.
        caption, src = src
        # Both caption styles count: ehhnew writes the bare stem
        # ("ehd2p"), the manual's vignette writes it out ("qr1.cf").  What
        # does NOT count is a name written with another extension --
        # "gd1.dat" is the data the plot reads, not the command file it
        # reproduces, and matching the bare token there pairs the plot
        # with the wrong file.
        toks = re.findall(r"\b([a-z]{1,5}\d+[a-z]?)(\.[a-z]+)?\b", caption)
        for c, ext in toks:
            if ext and ext != ".cf":
                continue
            if c not in seen:
                seen.add(c)
                cfs.append(c)
        cands = [c for c in cfs if os.path.exists(os.path.join(EX, c + ".cf"))
                 and "psfile" in open(os.path.join(EX, c + ".cf"),
                                      errors="replace").read()]
        if not cands:
            continue
        # A box can hold several plots.  Compared as one specification the
        # first xlim wins and the labels of every picture merge, which hides
        # exactly the kind of fault this is looking for -- so the source is
        # cut at each tda_ps() call and the pieces paired with the command
        # files in order.
        pieces = re.split(r"(?=tda_ps\()", src)[1:]
        if len(pieces) != len(cands):
            bad += 1
            print("Box %-3d %-8s  NOT COMPARED: %d tda_ps() call(s) against %d "
                  "command file(s) -- split the box, one plot per .cf"
                  % (n, ",".join(cands)[:8], len(pieces), len(cands)))
            continue
        for c, piece in zip(cands, pieces):
            v = vig_spec(piece)
            src_piece = piece
            k = cf_spec(os.path.join(EX, c + ".cf"))
            diffs = []
            for key in ("xlim", "ylim", "plxa", "plya", "frame"):
                if key not in k:
                    continue
                if key not in v:
                    diffs.append("%s UNREADABLE in the vignette source "
                                 "(cf has %s)" % (key, k[key]))
                elif key == "frame" and "FRAME" in allowed.get(c, ()):
                    pass
                elif (key in ("plxa", "plya") and k[key][1] is None
                      and k[key][0] == v[key][0]):
                    pass
                elif k[key] != v[key]:
                    diffs.append("%s %s vs cf %s" % (key, v[key], k[key]))
            if (v.get("countable", True) and k["series"] and
                    k["series"] != v["series"]):
                diffs.append("series %d vs cf %d" % (v["series"], k["series"]))
            # a label passed as a variable rather than a literal cannot be
            # read out of the source; those calls are left uncompared
            kt = sorted(t[0] for t in k["text"])
            vt = sorted(t[0] for t in v["text"])
            literal = not re.search(r"tda_pl_text\([^,]+,\s*[a-zA-Z.]", src_piece)
            vt = [t for t in vt if t not in allowed.get(c, ())]
            if literal and kt != vt:
                diffs.append("labels %s vs cf %s" % (vt, kt))
            if diffs:
                bad += 1
                print("Box %-3d %-8s  %s" % (n, c, "; ".join(diffs)))
    print("\n%d plot(s) disagree with their command file" % bad)


main()
