#!/usr/bin/env python3
"""Freeze what the vignettes print, and check a later render against it.

    python3 tools/vignette_fingerprint.py write <html> ... > doc/...tsv
    python3 tools/vignette_fingerprint.py check doc/...tsv <html> ...

One line per box: the vignette, the box's caption, how many numbers it
prints, and a hash of those numbers in order. Captions are the key --
box numbers shift whenever a box is added, and have done repeatedly.

What is hashed is the NUMBERS, not the text, so rewording prose or
renaming a column does not trip it, while a changed value does. That is
the point: the vignettes are the largest body of checked output in the
project, and until now nothing noticed if a change somewhere else in the
package altered one of them.

A failure is not necessarily a fault -- an intended improvement changes
values too. It means: look, and re-freeze deliberately.
"""
import hashlib, html as H, re, sys

NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def boxes(path):
    doc = open(path, errors="replace").read()
    out = []
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
        cap = " ".join(H.unescape(re.sub("<[^>]+>", " ", m.group(2))).split())
        body = H.unescape(re.sub("<[^>]+>", " ", doc[st:i]))
        # the source is echoed inside the box; hash what it PRINTS, so
        # that reformatting the code does not count as a change
        body = re.sub(r"```.*?```", " ", body, flags=re.S)
        nums = NUM.findall(body)
        out.append((cap, len(nums),
                    hashlib.sha1(" ".join(nums).encode()).hexdigest()[:16]))
    return out


def main():
    mode = sys.argv[1]
    if mode == "write":
        for p in sys.argv[2:]:
            name = p.rsplit("/", 1)[-1].rsplit(".", 1)[0]
            for cap, n, h in boxes(p):
                print("%s\t%s\t%d\t%s" % (name, cap, n, h))
        return 0
    frozen = {}
    for line in open(sys.argv[2]):
        v, cap, n, h = line.rstrip("\n").split("\t")
        frozen[(v, cap)] = (int(n), h)
    seen, bad = set(), 0
    for p in sys.argv[3:]:
        name = p.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        for cap, n, h in boxes(p):
            k = (name, cap)
            seen.add(k)
            if k not in frozen:
                print("NEW      %s: %s (%d numbers)" % (name, cap[:56], n))
                bad += 1
            elif frozen[k] != (n, h):
                print("CHANGED  %s: %s" % (name, cap[:56]))
                print("         was %d numbers %s, now %d numbers %s"
                      % (frozen[k][0], frozen[k][1], n, h))
                bad += 1
    for k in sorted(frozen):
        if k not in seen:
            print("GONE     %s: %s" % (k[0], k[1][:56]))
            bad += 1
    print("\n%d box(es) differ from the frozen record (%d frozen)"
          % (bad, len(frozen)))
    return 1 if bad else 0


sys.exit(main())
