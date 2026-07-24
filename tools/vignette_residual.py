"""Group the book-comparison misses by box, for working through by hand.

Writes /tmp/vignette-misses.tsv, OUTSIDE the repository on purpose: each
row carries the line of Blossfeld and Rohwer the number came from, and
that text is not ours to distribute.  Nothing derived from the scan is
committed or shipped.
"""
import re, csv, html, sys, collections
sys.path.insert(0, 'TDA/tools')
src = open('TDA/tools/check_vignette_book.py').read()
ns = {}
exec(src.split('bb = book_boxes')[0].replace('BOOK = sys.argv[1]', 'BOOK=None')
        .replace('HTML = sys.argv[2]', 'HTML=None')
        .replace('CSV = sys.argv[3]', 'CSV=None'), ns)
book_boxes, boxes_of = ns['book_boxes'], ns['boxes_of']
numbers, found, decimals = ns['numbers'], ns['found'], ns['decimals']
body_lines = ns['body_lines']

bb_raw = {}
txt = open('book.txt', errors='replace').read().split('\n')
CAP, SKIP = ns['CAP'], ns['SKIP']
starts = []
for i, l in enumerate(txt):
    m = CAP.match(l)
    if m and not SKIP.match(m.group(3).strip()):
        starts.append((i, ("Box " if m.group(1).lower() == "box" else "Fig. ") + m.group(2)))
for k, (i, name) in enumerate(starts):
    j = starts[k + 1][0] if k + 1 < len(starts) else len(txt)
    bb_raw.setdefault(name, body_lines(txt[i + 1:j]))

vb = boxes_of(open('vout/ehhnew.html', errors='replace').read())
rows = list(csv.reader(open('TDA/tests/vignette-qa/vignette-qa-book.csv')))[1:]

out = []
for r in rows:
    name = r[0] if r[0].startswith("Fig.") else "Box " + r[0]
    if name not in bb_raw or not r[4]:
        continue
    have = []
    for n in re.findall(r"\d+", r[4]):
        have += numbers(vb.get(int(n), ""))
    for line in bb_raw[name]:
        for t, d in numbers(line, drop_refs=True):
            if not found(t, d, have):
                out.append((name, t, line.strip()))
print("total misses:", len(out))
with open('/tmp/vignette-misses.tsv', 'w') as f:
    for name, t, line in out:
        f.write("%s\t%s\t%s\n" % (name, t, line[:120]))
c = collections.Counter(n for n, _, _ in out)
for n, k in c.most_common(18):
    print("  %-12s %d" % (n, k))
