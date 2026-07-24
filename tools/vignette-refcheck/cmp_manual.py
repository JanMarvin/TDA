import re, sys
NUM = re.compile(r'(?<![\w.])-?\d+\.\d+(?![\w.])')   # decimals only: the informative numbers
man = open('/home/claude/man/manual.txt').read()
# section blocks: a heading line "6.5.3 Title" at line start (not the running header, which has extra spacing)
secs = {}
cur = None
for line in man.split('\n'):
    m = re.match(r'^(\d+(?:\.\d+){1,3})\s{1,3}([A-Z][^\d].{2,60})$', line.rstrip())
    if m and not re.search(r'\s{4,}', line.strip()):
        cur = m.group(1); secs.setdefault(cur, [])
    if cur: secs[cur].append(line)
def nums(txt):
    out=[]
    for l in txt:
        if 'tex' in l and re.search(r'\d{4}$', l.strip()): continue   # file/date footers
        out += NUM.findall(l)
    return out
vig = open('/tmp/vig/t.md').read()
vsecs = re.split(r'\n## (\d+(?:\.\d+)+) ', vig)
res=[]
for i in range(1, len(vsecs), 2):
    s = vsecs[i]; body = vsecs[i+1]
    body_nums = set(NUM.findall(body)) | set(re.findall(r'(?<![\w.])-?\d+(?![\w.])', body))
    mn = nums(secs.get(s, []))
    # dedupe, keep those with >=3 significant digits
    mn = [x for x in dict.fromkeys(mn) if len(x.replace('-','').replace('.','').lstrip('0')) >= 3 and x not in ('10.4','12.4','13.10','12.8','24.16')]
    hit = [x for x in mn if x in body_nums or any(abs(float(x)-float(v)) <= 0.51*10**-len(x.split('.')[1]) for v in body_nums)]
    res.append((s, len(hit), len(mn), [x for x in mn if x not in hit][:8]))
for s,h,n,miss in res:
    flag = 'ok ' if n and h==n else ('-- ' if n==0 else '   ')
    print(f"{flag}{s:<10} {h:>3}/{n:<3} missed: {' '.join(miss)}")
