import re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boxes import boxes

NUM = re.compile(r'-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+|-?\d+')
REFDIR = os.environ.get('REFDIR', 'examples/exam')

# Every line of a ref that looks like a result row: two or more numbers,
# and not one of the command-file echo's own rows (a variable definition,
# a memory figure, a timestamp).
DEF = re.compile(r'\bc\d+\b|(ge|le|eq|lt|gt|if)\(|PFmt')
ROW = re.compile(r'^[ \t]*[-0-9\[]')
LOGL = re.compile(r'(log likelihood|Log likelihood|Minimum of function|'
                  r'Maximum of log likelihood|Chi2|chi2)', re.I)

def ref_results(path):
    lines = open(path, errors='replace').read().splitlines()
    keep = []
    SKIP = re.compile(r'memory|TDA\. Analysis|Read records|Reading|'
                      r'Maximum number of|Number of (cases|variables)|'
                      r'Current|bytes|End of|Iter|Function Value|'
                      r'iterations|evaluations|Tolerance|Algorithm|'
                      r'criterion|Armijo|step size|Scaling|version')
    for ln in lines:
        if not ln.strip() or SKIP.search(ln) or DEF.search(ln):
            continue
        if re.search(r'\(\d+,\d+\)', ln):   # the minimiser's own protocol
            continue
        if LOGL.search(ln):
            keep.append(ln); continue
        if not ROW.match(ln):
            continue
        if len(NUM.findall(ln)) >= 2:
            keep.append(ln)
    return numbers('\n'.join(keep))

def numbers(txt):
    out = []
    for m in NUM.finditer(txt):
        s = m.group(0)
        try: v = float(s)
        except ValueError: continue
        d = len(s.split('.')[1]) if '.' in s and 'e' not in s.lower() else None
        out.append((v, d, s))
    return out

def compare(refnums, vnums):
    vals = [v for v, _, _ in vnums]
    hit, missed = 0, []
    seen = set()
    for v, d, s in refnums:
        if s in seen: continue
        seen.add(s)
        ok = False
        for w in vals:
            if d is not None:
                if abs(round(w, d) - v) < 10**(-d)/2 + 1e-9: ok = True; break
            elif v and abs(w - v) <= abs(v)*1e-3: ok = True; break
            elif v == 0 and w == 0: ok = True; break
        if ok: hit += 1
        else: missed.append(s)
    return hit, len(seen), missed

md = sys.argv[1]
B = boxes(md)
for spec in sys.argv[2:]:
    bns, ref = spec.split(':')
    ids = [int(i) for i in bns.split(',')]
    cap = B[ids[0]][0]
    body = '\n'.join(B[i][1] for i in ids)
    r = ref_results(os.path.join(REFDIR, ref))
    hit, tot, missed = compare(r, numbers(body))
    flag = "ok" if hit == tot else "  "
    label = ','.join(str(i) for i in ids)
    print(f"{flag} Box {label:<10} {cap[:38]:<38} {ref:<13} {hit}/{tot}" +
          (f"   missed: {', '.join(missed[:8])}" if missed else ""))
