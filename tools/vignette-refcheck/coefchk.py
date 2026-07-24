import re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boxes import boxes
NUM = re.compile(r'-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+|-?\d+')

def ref_coefs(path):
    """Coefficient column of every estimate table."""
    lines = open(path, errors='replace').read().splitlines()
    out, intab = [], False
    for ln in lines:
        if re.search(r'Coeff\s+Error|Coeff\s+/E|Coeff\s+Coeff', ln):
            intab = True; continue
        if intab:
            if not ln.strip(): intab = False; continue
            if set(ln.strip()) <= set('-='): continue
            n = NUM.findall(ln)
            if len(n) >= 4:
                out += [float(n[-4])]   # Coeff only
    return out

def vals(txt):
    return [float(x) for x in NUM.findall(txt)]

refdir = os.environ.get('REFDIR', 'examples/exam')
B = boxes(sys.argv[1])
for spec in sys.argv[2:]:
    bns, ref = spec.split(':')
    ids = [int(i) for i in bns.split(',')]
    body = '\n'.join(B[i][1] for i in ids)
    v = vals(body)
    r = ref_coefs(os.path.join(refdir, ref))
    missed = [c for c in r if not any(abs(round(w, 4) - c) < 5e-5 for w in v)]
    hit = len(r) - len(missed)
    flag = "ok" if not missed else "  "
    print(f"{flag} {ref:<12} {bns:<8} coef    {hit}/{len(r)}" +
          (f"   missed: {missed[:6]}" if missed else ""))
