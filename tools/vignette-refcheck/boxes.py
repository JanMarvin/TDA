import re, sys, json
def boxes(md):
    txt = open(md, encoding='utf-8').read()
    parts = re.split(r'<p class="boxcap[^"]*"><span class="lab">Box (\d+)</span>', txt)
    out = {}
    for i in range(1, len(parts), 2):
        num = int(parts[i]); body = parts[i+1]
        cap = body.split('</p>')[0].strip()
        blk = body.split('</p>',1)[1]
        blk = blk.split('</div>')[0]
        blk = re.sub(r'!\[.*?\]\(.*?\)', '', blk, flags=re.S)
        out[num] = (cap, blk)
    return out
if __name__ == '__main__':
    b = boxes(sys.argv[1])
    for k in sorted(b):
        print(k, '|', b[k][0])
