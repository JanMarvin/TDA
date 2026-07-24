#!/bin/bash
# Rebuild from .orig files (requires .orig files alongside sources)
set -e
cd "$(dirname "$0")/../src"

# restore originals
for f in *.c.orig *.h.orig; do
    cp "$f" "${f%.orig}"
done

# patch tda.h
python3 - << 'EOF'
from pathlib import Path; import re
h = Path("tda.h").read_text()
lines = h.splitlines(keepends=True)
last_sys = max((i for i,l in enumerate(lines) if l.strip().startswith("#include <")), default=-1)
lines.insert(last_sys+1,'#include "tda_const.h"\n')
lines.insert(last_sys+2,'#include "tda_context.h"\n')
Path("tda.h").write_text("".join(lines))
EOF

python3 ../tools/make_context.py

python3 - << 'EOF'
from pathlib import Path; import re
h = Path("tda_context.h").read_text()
h = h.replace('#include "tda.h"\n','#include "tda_const.h"\n',1)
h = re.sub(r"/\* Constants defined.*?#endif\n\n","",h,flags=re.DOTALL)
for name in ["BHHMPC","BHHGPC","S1","S2","S3","S4","S5","S6","S7","DIG30"]:
    h = h.replace(f"    char {name}[];",f"    char *{name};")
h = re.sub(r"\n/\* include after.*","",h,flags=re.DOTALL)
h = h.rstrip()+"\n\n#endif /* TDA_CONTEXT_H */\n"
Path("tda_context.h").write_text(h)
EOF

python3 ../tools/thread_ctx.py
python3 ../tools/post_fixes.py
