#!/usr/bin/env python3
import re, sys
from pathlib import Path

SKIP_H = {"t_xwin.h", "t_xlib.h", "t_xnt.h"}
SKIP_C = {"t_xwin.c", "t_xlib.c"}
src    = Path("/tmp/tda")

def strip_comment(s):
    return re.sub(r"/\*.*?\*/", "", s).strip()

def extract_name(decl):
    m = re.search(r"([A-Za-z_][A-Za-z0-9_]*)(\[.*\])?$", decl)
    return m.group(1) if m else None

# ── collect #define constants from all sources (for resolving array sizes) ────
defines = {}
DEFINE_RE = re.compile(r"^#define\s+([A-Za-z_][A-Za-z0-9_]*)\s+(\S.*)")


for f in list(src.glob("*.h")) + list(src.glob("*.c")):
    for raw in f.read_text(errors="replace").splitlines():
        m = DEFINE_RE.match(raw.strip())
        if m and m.group(1) not in defines:
            defines[m.group(1)] = m.group(2).split()[0]

# ── collect top-level variable definitions from .c files ─────────────────────
# Must start at column 0 (no leading whitespace)
VAR_DEF = re.compile(
    r"^((?:(?:unsigned|signed|short|long|const)\s+)*"
    r"(?:int|double|float|char|short|long|FILE|size_t|void))"
    r"(\s+\*{0,3})"
    r"([A-Za-z_][A-Za-z0-9_]*)"
    r"(\[[^\]]*\](?:\[[^\]]*\])?)?"
    r"\s*(?:=|;)"
)
c_defs = {}
BLACKLIST = {"NULL","EOF","TRUE","FALSE","stdin","stdout","stderr",
             "SEEK_SET","SEEK_CUR","SEEK_END","EXIT_SUCCESS","EXIT_FAILURE"}
for c in sorted(src.glob("*.c")):
    if c.name in SKIP_C:
        continue
    depth = 0
    for raw in c.read_text(errors="replace").splitlines():
        s = strip_comment(raw).strip()
        # only collect at file scope (depth == 0)
        if depth == 0 and raw and raw[0] not in (' ', '\t', '\n'):
            if s and not s.startswith("#") and not s.startswith("//") and not s.startswith("static"):
                m = VAR_DEF.match(s)
                if m:
                    base  = m.group(1).strip()
                    stars = m.group(2).strip()
                    name  = m.group(3)
                    dims  = m.group(4) or ""
                    decl  = f"{base} {stars}{name}{dims}".strip()
                    if name not in c_defs and name not in BLACKLIST:
                        c_defs[name] = decl
                    # also extract additional names from multi-var decl "int a,b,c;"
                    rest = s[m.end():].rstrip(";").strip()
                    for chunk in rest.split(","):
                        chunk = chunk.strip().split("=")[0].strip()  # drop initializer
                        nm = re.match(r"\*{0,2}([A-Za-z_][A-Za-z0-9_]*)(?:\[[^\]]*\])?\s*$", chunk)
                        if nm and nm.group(1) and nm.group(1) not in c_defs and nm.group(1) not in BLACKLIST:
                            c_defs[nm.group(1)] = f"{base} {nm.group(1)}".strip()
        depth = max(0, depth + s.count('{') - s.count('}'))

# ── collect extern declarations from headers ──────────────────────────────────
BAD = {"XSynchronize", "XSetAfterFunction", '"C"', "Region"}
seen, decls = set(), []

for h in sorted(src.glob("*.h")):
    if h.name in SKIP_H:
        continue
    for raw in h.read_text(errors="replace").splitlines():
        line = strip_comment(raw.strip()).rstrip(";").strip()
        if not line.startswith("extern "):
            continue
        decl = line[len("extern "):].strip()
        if any(b in decl for b in BAD):
            continue
        name = extract_name(decl)
        if not name or name in seen:
            continue
        seen.add(name)
        # replace [] with concrete definition from .c files
        if "[]" in decl:
            if name in c_defs:
                decl = c_defs[name]
            else:
                # skip - no concrete definition found (shouldn't happen)
                print(f"WARNING: no concrete def for {name}", file=sys.stderr)
                continue
        decls.append(decl)

# ── file-local globals not extern'd ──────────────────────────────────────────
for name, decl in c_defs.items():
    if name not in seen:
        seen.add(name)
        decls.append(decl)

# ── rename 'ctx' to avoid clash with our parameter ────────────────────────────
decls = [
    re.sub(r'\bctx\b', 'tda_ctx_var', d) if extract_name(d) == "ctx" else d
    for d in decls
]

print(f"// {len(decls)} globals", file=sys.stderr)

# ── emit tda_context.h ────────────────────────────────────────────────────────
lines = [
    "#ifndef TDA_CONTEXT_H",
    "#define TDA_CONTEXT_H",
    "",
    "/* Zoo archive local constants needed for struct members */",
    "#ifndef SIZ_TEXT",
    "#  define SIZ_TEXT    20",
    "#  define PATHSIZE   256",
    "#  define LFNAMESIZE 256",
    "#endif",
    "",
    '#include "tda.h"',
    "#include <stdio.h>",
    "",
    "typedef struct TDAContext {",
]
for d in decls:
    lines.append(f"    {d};")
lines += [
    "} TDAContext;",
    "",
    "TDAContext *tda_context_new(void);",
    "void        tda_context_free(TDAContext *ctx);",
    "",
    "/* include after all system headers in each .c file */",
]
for d in decls:
    name = extract_name(d)
    if name:
        alias = "tda_ctx_var" if name == "ctx" else name
        lines.append(f"#define {alias} (ctx->{name})")
lines += ["", "#endif /* TDA_CONTEXT_H */", ""]

Path(src / "tda_context.h").write_text("\n".join(lines))

# ── emit tda_context.c ────────────────────────────────────────────────────────
Path(src / "tda_context.c").write_text(
    '#include <stdlib.h>\n'
    '#include "tda_context.h"\n'
    '\n'
    'TDAContext *tda_context_new(void) {\n'
    '    return calloc(1, sizeof(TDAContext));\n'
    '}\n'
    '\n'
    'void tda_context_free(TDAContext *ctx) {\n'
    '    free(ctx);\n'
    '}\n'
)
print("done", file=sys.stderr)
