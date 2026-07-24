#!/usr/bin/env python3
"""
Full context-threading pass on TDA sources. Run once on unmodified files.

For each .c file:
  1. Add #include "tda_context.h" after last system include
  2. Add TDAContext *ctx as first param to all function definitions
  3. Add ctx as first arg to all call sites (non-declaration lines only)
  4. Remove file-scope global definitions now owned by TDAContext
  5. Replace bare global names with ctx->name inside function bodies

For each .h file:
  6. Add TDAContext *ctx to function declarations
  7. Add #include "tda_ctx_fwd.h" where needed
  8. Strip extern declarations for globals now in struct
"""

import re, sys, shutil
from pathlib import Path

SKIP_C = {"t_xwin.c", "t_xlib.c", "tda_context.c"}
SKIP_H = {"t_xwin.h", "t_xlib.h", "t_xnt.h",
          "tda_context.h", "tda_const.h", "tda_ctx_fwd.h"}
src = Path("/tmp/tda")

def backup(p):
    bak = p.with_suffix(p.suffix + ".orig")
    if not bak.exists():
        shutil.copy2(p, bak)

# ── load globals from struct ─────────────────────────────────────────────────
ctx_h = (src / "tda_context.h").read_text()
BLACKLIST = {"NULL","EOF","TRUE","FALSE","stdin","stdout","stderr",
            "SEEK_SET","SEEK_CUR","SEEK_END","EXIT_SUCCESS","EXIT_FAILURE",
            "ctx"}
# extract ALL names from struct member lines incl. multi-name lines
# collect ALL #define constants from all headers and .c files
DEFINED_CONSTS = set()
for _hf in list(src.glob('*.h')) + list(src.glob('*.c')):
    for _dl in _hf.read_text(errors='replace').splitlines():
        _dm = re.match(r'^#define\s+([A-Za-z_][A-Za-z0-9_]*)', _dl)
        if _dm: DEFINED_CONSTS.add(_dm.group(1))
global_names = set()
_TW = {"int","double","float","char","short","long","void","FILE",
       "size_t","unsigned","signed","const","TDAContext","struct"}
for _mm in re.finditer(r"^\s+(\S[^{};]*);", ctx_h, re.MULTILINE):
    for _nm in re.findall(r"\b([A-Za-z_][A-Za-z0-9_]+)\b", _mm.group(1)):
        if _nm not in _TW and _nm not in DEFINED_CONSTS:
            global_names.add(_nm)
sorted_globals = sorted(global_names, key=len, reverse=True)
GLOBAL_RE = re.compile(r'\b(' + '|'.join(re.escape(n) for n in sorted_globals) + r')\b')
print(f"{len(global_names)} globals", file=sys.stderr)

# ── Pass 1: function names ───────────────────────────────────────────────────
FUNC_DEF_RE = re.compile(
    r"^((?:unsigned\s+|signed\s+|short\s+|long\s+|struct\s+[A-Za-z_][A-Za-z0-9_]*\s+)?[A-Za-z_][A-Za-z0-9_ *]*?)\s+(\*?)([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)\s*$"
)
func_names = set()
for c in sorted(src.glob("*.c")):
    if c.name in SKIP_C:
        continue
    lines = c.read_text(errors="replace").splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i]
        if not raw or raw[0] in (' ', '\t', '#', '/', '*'):
            i += 1; continue
        if ';' in raw:
            i += 1; continue
        # join continuation lines (sig spans multiple lines when ) is not on first line)
        joined = raw
        j = i + 1
        while ')' not in joined and j < len(lines) and j < i + 8:
            joined += " " + lines[j].strip()
            j += 1
        joined_nc = re.sub(r'\s*/\*.*?\*/\s*$', '', joined.strip())
        joined_nc = re.sub(r'\s*//.*$', '', joined_nc)
        m = FUNC_DEF_RE.match(joined_nc)
        if m:
            # check next non-blank line for {
            for k in range(j, min(j+3, len(lines))):
                nxt = lines[k].strip()
                if nxt:
                    if nxt.startswith('{'):
                        func_names.add(m.group(3))
                    break
        i += 1

print(f"{len(func_names)} functions", file=sys.stderr)

CALL_RE = re.compile(
    r'\b(' + '|'.join(re.escape(n) for n in sorted(func_names, key=len, reverse=True)) + r')\s*\('
)

def add_ctx(params):
    p = params.strip()
    return "TDAContext *ctx" if p in ("", "void") else "TDAContext *ctx, " + p

# ── file-scope global definition pattern ────────────────────────────────────
# matches: TYPE [*] NAME [= ...] ; at column 0
GDEF_RE = re.compile(
    r"^(?:(?:unsigned|signed|short|long|const)\s+)*"
    r"(?:int|double|float|char|short|long|FILE|size_t|void)\s+"
    r"\*{0,3}([A-Za-z_][A-Za-z0-9_]*)"
    r"(?:\[[^\]]*\])*"
    r"(?:\s*=\s*[^;]*)?\s*[;{]"   # also matches initialiser block start
)
# also match "unsigned\nchar NAME" split across two lines - handle via lookahead in step B

# protect strings/comments/preprocessor for substitution
PROTECT_RE = re.compile(
    r'("(?:[^"\\]|\\.)*")'
    r"|('(?:[^'\\]|\\.)*')"
    r"|(//[^\n]*)"
    r"|(/\*.*?\*/)"
    r"|(#[^\n]*)",
    re.DOTALL
)

LOCAL_DECL_RE = re.compile(
    r"^\s*(?:register\s+|volatile\s+|const\s+)*"
    r"(?:unsigned\s+|signed\s+|short\s+|long\s+)*"
    r"(?:int|double|float|char|short|long|FILE|size_t|void)\b"
    r".*;"
)

def subst_globals_in_body(text: str) -> str:
    placeholders = {}
    ctr = [0]
    def protect(m):
        key = f"\x00{ctr[0]:07d}\x00"
        ctr[0] += 1
        placeholders[key] = m.group(0)
        return key + '\x01' * (len(m.group(0)) - len(key))
    protected = PROTECT_RE.sub(protect, text)
    # protect local variable declaration LHS (before '=') from substitution
    # by temporarily replacing them with placeholders
    def protect_decl_lhs(m):
        key = f"\x00D{ctr[0]:06d}\x00"
        ctr[0] += 1
        placeholders[key] = m.group(0)
        return key.ljust(len(m.group(0)), '\x01')

    FUNC_SIG_RE = re.compile(
        r"^(?!#)([A-Za-z_][A-Za-z0-9_ *]*\s+\*?[A-Za-z_][A-Za-z0-9_]*\s*\([^)]*\)\s*)$",
        re.MULTILINE
    )
    LOCAL_DECL_FULL = re.compile(
        r"^([ \t]*(?:register\s+|volatile\s+|const\s+)*"
        r"(?:unsigned\s+|signed\s+|short\s+|long\s+)*"
        r"(?:int|double|float|char|short|long|FILE|size_t|void)\b"
        r"[^=;\n]*?)(?==)",  # LHS up to but not including '='
        re.MULTILINE
    )
    protected = LOCAL_DECL_FULL.sub(protect_decl_lhs, protected)
    # protect function signature lines entirely
    def protect_sig(m):
        key = f"\x00S{ctr[0]:06d}\x00"
        ctr[0] += 1
        placeholders[key] = m.group(0)
        return key.ljust(len(m.group(0)), '\x01')
    protected = FUNC_SIG_RE.sub(protect_sig, protected)

    # substitute all globals in one pass
    result = GLOBAL_RE.sub(lambda m: 'ctx->' + m.group(1), protected)
    for key, val in placeholders.items():
        padded = key + '\x01' * (len(val) - len(key))
        result = result.replace(padded, val)
    return result

# ── Pass 2: rewrite .c files ─────────────────────────────────────────────────
DEF_RE = re.compile(
    r"^((?:unsigned\s+|signed\s+)?[A-Za-z_][A-Za-z0-9_ *]*?)\s+(\*?)([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)",
    re.MULTILINE
)

def rewrite_c(path: Path):
    backup(path)
    raw_lines = path.read_text(errors="replace").splitlines(keepends=True)

    # Step A: add include
    last_inc = max((i for i, l in enumerate(raw_lines) if l.startswith("#include")), default=-1)
    if last_inc >= 0:
        raw_lines.insert(last_inc + 1, '#include "tda_context.h"\n')

    # Step B: remove file-scope global definitions
    # these are at col 0, not inside a function, not a typedef/struct
    cleaned = []
    depth = 0
    removing_block = False
    for l in raw_lines:
        s = l.strip()
        delta = s.count('{') - s.count('}')
        # removal check BEFORE updating depth, using current depth
        if depth == 0:
            if (not s.startswith('#') and not s.startswith('/') and
                    not s.startswith('*') and not s.startswith('typedef') and
                    not s.startswith('struct') and l[0:1] not in (' ', '\t', '\n')):
                m = GDEF_RE.match(s)
                if m and m.group(1) in global_names:
                    # if line opens a block (ends with {), skip until matching }
                    if s.rstrip().endswith('{') or (delta > 0):
                        depth = delta  # track depth to find closing }
                        if depth > 0:
                            removing_block = True
                    continue
        if removing_block:
            depth = max(0, depth + delta)
            if depth <= 0:
                removing_block = False
                depth = 0
            continue
        depth = max(0, depth + delta)
        cleaned.append(l)
    # Post-B: remove orphaned bare type-qualifier lines whose following definition was removed
    cleaned2 = []
    i2 = 0
    while i2 < len(cleaned):
        l2 = cleaned[i2]
        s2 = l2.strip()
        if (s2 in ("unsigned","signed","long","short") and
                (i2 + 1 >= len(cleaned) or not cleaned[i2+1].strip() or
                 cleaned[i2+1].strip().startswith('/*') or
                 not re.match(r'^[A-Za-z_*]', cleaned[i2+1].strip()))):
            i2 += 1
            continue  # drop orphaned qualifier
        cleaned2.append(l2)
        i2 += 1
    raw_lines = cleaned2

    text = "".join(raw_lines)

    # Pre-step: join multi-line function signatures so DEF_RE can match them
    # Pattern: line ending without ; or { that has unbalanced parens
    joined_lines = []
    i2 = 0
    text_lines = text.splitlines(keepends=True)
    while i2 < len(text_lines):
        line = text_lines[i2]
        raw = line.rstrip('\r\n')
        # if line starts at col 0, has '(' but no ')', not a preprocessor/comment/decl
        if (raw and raw[0] not in (' ','\t','#','/','*') and
                '(' in raw and ')' not in raw and ';' not in raw):
            # join next lines until ) found
            combined = raw
            i2 += 1
            while i2 < len(text_lines) and ')' not in combined:
                combined += ' ' + text_lines[i2].strip()
                i2 += 1
            joined_lines.append(combined + '\n')
        else:
            joined_lines.append(line)
            i2 += 1
    text = "".join(joined_lines)

    # Step C: add ctx to function definitions
    def replace_def(m):
        name = m.group(3)
        if name not in func_names or "TDAContext" in m.group(4):
            return m.group(0)
        after = text[m.end():m.end()+200].lstrip('\r\n \t')
        after = re.sub(r'\s*/\*[^*]*\*/\s*', ' ', after).lstrip('\r\n \t')
        if not after.startswith('{'):
            return m.group(0)
        start = m.start()
        if start > 0 and text[start-1] not in ('\n', '\r'):
            return m.group(0)
        return f"{m.group(1)} {m.group(2)}{name}({add_ctx(m.group(4))})"
    text = DEF_RE.sub(replace_def, text)

    # Step D: add ctx to call sites (skip declaration lines ending with ;)


    out_lines = []
    for line in text.splitlines(keepends=True):
        s = line.rstrip('\r\n')
        stripped = s.strip()
        # skip: preprocessor, comments, and declaration lines
        # do NOT skip call-statement lines even though they end with ;
        is_decl = bool(re.match(
            r"^(?:(?:extern|static|register|volatile|const|unsigned|signed|short|long)\s+)*"
            r"(?:int|double|float|char|short|long|void|FILE|size_t|struct|union|enum)\b.*\(",
            stripped
        ))
        if (stripped.startswith('#') or
                stripped.startswith('//') or stripped.startswith('/*') or
                (stripped.startswith('*') and stripped[1:2] in ('/', ' ', '\t', '\n', '')) or
                re.search(r'TDAContext\s*\*\s*ctx', s) or
                'TDAContext' in s or
                is_decl):
            out_lines.append(line)
            continue
        def _call_sub(m):
            # skip matches inside string literals
            before = line[:m.start()]
            in_str = False
            i2 = 0
            while i2 < len(before):
                if before[i2] == '\\': i2 += 2; continue
                if before[i2] == '"': in_str = not in_str
                i2 += 1
            return m.group(1) + '(ctx, ' if not in_str else m.group(0)
        line = CALL_RE.sub(_call_sub, line)
        out_lines.append(line)
    text = "".join(out_lines)
    text = re.sub(r'\(ctx,\s*\)', '(ctx)', text)

    # Step C2: update in-file forward declarations (end with ;)
    DECL_RE2 = re.compile(
        r"^([A-Za-z_][A-Za-z0-9_ *]*?\s+\*?)([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)\s*;",
        re.MULTILINE
    )
    def replace_indecl(m):
        name = m.group(2)
        if name not in func_names or "TDAContext" in m.group(3):
            return m.group(0)
        return m.group(1) + name + "(" + add_ctx(m.group(3)) + ");"
    text = DECL_RE2.sub(replace_indecl, text)

    # Step E: substitute bare globals with ctx-> inside function bodies
    # split into: file-scope segments vs function-body segments
    segments = []
    depth = 0
    in_body = False
    buf = []
    prev_sig = False

    for line in text.splitlines(keepends=True):
        s = line.rstrip('\r\n')
        opens  = s.count('{')
        closes = s.count('}')

        if not in_body and opens > closes:
            on_same_line = '(' in s and '{' in s and s.index('(') < s.index('{')
            if not on_same_line and not prev_sig:
                buf.append(line)
                depth = max(0, depth + opens - closes)
                prev_sig = False
                continue
            segments.append(('file', ''.join(buf)))
            buf = [line]
            in_body = True
            depth = opens - closes
            prev_sig = False
            continue
        elif in_body:
            buf.append(line)
            depth += opens - closes
            if depth <= 0:
                segments.append(('body', ''.join(buf)))
                buf = []
                in_body = False
                depth = 0
        else:
            buf.append(line)

        if s and s[0] not in (' ', '\t') and '(' in s and ')' in s:
            if not s.rstrip().endswith(';') and '{' not in s and '}' not in s:
                prev_sig = True
            else:
                prev_sig = False
        elif s.strip():
            prev_sig = False
    if buf:
        segments.append(('file' if not in_body else 'body', ''.join(buf)))

    result_parts = []
    for kind, chunk in segments:
        if kind == 'body':
            result_parts.append(subst_globals_in_body(chunk))
        else:
            result_parts.append(chunk)

    path.write_text(''.join(result_parts))

import time as _time
print("starting C rewrite loop", file=sys.stderr)
for c in sorted(src.glob("*.c")):
    if c.name in SKIP_C:
        continue
    _t0=_time.time(); rewrite_c(c); _e=_time.time()-_t0
    sys.stderr.write(f"  {c.name} {_e:.1f}s\n")
print("C done", file=sys.stderr)

# ── Pass 3: rewrite .h files ─────────────────────────────────────────────────
DECL_RE = re.compile(
    r"^([A-Za-z_][A-Za-z0-9_ *]*?\s+\*?)([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)\s*;",
    re.MULTILINE
)
EXTERN_RE = re.compile(
    r"^\s*extern\s+\S[^(]*?(\b[A-Za-z_][A-Za-z0-9_]*)\s*(?:\[[^\]]*\])*\s*;"
)

def rewrite_h(path: Path):
    backup(path)
    text = path.read_text(errors="replace")

    # strip extern declarations for globals now in struct
    lines = text.splitlines(keepends=True)
    new_lines = []
    for l in lines:
        m = EXTERN_RE.match(l)
        if m:
            # check all identifiers in the extern line
            all_names = re.findall(r"\b([A-Za-z_][A-Za-z0-9_]+)\b", m.group(0))
            _TW2 = {"extern","int","double","float","char","short","long","void",
                    "FILE","size_t","unsigned","signed","const","struct"}
            var_names = [n for n in all_names if n not in _TW2]
            if var_names and all(n in global_names for n in var_names):
                continue  # skip entire extern line
        new_lines.append(l)
    text = "".join(new_lines)

    # update function declarations
    def replace_decl(m):
        name = m.group(2)
        if name not in func_names or "TDAContext" in m.group(3):
            return m.group(0)
        return m.group(1) + name + "(" + add_ctx(m.group(3)) + ");"
    text = DECL_RE.sub(replace_decl, text)

    # add forward decl if needed
    if "TDAContext" in text and 'tda_ctx_fwd.h' not in text:
        hlines = text.splitlines(keepends=True)
        insert_at = 0
        for i, l in enumerate(hlines):
            s = l.strip()
            if s.startswith("#define _") or (s.startswith("#define ") and len(s.split()) == 2):
                insert_at = i + 1
                break
        hlines.insert(insert_at, '#include "tda_ctx_fwd.h"\n')
        text = "".join(hlines)

    path.write_text(text)

for h in sorted(src.glob("*.h")):
    if h.name in SKIP_H:
        continue
    rewrite_h(h)
print("H done", file=sys.stderr)
