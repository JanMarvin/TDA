#!/usr/bin/env python3
"""Option sweep: run every documented enumerated option value of every
covered command, and report anomalies.

For each command that already has a working .cf case, the sweep takes
that case, and for each integer option whose syntax box enumerates its
values (lines like "3 = ..." under "opt=...," / "alg=...,"), re-runs the
case with the option forced to each documented value in turn.  One
variant per run, because TDA halts a command file at the first error.

The output is a report, not reference files: variants that abort, print
new errors, or collapse to far less output than the base run are the
ones worth a human look.  This is the mechanical version of the sweep
that found the projection-20 defect.

    python3 tools/option_sweep.py [command ...]     sweep some or all
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TDA = ROOT / "src" / "tda"
SUITES = [ROOT / "examples" / s for s in ("tests", "coverage", "exam", "ehhnew")]

def dispatched():
    src = (ROOT / "src" / "t_cmd.c").read_text(errors="replace")
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.DOTALL)
    return [m.group(1) for m in
            re.finditer(r'strncmp\(p,"([a-z0-9_#]+)",\d+\)', src)]

def syntax_box(name):
    for f in sorted((ROOT / "src").glob("t_*.c")):
        s = f.read_text(errors="replace")
        m = re.search(r"int %s\(TDAContext \*ctx\)\s*\n\{" % re.escape(name), s)
        if not m:
            continue
        lines, box = s[:m.start()].splitlines(), []
        for line in reversed(lines):
            t = line.strip()
            if not t or t.startswith("/*") or t.endswith("*/"):
                box.append(t.strip("/*").rstrip("*/"))
                if re.match(r"-{10,}", t.strip("/*# ")):
                    break
            else:
                break
        return list(reversed(box))
    return []

def enum_options(box):
    """option name -> sorted list of documented integer values"""
    opts, cur = {}, None
    for line in box:
        m = re.match(r"\s*([a-z][a-z0-9]*)\s*=\s*\.\.\.", line)
        if m:
            cur = m.group(1)
            opts.setdefault(cur, set())
            d = re.search(r"def\.?\s*=?\s*(\d+)", line)
            if d:
                opts[cur].add(int(d.group(1)))
            continue
        if cur:
            # value lines: a small integer, two or more spaces, text
            v = re.match(r"\s*(\d+)\s*(?:[=:]\s*\S|\s{2,}\S)", line)
            if v:
                opts[cur].add(int(v.group(1)))
            elif not line.strip():
                cur = None
    return {k: sorted(v) for k, v in opts.items() if len(v) >= 2}

def find_case(name):
    """first .cf in the corpus whose text invokes the command"""
    pat = re.compile(r"(?m)^\s*%s\s*[(;=]" % re.escape(name))
    for suite in SUITES:
        for cf in sorted(suite.glob("*.cf")):
            if cf.name.startswith("err"):
                continue        # refusal cases are no base for a sweep
            t = re.sub(r"#[^\n]*", "", cf.read_text(errors="replace"))
            if pat.search(t):
                return cf
    return None

def force_option(text, name, opt, val):
    """rewrite the FIRST invocation of the command to carry opt=val"""
    m = re.search(r"(?m)^(\s*)(%s)\s*(\(|;|=)" % re.escape(name), text)
    if not m:
        return None
    if m.group(3) == "(":
        i = text.index("(", m.start())
        depth, j = 0, i
        while j < len(text):
            if text[j] == "(":
                depth += 1
            elif text[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        inner = text[i + 1:j]
        # replace existing option or append
        new_inner, n = re.subn(r"(?<![a-z0-9])%s\s*=\s*[^,()]+" % re.escape(opt),
                               "%s=%d" % (opt, val), inner, count=1)
        if n == 0:
            new_inner = inner.rstrip()
            if new_inner and not new_inner.endswith(","):
                new_inner += ", "
            new_inner += "%s=%d" % (opt, val)
        return text[:i + 1] + new_inner + text[j:]
    # no parentheses: add them
    return (text[:m.end(2)] + "(%s=%d)" % (opt, val) + text[m.end(2):])

def run_cf(workdir, cf_name):
    r = subprocess.run([str(TDA), "cf=" + cf_name], cwd=workdir,
                       capture_output=True, text=True, timeout=120)
    out = r.stdout + r.stderr
    return r.returncode, out

def main():
    only = set(sys.argv[1:])
    report = []
    for name in dispatched():
        if only and name not in only:
            continue
        box = syntax_box(name)
        opts = enum_options(box)
        if not opts:
            continue
        cf = find_case(name)
        if cf is None:
            continue
        with tempfile.TemporaryDirectory() as td:
            for f in cf.parent.iterdir():
                if f.is_file() and f.suffix != ".ref":
                    shutil.copy(f, td)
            base_rc, base_out = run_cf(td, cf.name)
            base_err = [l for l in base_out.splitlines() if "Error" in l]
            for opt, vals in sorted(opts.items()):
                for val in vals:
                    text = cf.read_text(errors="replace")
                    variant = force_option(text, name, opt, val)
                    if variant is None:
                        continue
                    vname = "sweep_%s_%s_%d.cf" % (name, opt, val)
                    (Path(td) / vname).write_text(variant)
                    try:
                        rc, out = run_cf(td, vname)
                    except subprocess.TimeoutExpired:
                        report.append((name, opt, val, "TIMEOUT", ""))
                        continue
                    errs = [l for l in out.splitlines()
                            if "Error" in l and l not in base_err]
                    tag = None
                    if rc not in (0, base_rc):
                        tag = "rc=%d" % rc
                    elif errs:
                        tag = "new-error"
                    elif len(out.splitlines()) < len(base_out.splitlines()) // 3:
                        tag = "short-output"
                    if tag:
                        report.append((name, opt, val, tag,
                                       errs[0].strip() if errs else ""))
    for name, opt, val, tag, detail in report:
        print("%-10s %-6s=%-3d  %-12s %s" % (name, opt, val, tag, detail))
    print("%d anomalies" % len(report))

if __name__ == "__main__":
    main()
