#!/usr/bin/env python3
"""TDA regression test runner.

Usage: python3 tests/check.py src/tda SUITE_DIR [SUITE_DIR ...] [--regen]

Why this is not a plain text diff
---------------------------------
The .ref files were captured with stdout and stderr merged into one pipe.
stdout was fully buffered and stderr was not, so a reference encodes the flush
boundaries of the machine that produced it: change the stdout stream by one
character and every later boundary shifts, and a stderr write can end up
spliced into the middle of a stdout line.  Comparing the merged streams
therefore reports failures that have nothing to do with the program.

So this runner captures the two streams separately, reconstructs the reference
stdout stream by deleting the stderr writes from it, and compares that against
stdout.  Numbers are compared with a relative tolerance rather than as text,
because the references are one machine's floating point rather than analytic
truth, and a few tests are ill-conditioned by construction (see the tables
below).
"""
import subprocess, sys, re
from pathlib import Path

binary = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./tda")
if not binary.exists():
    sys.exit("ERROR: binary %s not found" % binary)

REGEN = "--regen" in sys.argv
argv = [a for a in sys.argv if a != "--regen"]

SKIP_PREFIXES = (
    "TDA. Analysis",
    "Current time:",
    "Current memory",
    "Max memory",
    "End of program.",
)
SKIP_CONTAINS = (
    "Current memory:",
    "Function Value     Norm of Gradient",
    "Norm of final gradient vector:",
    "Last absolute change of function value:",
    "Last relative change in parameters:",
    "Number of function evaluations:",
    "Convergence reached in",
    "Final scaled gradient:",
    "Final scaled parameter change:",
)

_ITER_LINE = re.compile(r"^\s+\d+\s+[\d.]+e[+\-]\d{2}\s+[\d.]+e[+\-]\d{2}")
_NUM = re.compile(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?")
_NEGZERO = re.compile(r"-(0\.0+)(?![1-9])")
_CREATION = re.compile(rb"%%CreationDate:[^\n]*")

# Deleting stderr writes one by one keeps them in stream order, which drifts if
# the reference run made one more or one fewer write than this one.  These
# few printfe() formats are distinctive enough to delete on sight afterwards,
# which mops up whatever the ordered pass missed.
STDERR_RESIDUE = [re.compile(p) for p in (
    # printfe() writes these as "%7d     %c": a wide right-aligned field and
    # five trailing blanks before the carriage return.  printf1() also prints a
    # "Read records: N" line to stdout, but narrow and without the padding, so
    # the trailing blanks are what tells the two apart.
    r"(?:Read records|Records written|Component|Done|Node|Tree):\s+\d+[ \t]{4,}",
    r"Sorting keys \.\.\.",
    r"Sorting: \S+",
    r"Iter\s*\d+\s+Crit\s+[-\d.]+e[-+]\d+",
    r"Iteration\s*\d+\s+Crit:\s+[-\d.]+e[-+]\d+",
    # The ML iteration table is written to stderr too.  The row pattern keys on
    # the exact field widths of its "%3d  %20.13e %17.10e" format so it cannot
    # match a result table on stdout.
    r"\s*Iter\s+Function Value\s+Norm of Gradient\s+Par Change\s+FCall",
    r"\s*\d+\s+[-+]?\d\.\d{13}e[-+]\d{2}"
    r"(?:\s+[-+]?\d\.\d{10}e[-+]\d{2})*"
    r"(?:\s+(?:--|[-+]?\d\.\d{4}e[-+]\d{2}))?"
    r"(?:\s+\d+(?:\s*\(\d+,\d+\))?)?[ \t]*",
)]


def strip_residue(text):
    """Delete leftover stderr writes, including ones embedded in a stdout line.

    The write is removed together with the break that ended it, so a stdout
    line it was spliced into is put back together rather than left in two
    pieces.  Only a match that runs to the end of a line is removed, so a
    genuine stdout line is never cut in half.
    """
    text = text.replace("\r\n", "\n")
    for pat in STDERR_RESIDUE:
        out, pos = [], 0
        for m in pat.finditer(text):
            if m.start() < pos:
                continue
            j = m.end()
            if j < len(text) and text[j] not in "\r\n":
                continue
            out.append(text[pos:m.start()])
            pos = j + 1 if j < len(text) else j
        out.append(text[pos:])
        text = "".join(out)
    return text

# Relative tolerance for numeric fields.  Text outside numbers must match
# exactly; digits may drift because the references are snapshots.
RTOL = 1e-6

# Relative tolerance cannot compare a value against zero, so a result that
# should be exactly zero and comes out at 1e-8 needs an absolute floor as
# well.  Output carries 10 decimals, so this sits an order of magnitude above
# the print resolution and well below anything meaningful.
ATOL = 1e-9

# Two numbers one unit apart in their last printed digit count as equal, up to
# this relative bound.
#
# This is not about IEEE rounding.  A derived column can be printed to more
# significant digits than the quantities it derives from, so it is the first
# place a small genuine disagreement becomes visible.  qr7 is the case: it
# maximises a likelihood, stops on the scaled-gradient test (TOLSG 1e-5,
# reached at ~5e-8), and reports a final scaled parameter change of 6.3e-4 --
# the parameter test (TOLSP 1e-8) is never met, so the iterate is only located
# to ~1e-4 in scaled parameter space.  With a numerical gradient and a flat
# direction in the model (Sigma 2,1 is 0.4134 with a standard error of 0.6320)
# the estimates land 1e-5 apart on different libm implementations.  Coeff and
# Error still print identically at 4 decimals; C/Error, printed to 5
# significant digits, moves by one in its last place.  Moving it that far takes
# a coefficient change of 3.5e-6 or a standard error change of 6.8e-7, both
# ~100x finer than one unit in their own last printed place.  The log
# likelihood agrees to all 6 printed digits.
ULP_RTOL = 1e-3

# Tests that are ill-conditioned by construction, with the tolerance each needs
# and the reason it needs it.
RTOL_OVERRIDE = {
    # Binary model fitted by ML.  Iteration 0 agrees to 11 digits, but the line
    # search amplifies the last bit, so standard errors drift in the 4th digit.
    "qr6.cf": 5e-3,
}

# An estimate row from an ML fit: index, name, coefficient, standard error,
# coefficient/error, significance.
_EST_ROW = re.compile(
    r"^\s*\d+\s+(\S+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")

# A coefficient smaller than a tenth of its own standard error sits in a flat
# direction of the likelihood: the optimiser stops wherever its convergence
# test happens to fire, so the value carries no information and cannot be
# compared across builds.  cd2 (and test, which repeats it) fit a Negbin II
# whose overdispersion parameter is unidentified this way -- the maximised log
# likelihood and every other coefficient still have to match exactly.
FLAT_RATIO = 0.1

# Eigenvectors are defined only up to sign, and for the repeated eigenvalues in
# these matrices only up to a rotation of the eigenspace, so which sign
# convention a solver lands on is arbitrary.  For these a line also counts as
# matching when every number in it is negated.
SIGN_FREE = {"mat5.cf", "test.cf"}


def lines_of(text):
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    out = []
    for line in text.split("\n"):
        line = _NEGZERO.sub(r" \1", line.rstrip())
        if not line:
            continue
        if any(line.startswith(s) for s in SKIP_PREFIXES):
            continue
        if any(s in line for s in SKIP_CONTAINS):
            continue
        if _ITER_LINE.match(line):
            continue
        out.append(line)
    return out


def shape_re(line):
    """Regex matching `line` with every numeric run wildcarded."""
    parts, last = [], 0
    for m in _NUM.finditer(line):
        parts.append(re.escape(line[last:m.start()]))
        parts.append(r"[-+]?[\d.]+(?:[eE][-+]?\d+)?")
        last = m.end()
    parts.append(re.escape(line[last:]))
    return re.compile("".join(parts))


_MERGED_MARK = re.compile(
    r"(?:Read records|Records written):\s+\d+[ \t]{4,}|^Sorting keys"
    r"|^\s*Iter\s+Function Value\s*$", re.M)


def is_merged(ref_text):
    """True if the reference still has stderr writes mixed into it."""
    return bool(_MERGED_MARK.search(ref_text))


def recover_stdout(ref_text, err_text):
    """Delete the stderr writes from a merged reference, in stream order.

    Each write is matched at the leftmost position it can occupy, preferring a
    literal match over a shape match (numeric runs wildcarded, because
    iteration traces carry platform-dependent digits) and a line start over a
    mid-line position when both start in the same place.  The terminating line
    break goes with the write: it belongs to stderr, and leaving it behind
    splits whichever stdout line the write landed inside.
    """
    ref = ref_text.replace("\r\n", "\n")
    pos = 0
    for raw in re.split(r"[\r\n]", err_text.replace("\r\n", "\n")):
        if not raw.strip():
            if pos < len(ref) and ref[pos] in "\r\n":
                ref = ref[:pos] + ref[pos + 1:]
            continue
        body = shape_re(raw).pattern + r"(?=[\r\n])"
        cands = []
        for rank, pat in enumerate((
                r"(?:^|(?<=\n))" + re.escape(raw) + r"\s*(?=[\r\n])",
                r"(?:^|(?<=\n))" + body,
                re.escape(raw) + r"\s*(?=[\r\n])",
                body)):
            m = re.compile(pat).search(ref, pos)
            if m:
                cands.append((m.start(), rank, m.span()))
        if not cands:
            continue
        i, j = min(cands)[2]
        if j < len(ref) and ref[j] in "\r\n":
            j += 1
        ref = ref[:i] + ref[j:]
        pos = i
        # A stderr write can land in the middle of a stdout line that has
        # not reached its own newline yet -- stdio buffers can flush
        # mid-write, so the boundary the merged capture shows there is not
        # a real stdout newline at all, just where the buffer happened to
        # end. Deleting the write (just above) correctly removes it, but
        # leaves that not-really-there newline behind, splitting one
        # stdout line (a boilerplate dash separator, in every case seen)
        # into two consecutive dash-only lines. Confirmed directly: TDA's
        # own actual stdout for the same run has this as a single
        # unbroken line. Checked only right at the position the deletion
        # just happened, not swept across the whole file -- an earlier,
        # unanchored version of this fix merged unrelated dash lines
        # elsewhere in the file that happened to end up adjacent for
        # other reasons, breaking two different tests it had no business
        # touching.
        before_m = re.search(r"(-+)\n\Z", ref[:i])
        after_m = re.match(r"(-+)(?=[\r\n])", ref[i:])
        if before_m and after_m:
            start = i - len(before_m.group(0))
            end = i + len(after_m.group(1))
            ref = ref[:start] + before_m.group(1) + after_m.group(1) + \
                ref[end:]
            pos = start
    return ref


def negate(line):
    return _NUM.sub(lambda m: m.group()[1:] if m.group().startswith("-")
                    else "-" + m.group(), line)


def print_ulp(token):
    """One unit in the last printed digit of `token`, or None if it has none."""
    m = re.match(r"[-+]?\d*\.(\d+)$", token)
    return 10.0 ** -len(m.group(1)) if m else None


def _cmp(a, b, rtol):
    if a == b:
        return True
    na, nb = list(_NUM.finditer(a)), list(_NUM.finditer(b))
    if len(na) != len(nb) or _NUM.sub("#", a) != _NUM.sub("#", b):
        return False
    for ma, mb in zip(na, nb):
        try:
            x, y = float(ma.group()), float(mb.group())
        except ValueError:
            return False
        if x == y:
            continue
        if abs(x - y) <= max(rtol * max(abs(x), abs(y)), ATOL):
            continue
        # One unit in the last printed digit, bounded by ULP_RTOL so it
        # cannot excuse a value near zero moving from 0.0001 to 0.0002 --
        # also one printed unit, but a 67% change.  See ULP_RTOL above for
        # why a derived column moves when its inputs do not.
        u = print_ulp(ma.group())
        if (u is not None and abs(x - y) <= u * 1.000001 and
                abs(x - y) <= ULP_RTOL * max(abs(x), abs(y))):
            continue
        return False
    return True


def unidentified(a, b):
    ma, mb = _EST_ROW.match(a), _EST_ROW.match(b)
    if not (ma and mb) or ma.group(1) != mb.group(1):
        return False
    try:
        return (abs(float(ma.group(4))) < FLAT_RATIO and
                abs(float(mb.group(4))) < FLAT_RATIO)
    except ValueError:
        return False


def line_match(a, b, rtol, sign_free=False):
    if _cmp(a, b, rtol):
        return True
    if unidentified(a, b):
        return True
    if sign_free:
        # Flipping a sign also shifts the column padding.
        ws = lambda t: re.sub(r"\s+", " ", t).strip()
        return _cmp(ws(a), ws(negate(b)), rtol)
    return False


def lines_match(want, got, rtol, sign_free=False):
    """Walk both lists, allowing consecutive `want` lines to be joined.

    Reconstruction can still leave a stdout line split where a burst of stderr
    writes landed inside it, so a `got` line is also accepted when it equals
    the concatenation of the next few `want` lines.
    """
    i = j = 0
    while i < len(want) and j < len(got):
        if line_match(want[i], got[j], rtol, sign_free):
            i += 1; j += 1
            continue
        # Trailing blanks were stripped from each fragment, so the rejoined
        # line is compared with whitespace collapsed.
        ws = lambda t: re.sub(r"\s+", " ", t).strip()
        frags, k, hit = [want[i]], i + 1, False
        while k < len(want) and len("".join(frags)) < len(got[j]) + 2:
            frags.append(want[k])
            k += 1
            for cand in ("".join(frags), " ".join(frags)):
                if (line_match(cand, got[j], rtol, sign_free) or
                        line_match(ws(cand), ws(got[j]), rtol, sign_free)):
                    hit = True
                    break
            if hit:
                break
        if not hit:
            return False
        i, j = k, j + 1
    return i == len(want) and j == len(got)


def first_diff(want, got, rtol, sign_free):
    d = next(((n, a, b) for n, (a, b) in enumerate(zip(want, got))
              if not line_match(a, b, rtol, sign_free)), None)
    if d is None:
        d = (min(len(want), len(got)),
             "<%d lines>" % len(want), "<%d lines>" % len(got))
    return d


def write_diff(binary, name, want, got, rtol, sign_free):
    """Write the first differing lines to tda_diffs/ for the CI artifact."""
    d = binary.parent / "tda_diffs"
    d.mkdir(exist_ok=True)
    rows = ["--- want (reference, stderr removed)",
            "+++ got  (this build, stdout only)",
            "reference: %d lines, output: %d lines" % (len(want), len(got)), ""]
    shown = 0
    for n, (a, b) in enumerate(zip(want, got)):
        if line_match(a, b, rtol, sign_free):
            continue
        rows += ["line %d" % n, "  want %r" % a, "  got  %r" % b]
        shown += 1
        if shown == 20:
            rows.append("... further differences suppressed")
            break
    if not shown:
        rows.append("lines all match pairwise; the two differ in length")
    (d / (name + ".diff")).write_text("\n".join(rows) + "\n")


def run_suite(name, suite_dir):
    suite = Path(suite_dir)
    if not suite.exists():
        print("  (skipping %s: %s not found)" % (name, suite_dir))
        return 0, []
    # Every case runs with cwd set to the suite directory, so a `df=`
    # output lands in the SOURCE TREE.  Some of that is deliberate --
    # cl1.df, cl3.df and ehi12.bl are written by one case and read by a
    # later one -- so those files cannot be deleted as they appear.
    #
    # Leaving them behind, though, makes the suite history-dependent: a
    # residue file from an earlier run (or an older binary) can be read
    # by a case before its producer runs.  That is the likeliest cause
    # of the one-off 490/491 with "Error in checking transitions" seen
    # once and never reproduced.
    #
    # So snapshot what is here first and delete anything created during
    # the run, once the suite finishes.  Dependencies still work while
    # it runs; the tree is left as it was found.
    before = {f.name: f.read_bytes() for f in suite.iterdir() if f.is_file()}

    ok, fails, skipped, regenerated = 0, [], 0, []
    for ref in sorted(suite.glob("*.ref")):
        cf = ref.with_suffix(".cf")
        if not cf.exists():
            continue

        # Skip a case whose input is missing.  The third-party SPSS, Excel,
        # shapefile and ArcInfo samples are not distributed with the sources,
        # so the cases that read them only run when those files are present.
        missing = False
        for tok in re.findall(r"[=\s]([\w.]+\.(?:sav|por|xls|dta|dbf|shp|e00|zoo|zad))",
                              cf.read_text(errors="replace")):
            if not (suite / tok).exists():
                missing = True
        for tok in re.findall(r"sdshp\([^)]*\)\s*=\s*(\w+)",
                              cf.read_text(errors="replace")):
            if not (suite / (tok + ".shp")).exists():
                missing = True
        if missing:
            skipped += 1
            continue
        try:
            r = subprocess.run([str(binary.resolve()), "cf=" + cf.name],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               cwd=str(suite.resolve()), timeout=120)
        except subprocess.TimeoutExpired:
            fails.append((cf.name, (0, "timeout", "timeout")))
            continue

        got = lines_of(strip_residue(r.stdout.decode(errors="replace")))
        ref_text = ref.read_text(errors="replace")
        # References captured with the two streams merged need the stderr
        # writes taken back out; ones captured from stdout alone must be left
        # exactly as they are, or the same deletion would eat real output.
        if is_merged(ref_text):
            ref_text = recover_stdout(ref_text, r.stderr.decode(errors="replace"))
        want = lines_of(strip_residue(ref_text))
        rtol = RTOL_OVERRIDE.get(cf.name, RTOL)
        sf = cf.name in SIGN_FREE

        if REGEN:
            # Rewrite a reference only when its compared content moved.
            # Every run stamps the date and the memory figures into the
            # header and footer, which the comparison ignores; writing
            # those back unconditionally touched all 650 files per regen.
            if not lines_match(want, got, rtol, sf):
                ref.write_bytes(r.stdout)
                regenerated.append(ref.name)
            ok += 1
            continue

        if lines_match(want, got, rtol, sf):
            ok += 1
        else:
            fails.append((cf.name, first_diff(want, got, rtol, sf)))
            write_diff(binary, cf.stem, want, got, rtol, sf)
    # Restore the tree: delete anything this run created, and put back
    # what it overwrote.  See the note at the top of run_suite -- the
    # dependencies between cases work while it runs, and nothing is left
    # behind to influence the next.  Sixty of the suite's PostScript
    # files are committed outputs that every run rewrites, and TDA stamps
    # %%CreationDate into each; a checking run restores every file it
    # touched, a --regen run keeps a file whose drawing changed and
    # restores one that differs by the date line alone.
    for f in sorted(suite.iterdir()):
        if not f.is_file():
            continue
        if f.name not in before:
            try:
                f.unlink()
            except OSError:
                pass
            continue
        old = before[f.name]
        new = f.read_bytes()
        if new == old:
            continue
        if REGEN and _CREATION.sub(b"", new) != _CREATION.sub(b"", old):
            if f.name not in regenerated:
                regenerated.append(f.name)
            continue
        f.write_bytes(old)

    print("%s: %d/%d%s" % (name, ok, ok + len(fails),
                           "  (%d skipped, input not present)" % skipped if skipped else ""))
    if REGEN:
        print("  regenerated %d reference(s)%s"
              % (len(regenerated), ": " + " ".join(regenerated) if regenerated else ""))
    for nm, d in fails:
        print("  FAIL %-12s line %d\n    want %r\n    got  %r"
              % (nm, d[0], d[1], d[2]))
    return ok, fails


# Any number of suite directories, named after the directory rather than
# assumed to be exam and ehhnew.
suites = argv[2:] if len(argv) > 2 else ["exam", "ehhnew"]

t_ok = t_fail = 0
for d in suites:
    o, f = run_suite(Path(d).name, d)
    t_ok += o
    t_fail += len(f)
print("Total: %d/%d" % (t_ok, t_ok + t_fail))
sys.exit(1 if t_fail else 0)
