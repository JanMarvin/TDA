#!/usr/bin/env python3
"""Run every TDA command against real input and report crashes and hangs.

smoke.py calls each command with no arguments, so most of them stop at the
first check ("need right-hand side", "no data matrix") and almost no code runs.
This gives each command something to work on -- a data matrix, episode data, a
graph, a PostScript page -- and escalates through argument forms until one gets
past the input checks.  It still does not test what anything computes; it tests
that the command survives being handed plausible input.

That is the bar `rfit`/`rfit1` failed, and unlike the 350 worked examples it
covers every command rather than the 93 the textbook happened to use.

Usage: python3 tests/smoke2.py src/tda [-v] [--full] [--only cmd,cmd]
"""
import json, os, re, subprocess, sys, tempfile
from pathlib import Path

# make_tests.py imports this module for its fixtures and preludes, so the
# binary is only required when the sweep actually runs.
binary = Path(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-")
              else "./tda").resolve()
if __name__ == "__main__" and not binary.exists():
    sys.exit("ERROR: binary %s not found" % binary)
VERBOSE = "-v" in sys.argv
FULL = "--full" in sys.argv

# Combinations that worked last time.  Without it every pass rediscovers them
# from scratch, and the sweep costs five minutes instead of one.
CACHE = Path(__file__).with_name("smoke2_cache.json")
try:
    known = json.loads(CACHE.read_text())
except Exception:
    known = {}

# No attempt cap: a working combination can sit anywhere in the list, and
# capping at 250 lost 74 commands whose form happened to be late.  The cache is
# what makes a pass fast; a cold run still costs a few minutes.
only = None
if "--only" in sys.argv:
    only = set(sys.argv[sys.argv.index("--only") + 1].split(","))

src = (Path(__file__).parent.parent / "src" / "t_cmd.c").read_text(errors="replace")
lines, keep, xwin = src.split("\n"), [], 0
for l in lines:
    if re.match(r"\s*#if\s+S_XWIN", l):
        xwin = 1
    elif re.match(r"\s*#endif", l) and xwin:
        xwin = 0
    elif not xwin and not re.match(r"\s*/\*", l):
        keep.append(l)
CMDS = sorted(set(re.findall(r'strn?cmp\(p,"([a-z][a-z0-9_]*)"', "\n".join(keep))))

work = tempfile.mkdtemp(prefix="tdasmoke2")

# The help command looks for tda.hlp in the working directory, and the format
# readers need real files -- there is no way to synthesise an SPSS system file
# from a command line.  fixtures/ carries an .sav and a .por (from the readspss
# project), a small .dta and a hand-built .dbf.
_hlp = Path(__file__).parent.parent / "src" / "tda.hlp"
if _hlp.exists():
    Path(work, "tda.hlp").write_bytes(_hlp.read_bytes())
_fix = Path(__file__).parent / "fixtures"
if _fix.is_dir():
    for f in _fix.iterdir():
        Path(work, f.name).write_bytes(f.read_bytes())

# Third-party input files -- the SPSS, Excel, shapefile and ArcInfo samples --
# are not ours and are not GPL, so they are not distributed with the sources.
# Point TDA_EXT_INPUT at a directory holding them, or drop them in
# ../TDA_ext_input next to the repository, and the commands that read them are
# exercised; otherwise those commands are reported separately rather than
# counted as failures.
EXT = Path(os.environ.get(
    "TDA_EXT_INPUT", Path(__file__).parent.parent.parent / "TDA_ext_input"))
HAVE_EXT = EXT.is_dir() and any(EXT.iterdir())
if HAVE_EXT:
    for f in EXT.iterdir():
        if f.is_file():
            Path(work, f.name).write_bytes(f.read_bytes())

# Commands that cannot run without those files.
EXT_ONLY = {"rspss", "rxls", "sdshp", "sde00", "sdinf", "sdsel", "sdvd",
            "sdclip", "sdencl", "sdppol", "sdlpol", "sdipol", "sdplot",
            "sdpmap", "sdpgeo", "sdnvar", "sdpdata", "sdrel", "sdnl",
            "sddcwp", "sdgshhs", "sdplot31", "sdplot32", "sdplot33"}

# A small data file: an id, two covariates, a group, a duration and a status.
# 20 cases, not 200: ghd1 searches for a maximal hierarchy and goes
# exponential in the number of cases somewhere between 22 and 30, which is
# indistinguishable from a hang.
# Columns 7-10 are two censoring intervals, lower then upper.  The interval
# regressions in t_ireg.c read four variables as two (lower, upper) pairs and
# reject any record where upper < lower, so the fixture has to carry properly
# ordered pairs rather than arbitrary columns.
rows = []
for i in range(1, 21):
    yl, xl = (i * 7) % 23 + 1, (i * 13) % 17 + 1
    # Columns 11-13 hold a distribution: subm compares distributions and
    # rejects a variable that is not one.
    rows.append("%3d %6.2f %6.2f %d %3d %d %6.2f %6.2f %6.2f %6.2f"
                " %6.3f %6.3f %6.3f" %
                (i, yl, xl, i % 3 + 1, (i * 11) % 60 + 5, i % 2,
                 yl, yl + 2.0, xl, xl + 3.0, 0.2, 0.3, 0.5))
Path(work, "d.dat").write_text("\n".join(rows) + "\n")

# Upper triangle of a 6x6 dissimilarity matrix.  Small on purpose: gnst
# enumerates every labelled spanning tree, n^(n-2) by Cayley's formula, so ten
# nodes is 10^8 of them and looks exactly like a hang.
tri = []
for i in range(6):
    for j in range(i + 1, 6):
        tri.append("%d" % (abs(i - j) * 3 + 1))
Path(work, "g.dat").write_text("\n".join(tri) + "\n")

# a small connected digraph with a cycle, on six nodes
edges = [(1, 2, 1), (2, 3, 2), (3, 4, 3), (4, 5, 4), (5, 6, 5), (6, 1, 6),
         (2, 5, 7), (1, 4, 8), (3, 6, 9)]
Path(work, "e.dat").write_text(
    "\n".join("%d %d %d" % e for e in edges) + "\n")

# A six-node tree.  t.dat points away from the root for the undirected uses;
# r.dat points towards it, because hclsp wants a directed tree whose root is
# the single node with outdegree 0.
Path(work, "t.dat").write_text(
    "\n".join("%d %d 1" % e for e in
              [(1, 2), (1, 3), (2, 4), (2, 5), (3, 6)]) + "\n")
Path(work, "r.dat").write_text(
    "\n".join("%d %d 1" % e for e in
              [(2, 1), (3, 1), (4, 2), (5, 2), (6, 3)]) + "\n")

NVAR = """nvar(
    dfile = d.dat,
    ID [3.0] = c1,
    X  [6.2] = c2,
    Y  [6.2] = c3,
    G  [1.0] = c4,
    T  [3.0] = c5,
    D  [1.0] = c6,
    YL [6.2] = c7,
    YU [6.2] = c8,
    XL [6.2] = c9,
    XU [6.2] = c10,
    P1 [6.3] = c11,
    P2 [6.3] = c12,
    P3 [6.3] = c13,
    TS [1.0] = 0,
    ORG[1.0] = 0,
    TL [3.0] = T,
    TU [3.0] = T + 5,
    G1 = G[1],
    G2 = G[2],
);
"""
EDEF = "edef(ts = 0, tf = T, org = 0, des = D);\n"
# psfile= opens the output; psetup() defines the coordinate system and takes
# no right-hand side.
PSETUP = """psfile = p.ps;
psetup(
    pxlen = 80, pylen = 40,
    pxa = 0,100, pya = 0,100,
);
"""
# gdd(opt=4) reads the upper triangle of a dissimilarity matrix from a single
# variable, so the graph fixture needs its own data file of n(n-1)/2 values and
# its own nvar -- a second nvar replaces the first, so it cannot be combined
# with the episode fixture.
# psetup3 takes pxlen only: pylen is commented out in t_psf.c and there is no
# pzlen at all.  The extents come from pxa/pya/pza.
# psetup3 parses its own options inline with sscanf("pxlen=%lg", ...) instead
# of going through parm(), so unlike psetup it does not tolerate spaces around
# the '='.  "pxlen = 80" fails with "Error: in number of coordinates".
PSETUP3 = """psfile = p3.ps;
psetup3(pxlen=80, pxa=0,100, pya=0,100, pza=0,100,);
"""

# An edge list: from-node, to-node, value.  gdd() without gt= builds a graph
# from one, which is what the g* commands mostly want -- the dissimilarity
# form below feeds the clustering side instead.
EDGES = """nvar(
    dfile = e.dat,
    I = c1,
    J = c2,
    V = c3,
);
gdd = I,J,V;
"""

EDGESD = EDGES.replace("gdd = I,J,V;", "gdd(gt=2) = I,J,V;")

# Three matrices for the matrix operators.
MATS = """mdef(A,3,3) = 1,2,3, 4,5,6, 7,8,9;
mdef(B,3,3) = 9,8,7, 6,5,4, 3,2,1;
mdef(C,3,3) = 1,0,0, 0,1,0, 0,0,1;
"""

GRAPH = """nvar(
    dfile = g.dat,
    D = c1,
);
gdd(opt=4) = D;
"""

# Several commands need a graph with particular properties: pltree wants an
# undirected tree, hclsp wants one built with gdd option 1, gcni wants
# undirected.  gt=2 is what makes gdd produce an undirected graph.
TREE = """nvar(
    dfile = t.dat,
    I = c1,
    J = c2,
    V = c3,
);
gdd(opt=1, gt=2) = I,J,V;
"""

UNDIR = EDGES.replace("gdd = I,J,V;", "gdd(opt=1, gt=2) = I,J,V;")

ROOTED = """nvar(
    dfile = r.dat,
    I = c1,
    J = c2,
    V = c3,
);
gdd(opt=1) = I,J,V;
"""

# psetupg establishes a geographical coordinate system, which sdpgrat and the
# other graticule commands draw into.  region= has no default and the command
# fails without it.
GEO = """psfile = g.ps;
psetupg(proj=10, region=30,30, pxlen=120);
"""

# xreg, xconh and the other "x" plot decorations attach themselves to the most
# recent scatterplot, so a plot has to exist before they will do anything.
SCATTER = PSETUP + "plot = X,Y;\n"

# Preludes in increasing order of what they provide.
PRELUDES = [
    ("bare",     ""),
    ("data",     NVAR),
    ("episodes", NVAR + EDEF),
    ("graph",    GRAPH),
    ("edges",    EDGES),
    ("digraph",  EDGESD),
    ("plot",     NVAR + PSETUP),
    ("all",      NVAR + EDEF + PSETUP),
    ("matrix",   NVAR + MATS),
    ("scatter",  NVAR + SCATTER),
    ("tree",     TREE),
    ("tree+ps",  TREE + PSETUP),
    ("undir",    UNDIR),
    ("matrix+ps", NVAR + MATS + PSETUP),
    ("dblock",   NVAR + "dblock(id=ID);\n"),
    # arcc and arcv need an open data archive; fixtures/ carries one built by
    # make_archive.py.
    ("archive",  "arcd = tda.zad;\n"),
    ("coverage", "sde00() = sample.e00;\n"),
    ("geo",      GEO),
    ("rooted",   ROOTED),
    ("plot3",    NVAR + PSETUP3),
    ("graph+ps", GRAPH + PSETUP),
    ("edges+ps", EDGES + PSETUP),
]

# Argument forms to try for the command itself.
# Escalating argument shapes.  Many commands want a variable list either
# inside the parentheses or on the right, and some want a number or an
# expression rather than a file name.
FORMS = [
    "%s();",
    "%s() = out.txt;",
    "%s() = 1;",
    "%s;",
    "%s() = X;",
    "%s() = X,Y;",
    "%s(X) = out.txt;",
    "%s(v=X) = out.txt;",
    "%s(tp = 0 (10) 100) = out.txt;",
    "%s(sc=50) = out.txt;",
    "%s() = D;",
    # A regression-shaped right-hand side: dependent first, then regressors.
    # lsreg is written this way in the examples ("lsreg = LogSurv,Time;").
    "%s = Y,X;",
    "%s() = Y,X;",
    "%s(df=df) = Y;",
    "%s() = I,J,V;",
    "%s(opt=1) = out.txt;",
    "%s(tp = 0 (10) 100) = X;",
    # fml and its relatives take a parameter definition on the right, not a
    # file; glm takes its variables through v=.
    "%s() = f = X * b1 + b2;",
    "%s(v=X,Y) = out.txt;",
    "%s(v=X,Y);",
    "%s() = 2;",
    "%s() = X + Y;",
    # Plot primitives take coordinates on the right ("plotp(a=1,1) = 5,3,8,3;")
    # and the instrumental-variable regressions take dependent, regressor and
    # instrument.
    "%s() = 5,3,8,3;",
    "%s(a=1,1) = 5,3,8,3;",
    "%s() = 10,10;",
    "%s() = Y,X,X;",
    "%s(tp = 0 (10) 100);",
    "%s() = T,D;",
    # dltb/dple take the duration and status as two variables and the output
    # file through df=, not on the right; diple wants four.
    "%s(df=out.txt) = T,D;",
    "%s(df=out.txt) = T,D,X,Y;",
    "%s(df=out.txt) = X;",
    # Shapes the commands' own diagnostics ask for: four variables for the
    # interval regressions, an interval pair for the censored ones, a set of
    # evaluation points for the non-parametric ones, a grouping variable.
    "%s() = Y,X,G,ID;",
    "%s(yl=Y,cen=D) = Y,X;",
    "%s(x=0 (10) 100, df=out.txt) = Y,X;",
    "%s(x=0 (10) 100) = Y,X;",
    "%s(grp=G, df=out.txt) = X;",
    "%s(cen=D, df=out.txt) = Y,X;",
    "%s() = YL,YU,XL,XU;",
    "%s(df=out.txt) = YL,YU,XL,XU;",
    # range, xf and int take a function of a formal argument rather than a
    # variable list; fmin/gmin want the function through fn= and a parameter
    # on the right.
    "%s() = x*x;",
    "%s(x=0 (1) 10) = x*x;",
    "%s(fn = b1*b1) = b1;",
    # Shapes taken from the syntax boxes in doc/tman1.pdf.  int wants its
    # interval through ab=, not x=, which is why the guessed forms above all
    # failed; the plot primitives place themselves with xy=; spl and xreg take
    # a smoothing factor; segr needs the grouping variable through g=.
    "%s(ab=0,10) = x*x;",
    "%s(sig=1, df=out.txt) = X,Y;",
    "%s(g=G, df=out.txt) = X,Y;",
    "%s(xy=5,5) = 3,0.3;",
    "%s(xy=5,5) = label;",
    "%s(xy=5,5, lt=1) = 3,3;",
    "%s(x=0 (1) 10, nn=10,10) = out.txt;",
    "%s(gn=1, gt=2) = out.txt;",
    "%s(lt=1, lw=0.2) = x*x;",
    "%s(cn=1,2) = out.txt;",
    "%s(nlev=2) = out.txt;",
    # Shapes found by reading each command's argument checks next to a working
    # sibling: ilsreg takes three variables plus the censoring pair, inpreg
    # four plus evaluation points, diple a discrete interval episode.
    "%s(yl=Y,cen=D) = YL,YU,X;",
    "%s(x=1,2,3, df=out.txt) = YL,YU,XL,XU;",
    "%s(df=out.txt) = TS,TL,TU,D;",
    # sdgen from tda.hlp, plotr from doc/history.txt.
    "%s(id=ID, xyv=X,Y) = out.txt;",
    "%s(gs=1,1) = A;",
    "%s(alg=1) = out.txt;",
    "%s(M,M,M);",
    "%s(L1 = X,Y);",
    # Block constructs have to be closed in the same command file.
    "%s (0);\nendwhile;",
    "%s (n = 1);\nendrepeat;",
    "%s(M1) = { mem; };",
    # More shapes from the syntax boxes: several commands take no right-hand
    # side at all (ndvar, pltree, psetup), ploto takes a radius, plotc and
    # xconh take a function, and the matrix operators take matrix names.
    "%s(pxlen=80, pxa=0,100, pya=0,100);",
    "%s(sig=0.5);",
    "%s(xy=5,5) = 3;",
    "%s(gn=1);",
    "%s(x=1,2,3) = x*y;",
    "%s(sm=[1,1,1]) = X,Y;",
    "%s(Y1 = dum(G),);",
    "%s(A,B,C);",
    "%s(A) = out.txt;",
    # spl interpolates over a range given by rx=; smd takes its sequence of
    # smoothing operations in brackets.
    "%s(sig=1, max=50, rx=0(1)10, df=out.txt) = X,Y;",
    "%s(id=ID, xyv=X,Y) = out.txt;",
    "%s rate;",
    "%s() = electric.sav;",
    "%s() = electric.por;",
    "%s() = electric.dta;",
    "%s() = electric.dbf;",
    "%s() = deaths.xls;",
    "%s() = nc;",
    "%s() = sample.e00;",
    # sdcpol builds polygons from line segments: four variables are the two
    # endpoints, and max= has to clear the segment count.
    "%s(df=out.txt, max=200) = X,Y,YL,YU;",
    "%s() = P1,P2,P3;",
    "%s(nlev=2) = out.txt;",
    # Shapes found by reading each command's argument checks next to a working
    # sibling: ilsreg takes three variables plus the censoring pair, inpreg
    # four plus evaluation points, diple a discrete interval episode.
    "%s(yl=Y,cen=D) = YL,YU,X;",
    "%s(x=1,2,3, df=out.txt) = YL,YU,XL,XU;",
    "%s(df=out.txt) = TS,TL,TU,D;",
    # sdgen from tda.hlp, plotr from doc/history.txt.
    "%s(id=ID, xyv=X,Y) = out.txt;",
    "%s(gs=1,1) = A;",
    "%s;",
    # Shapes read out of the sources for the commands that appear in neither
    # the manual nor tda.hlp.
    "%s(yl=Y,cen=D) = YL,YU,XL,XU;",
    "%s(xv=X,Y) = X,Y;",
    "%s(yw=D) = Y,X;",
    "%s() = YL,YU;",
    "%s(rx=0(1)10, ry=0(1)10, df=out.txt) = X,Y,T;",
    "%s(sm=[3], df=out.txt) = X;",
    "%s(sm=[3,3], df=out.txt) = X,Y;",
    # 3-D shapes.  Chapter 4.7 of the manual is a stub -- four sections of
    # identical placeholder text and no syntax boxes -- so these come from
    # reading t_plot3.c and t_psf.c.
    "%s(pxlen=80, pxa=0,100, pya=0,100, pza=0,100);",
    "%s = X,Y,T;",
    "%s() = 5,5,5;",
    # The 3-D primitives place themselves with xyz=, not xy=, and plotp3
    # needs at least two points to draw between.
    "%s(xyz=5,5,5) = label;",
    "%s(xyz=5,5,5) = 3;",
    "%s() = 5,5,5, 10,10,10;",
    "%s() = x*y;",
    "%s() = 5,5;",
]

# Any diagnostic at all.  A command that always errors has not really been
# exercised, however plausible the message -- so it is reported as blocked
# rather than quietly counted as working.
ERR = re.compile(r"^(?:Error|Warning|Syntax error|Unknown command)[^\n]*",
                 re.I | re.M)

# Errors about the command's own arguments, which no change of fixture can
# fix, so the form can be dropped after one attempt.  Everything else -- and in
# particular "undefined variables", since each fixture defines different names
# -- is retried under the other fixtures.
FORM_ERR = re.compile(
    r"need right-hand side|need an output file|Unknown command|"
    r"need (?:two|three|four|exactly one|at least one)\b|"
    r"function should contain|exactly one function argument", re.I)


def tail_after(out, cmd):
    """Only the output the command itself produced.

    The preludes print too, and an error from a fixture would otherwise be
    reported as the command's.
    """
    lines = out.split("\n")
    idx = [k for k, l in enumerate(lines) if re.match(r"%s\s*[(;]" % re.escape(cmd), l)]
    return "\n".join(lines[idx[-1]:]) if idx else out


def run(text):
    cf = os.path.join(work, "s.cf")
    Path(cf).write_text(text)
    try:
        # 60s, not 25: nmca takes ~28s on a nine-edge graph, and gnst and ghd1
        # are combinatorial (see the fixture sizes above).  A genuine hang is
        # still caught, it just takes longer to declare.
        r = subprocess.run([str(binary), "cf=" + cf], cwd=work, timeout=60,
                           stdin=subprocess.DEVNULL,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.TimeoutExpired:
        return None, ""
    return r.returncode, r.stdout.decode(errors="replace")


# The scatterplot state these commands decorate lives in t_xplot.c, the X11
# module.  The portable build sets S_XWIN=0 and does not link t_xwin.c, so
# there is no window for xplot to open and no current scatterplot for the rest
# to attach to.  They are reported separately rather than counted as failures.
XWIN_ONLY = {"xplot", "xlog", "xlog1", "xconh", "xreg", "xplotf",
             "xopen", "xdelete", "xdens"}

# Only meaningful as the closing half of a block the harness cannot open with a
# single command, so they can never be reached standalone.
BLOCK_ONLY = {"endwhile", "endrepeat"}

# Every form is tried against every command, and for most commands a name on
# the right-hand side means "write your output here".  So a command earlier in
# the alphabet can and does overwrite the fixture files that later commands
# read.  Snapshot them once and restore before each command.
FIXTURES = {p.name: p.read_bytes() for p in Path(work).iterdir() if p.is_file()}


def restore_fixtures():
    for name, data in FIXTURES.items():
        Path(work, name).write_bytes(data)


# Trying every form under every prelude is ~15 x 45 attempts per command and
# almost all of them fail.  Since a command stops at its first success,
# ordering the combinations by how often they have already worked turns a full
# pass from twenty minutes into a couple of minutes.  The order is learned
# during the run, so it needs no maintenance as forms are added.
hits = {}


def combos():
    pairs = [(pn, pre, f) for pn, pre in PRELUDES for f in FORMS]
    return sorted(pairs, key=lambda t: -hits.get((t[0], t[2]), 0))


PRELUDE = dict(PRELUDES)


def _sweep():
    """The sweep itself.  Kept in a function so make_tests.py can import this
    module for its fixtures without running anything."""
    crashed, hung, clean, blocked, xwin, ext_missing = [], [], [], {}, [], []
    for c in CMDS:
        if only and c not in only:
            continue
        restore_fixtures()
        ok, msgs, plausible = False, [], None
        tries = 0

        # Forms are the outer loop, fixtures the inner one, so that a form which
        # fails for a reason the fixtures cannot change -- a syntax error, the
        # wrong number of arguments -- is dropped after one attempt instead of
        # being retried against every one.  Whatever worked last time is tried
        # first.
        combos = [(p, f) for f in FORMS for p in PRELUDES]
        hit = known.get(c)
        if hit and not FULL:
            combos.sort(key=lambda pf: (pf[0][0], pf[1]) != tuple(hit))

        for (pname, pre), form in combos:
            tries += 1
            rc, out = run(pre + (form % c) + "\n")
            if rc is None:
                hung.append("%s (%s, %r)" % (c, pname, form))
                continue
            if rc >= 128:
                crashed.append("%s (rc=%d, %s, %r)" % (c, rc, pname, form))
                continue
            m = ERR.search(tail_after(out, c))
            if not m:
                ok = True
                known[c] = [pname, form]
                if not FULL:
                    break
                continue
            msg = m.group(0).strip()[:56]
            msgs.append(msg)
            if pname == "all" and form == "%s() = out.txt;":
                plausible = msg

        if ok:
            clean.append(c)
        elif c in XWIN_ONLY or c in BLOCK_ONLY:
            xwin.append(c)
        elif not HAVE_EXT and c in EXT_ONLY:
            ext_missing.append(c)
        else:
            blocked[c] = (plausible or
                          (max(set(msgs), key=msgs.count) if msgs else "?")) + \
                         ""

    print("commands tried      : %d" % (len(CMDS) if not only else len(only)))
    print("  ran clean at least once : %d" % len(clean))
    print("  always errored          : %d" % len(blocked))
    print("  not reachable alone     : %d" % len(xwin))
    if ext_missing:
        print("  need external input     : %d  (set TDA_EXT_INPUT)" % len(ext_missing))
    print("  crashed                 : %d" % len(crashed))
    print("  hung                    : %d" % len(hung))
    for c in crashed:
        print("CRASH: " + c)
    for c in hung:
        print("HANG:  " + c)
    if VERBOSE and blocked:
        print("\nalways errored -- needs a fixture this harness does not build:")
        for c in sorted(blocked):
            print("  %-12s %s" % (c, blocked[c]))

    try:
        CACHE.write_text(json.dumps(known, indent=0, sort_keys=True))
    except Exception:
        pass

    return 1 if crashed or hung else 0


if __name__ == "__main__":
    sys.exit(_sweep())
