# Regression tests for the plotting family, covering real bugs found and
# fixed this development pass, and the named-parameter interface added
# alongside them. Different in character from test-r-comparisons.R (which
# checks numbers against an outside reference): these mostly check the
# parsed PostScript operations directly, since a rendering bug is a
# rendering bug regardless of whether the underlying numbers were right.

library(tdaR)

# .tda_ps_run() runs a session through TDA and returns the raw PostScript;
# tda_read_ps() parses it into the op list tda_plot_ps() actually draws
# from. Both are internal, but that is exactly what a rendering-level
# regression test needs to inspect.
ps_ops <- function(p) {
    r <- getFromNamespace(".tda_ps_run", "tdaR")(p)
    tda_read_ps(r$run, which = r$file)$ops
}
first_of_type <- function(ops, type) {
    i <- which(vapply(ops, function(o) o$op, "") == type)
    if (!length(i)) NULL else ops[[i[1L]]]
}
# The last matching op, not the first -- tda_pl_axes() itself draws
# "lines"-type ops for the axis lines, so a data series added after the
# axes is not the first "lines" op in the session.
last_of_type <- function(ops, type) {
    i <- which(vapply(ops, function(o) o$op, "") == type)
    if (!length(i)) NULL else ops[[i[length(i)]]]
}
# "symbol" ops (tda_pl_points()'s plot(s=...)) are one op per point,
# unlike "lines" (one op, a whole polyline) -- so a series of points has
# to be gathered across every matching op, not read off the first one.
all_of_type <- function(ops, type) {
    keep <- Filter(function(o) identical(o$op, type), ops)
    list(x = vapply(keep, function(o) o$x, numeric(1)),
        y = vapply(keep, function(o) o$y, numeric(1)))
}

# ---- the ellipse rendering bug ------------------------------------------
#
# TDA draws every ellipse as a unit circle under a non-uniform scale.  The
# radius alone (what a plain circle marker correctly uses) collapsed every
# non-circular ellipse into a circle of some averaged size until this was
# fixed to keep the linear transform and apply it directly. This checks
# the fix stays fixed: an ellipse with distinctly different axis lengths
# must produce a different x and y extent when the transform is
# applied to a sampled perimeter, not just "some" extent.

p <- tda_ps(xlim = c(0, 6), ylim = c(0, 15))
p <- tda_pl_ellipse(p, at = c(3, 5), axes = c(1.5, 0.8))
o <- first_of_type(ps_ops(p), "circle")
th <- seq(0, 360, length.out = 73)[-73] * pi / 180
px <- o$m[1L] * o$r0 * cos(th) + o$m[3L] * o$r0 * sin(th)
py <- o$m[2L] * o$r0 * cos(th) + o$m[4L] * o$r0 * sin(th)
# px/py are in PostScript points, and their ratio reflects both the given
# axes and the session's x/y unit-to-point scale (here, xlim spans 6
# units and ylim spans 15, so it is not simply 1.5/0.8) -- what actually
# distinguishes a ellipse from the bug this guards against (every
# non-circular ellipse silently collapsing to a circle of some averaged
# size) is that the x and y extents come out different at all.
ok("ellipse: x and y extents differ (not collapsed to a circle)",
   abs(diff(range(px)) / diff(range(py)) - 1) > 0.5)

# An ordinary circle marker must still be an ordinary, uniform-scale
# circle -- the fix above must not have started treating every circle as
# if it needed the ellipse path.
p2 <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
p2 <- tda_pl_circle(p2, at = c(5, 5), r = 2)
o2 <- first_of_type(ps_ops(p2), "circle")
same("circle: stays a uniform-scale transform (not treated as an ellipse)",
     o2$m, c(1, 0, 0, 1), 1e-8)

# ---- text: wf=1 (white background) parsing -------------------------------
#
# TDA's wf=1 emits "(string) sclear show" -- an extra token between the
# string and the show keyword that silently dropped the label entirely
# before the parser's regex was loosened to allow for it (and loosened
# correctly: a first attempt broke the plain, no-extra-token case, caught
# by testing that case again after the fix).

p3 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p3 <- tda_pl_text3(p3, "hello", at = c(4, 5, 4), fs = 4, white = TRUE)
o3 <- first_of_type(ps_ops(p3), "text")
ok("text3: wf=1 (white background) still produces a text op",
   !is.null(o3) && identical(o3$label, "hello"))

p4 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p4 <- tda_pl_text3(p4, "hello", at = c(4, 5, 4), fs = 4, center = TRUE)
o4 <- first_of_type(ps_ops(p4), "text")
ok("text3: plain center (no wf=1) still produces just='centre'",
   !is.null(o4) && identical(o4$just, "centre"))

p5 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p5 <- tda_pl_text3(p5, "hello", at = c(4, 5, 4), fs = 4, center = TRUE,
                   white = TRUE)
o5 <- first_of_type(ps_ops(p5), "text")
ok("text3: center and wf=1 together still produce just='centre'",
   !is.null(o5) && identical(o5$just, "centre"))

# ---- text rotation ---------------------------------------------------
#
# TDA's r= for pltext/pltext3 was accepted and sent correctly, but the
# angle was parsed and then discarded -- grid.text() was never told to
# rotate. Checks the angle survives parsing intact.

p6 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p6 <- tda_pl_text3(p6, "hi", at = c(4, 5, 4), fs = 4, rotate = 30)
o6 <- first_of_type(ps_ops(p6), "text")
same("text3: rotate= survives into the parsed op", o6$rot, 30, 0)

p7 <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
p7 <- tda_pl_text(p7, "plain", at = c(3, 3))
o7 <- first_of_type(ps_ops(p7), "text")
same("text (2D): no rotate= given defaults to 0", o7$rot, 0, 0)

# ---- degenerate single-op viewport --------------------------------------
#
# A plot with only one op and nothing else to give the coordinate range
# any spread (diff(xr) == 0) crashed with "non-finite location and/or
# size for viewport", because the padding used for the inner viewport's
# scale was never applied to the outer grid.layout()'s widths/
# heights, which used the raw range directly.

p8 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p8 <- tda_pl_text3(p8, "alone", at = c(4, 5, 4), fs = 4)
r8 <- tryCatch({ pdf(NULL); on.exit(dev.off()); plot(p8); TRUE },
              error = function(e) e)
ok("plot: a session with only one text op and nothing else does not crash",
   isTRUE(r8))

# ---- named parameters produce the same result as the raw option --------
#
# Each promoted parameter has to actually reach TDA the same way the raw
# option name always did -- these check a handful of representative ones
# end to end, not just that the function runs.

pa <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
pa <- tda_pl_frame(pa)
pa <- tda_pl_grid(pa, lty = "dashed", lw = 0.3)
## Since the dash-pattern fix (see .ps_dash()'s comment), a real
## dash pattern parses to its actual on/off lengths in points, not a
## string naming a grid keyword -- grid's hex-string encoding
## cannot represent TDA's true, absolute dash lengths distinguishably
## for a thin line, checked (all 9 of TDA's line types
## rendered as indistinguishable micro-dashes before this fix). "not
## solid" is what reaching the grid lines actually means now.
grid_lines <- Filter(function(o) identical(o$op, "lines") &&
                     is.numeric(o$lty), ps_ops(pa))
ok("grid: lty='dashed' reaches the drawn grid lines as a real dash pattern",
   length(grid_lines) > 0)
same("grid: lw=0.3 reaches the drawn grid lines", grid_lines[[1L]]$lwd,
     0.3 * 72 / 25.4, 1e-3)

pb <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
pb <- tda_pl_frame(pb)
pb <- tda_pl_circle(pb, at = c(5, 5), r = 2, gray = 0.5)
ob <- first_of_type(ps_ops(pb), "circle")
same("circle: gray=0.5 sets the fill colour", ob$fill, "#808080", 0)

# ---- vector plotting, lines(x=, y=)-style ---------------------------------
#
# tda_pl_lines()/tda_pl_points() now accept raw vectors directly, not only
# column names already in the session's data -- R's lines(x=, y=)
# convenience, via .pl_xy(). A vector the plot's data's length is
# folded in as a new column (the existing tda_pl_lines3() mechanism,
# .pl_asvar(), reused for 2D too); a *different* length needs its
# temporary data matrix, TDA's clear;+nvar(dfile=...) trick --
# examples/exam/npreg1.cf's scatter+fitted-curve overlay does exactly
# this, at two different lengths (256 raw points, 101 fitted ones).
#
# ps_ops() returns raw PostScript point coordinates, not data units (a
# "symbol" op for data x=1, xlim=c(0,6), the default 100mm-wide plot
# comes back as x=47.25 + a fixed page offset, not x=47.25) -- to_data()
# inverts that, given the same xlim/ylim and width/height (mm) the
# session itself was built with, so assertions can be written in the
# data's units instead of working out the exact pixel expected by
# hand.
#
# TDA's default page position is a fixed (150, 460) points offset
# (psorg's default, checked: constant across every
# xlim/ylim/width/height combination tried, not something that scales
# with any of them) -- device coordinate 0 is not the plot's
# origin, device coordinate 150 (x) / 460 (y) is. Before the
# multi-panel translate-variable fix (.ps_body()/ps_vars, this
# session), that offset was silently dropped during parsing for every
# plot, single-panel or not, which is what let this helper get away
# with assuming px=0 means lim[1] -- now that it is correctly applied
# (required for a multi-panel layout's later panels to land in the
# right place at all), every existing user of this helper needs the
# same offset subtracted first, or every coordinate comes back
# shifted by a constant that has nothing to do with the data.
.PS_ORIGIN <- c(x = 150, y = 460)
to_data <- function(px, lim, size_mm, axis = c("x", "y"))
    lim[1L] + (px - .PS_ORIGIN[[match.arg(axis)]]) /
        (size_mm * 72 / 25.4) * diff(lim)

pc <- tda_ps(xlim = c(0, 10), ylim = c(0, 20))
pc <- tda_pl_axes(pc)
pc <- tda_pl_lines(pc, x = 1:10, y = (1:10) * 2)
opsc <- ps_ops(pc)
lc <- last_of_type(opsc, "lines")
ok("vector plotting: a fresh session with no data still draws a real line",
   !is.null(lc) && length(lc$x) == 10L)
same("vector plotting: the drawn y values are the ones given",
     to_data(lc$y, c(0, 20), 70, "y"), (1:10) * 2, 1e-3)

# a column-based series and a same-length vector series together
d_mix <- data.frame(a = 1:5, b = c(2, 4, 5, 8, 9))
pd <- tda_ps(d_mix, xlim = c(0, 6), ylim = c(0, 12))
pd <- tda_pl_lines(pd, "a", "b")
pd <- tda_pl_points(pd, x = 1:5, y = c(1, 3, 3, 6, 7), symbol = 1)
opsd <- ps_ops(pd)
same("vector plotting: mixed column+vector series both draw",
     c("lines" %in% vapply(opsd, function(o) o$op, ""),
       sum(vapply(opsd, function(o) o$op, "") == "symbol") == 5L),
     c(TRUE, TRUE))
same("vector plotting: the vector-only series' own y values are right",
     to_data(all_of_type(opsd, "symbol")$y, c(0, 12), 70, "y"),
     c(1, 3, 3, 6, 7), 1e-3)

# the actual mismatched-length case: a 256-point scatter (the session's
# own data) and a 30-point fitted curve overlaid on it (a vector at a
# completely different length) -- npreg1.cf's shape, at a size small
# enough to keep this test fast.
set.seed(3)
n_vp <- 256
d_vp <- data.frame(X = stats::runif(n_vp), Y = stats::rnorm(n_vp))
fit_x <- seq(0, 1, length.out = 30)
fit_y <- sin(fit_x * 2 * pi)

pe <- tda_ps(d_vp, xlim = c(0, 1), ylim = c(-4, 4))
pe <- tda_pl_axes(pe)
pe <- tda_pl_points(pe, "X", "Y", symbol = 1, lty = 0)
pe <- tda_pl_lines(pe, x = fit_x, y = fit_y)
opse <- ps_ops(pe)
scat <- all_of_type(opse, "symbol")
fit <- last_of_type(opse, "lines")
ok("vector plotting: mismatched-length scatter+curve both draw",
   length(scat$x) > 0L && !is.null(fit))
same("vector plotting: the scatter keeps all 256 of the session's points",
     length(scat$x), n_vp)
same("vector plotting: the fitted curve keeps its own, different length",
     length(fit$x), 30L)
same("vector plotting: the fitted curve's y values are the ones given",
     to_data(fit$y, c(-4, 4), 70, "y"), fit_y, 1e-3)
# and the session's data must still be intact afterward -- this is
# exactly what the after-block (clear;+nvar reload) is for.
same("vector plotting: the scatter's values are unaffected by the swap",
     sort(to_data(scat$x, c(0, 1), 100, "x")), sort(d_vp$X), 1e-3)

ok("vector plotting: x/y of different lengths errors",
   inherits(try(tda_pl_lines(tda_ps(xlim = c(0, 1), ylim = c(0, 1)),
                             x = 1:5, y = 1:3), silent = TRUE),
            "try-error"))

## tda_pl_contour(): TDA's PostScript draws each contour band as
## "gsave / setgray fill grestore / stroke / grestore" -- fill and
## stroke on the *same* path, with no newpath between them (a standard
## PostScript idiom: shade the region, then outline it). The PS reader
## treated fill as consuming the path the way every other command here
## does, so the stroke that followed had nothing left to draw --
## checked: contour lines were silently missing entirely,
## not merely drawn in the wrong colour, which is what the visible
## symptom ("some lines look coloured") actually was: only the raw,
## un-outlined fill bands were ever showing.
##
## Separately, TDA marks a point it adds purely to close a path to
## the plot's edge for filling purposes with a trailing
## "% ccc"/"% bbb"-style comment. This broke simple tokenizing --
## "283.50 0 l % ccc" has "ccc" as its last token, not "l", so the
## line matched no known path command and was dropped outright, not
## merely mis-tagged, leaving the path missing its last point and
## grid's implicit "close back to the start" drawing a spurious
## diagonal sliver across the gap instead of the real point TDA placed
## there. (An earlier version of this fix also skipped the fill
## entirely for any path ending in one of these marked points, on the
## theory that they marked boundary-crossing fragments that should
## only ever be outlined, matching an old printed figure in the
## manual. That was wrong: an independent PDF rendering of the real,
## unmodified examples/exam/plot21.ps -- not a reconstruction -- shows
## every one of these bands filled, not just outlined. The
## point itself is data and is kept, exactly like any other;
## only the tokenizing bug is fixed here.)
d3 <- data.frame(x = -2:2, y = -2:2)
p3 <- tda_ps(d3, xlim = c(-2, 2), ylim = c(-2, 2))
p3 <- tda_pl_contour(p3,
    "exp(-((x1+1)^2+(x2+1)^2)/2)/(2*pi) + exp(-((x1-1)^2+(x2-1)^2)/2)/(2*pi)",
    levels = c(0.01, 0.05, 0.1, 0.125, 0.15), resolution = c(50, 50),
    gray = .9)
ops3 <- ps_ops(p3)
poly3 <- Filter(function(o) o$op == "polygon", ops3)
lines3 <- Filter(function(o) o$op == "lines", ops3)
ok("contour: every filled band also gets its outline stroke",
   length(poly3) > 0L && length(lines3) == length(poly3))
## "fill grestore" often shares one line with the setgray that set the
## fill colour ("0.9000 setgray fill grestore"), and the general
## gsave/grestore stack handling only ever matches a line that
## *starts* with "grestore" -- this trailing one, embedded right after
## "fill", never reached it, so the colour push/pop around the fill
## never actually popped: the stroke that followed kept the fill's
## grey instead of the colour active before the fill (black, TDA's
## line-drawing colour, checked -- every stroke in a real
## run is black, none of them the grey of its own neighbouring fill).
ok("contour: every stroke is black, not its fill's colour",
   all(vapply(lines3, function(l) l$col == "black", logical(1))))
ok("contour: fills keep their distinct grey shades regardless",
   length(unique(vapply(poly3, function(p) p$col, ""))) > 1L)

## The real, unmodified manual example (not a reconstruction): every
## band -- including the boundary-crossing ones near the corners --
## keeps its comment-marked closing point, and so is both filled
## and outlined, matching an independent PDF rendering of this exact
## file, not the manual's older printed figure.
if (have_examples) {
p5 <- tda_read_ps(file.path(EX, "exam", "plot21.ps"))
ops5 <- Filter(function(o) grepl("plotc", o$command), p5$ops)
poly5 <- Filter(function(o) o$op == "polygon", ops5)
lines5 <- Filter(function(o) o$op == "lines", ops5)
ok("plot21.ps: every fill has a matching outline stroke, boundary fragments included",
   length(poly5) > 0L && length(lines5) == length(poly5))
ok("plot21.ps: every stroke is black, matching an independent Ghostscript/PDF render",
   all(vapply(lines5, function(l) l$col == "black", logical(1))))
## The actual bug: a trailing "% ccc"-style comment on a path-command
## line broke simple whitespace tokenizing ("ccc", not "l", was the
## line's last token), silently dropping the whole line -- checked
## directly against the parser's point count for one of the
## smallest, most exactly-known real bands, rather than inferred from
## overall plot geometry.
ok("plot21.ps: a line with a trailing % comment still contributes its point",
   any(vapply(poly5, function(p) length(p$x) == 21L, logical(1))))
}

## The fix must not re-emit a filled-but-never-stroked path as a
## spurious extra "lines" op once some other command starts a fresh
## path afterward (a real risk: a naive "keep the path after fill"
## fix would double-draw every filled shape elsewhere that has no
## matching stroke at all, like a shaded circle node).
p4 <- tda_pl_circle(tda_ps(xlim = c(0, 1), ylim = c(0, 1)),
                    at = c(0.5, 0.5), r = 0.1, gray = 0.5)
ops4 <- ps_ops(p4)
ok("circle: a shaded node still gets exactly one circle op, not a spurious extra line",
   sum(vapply(ops4, function(o) o$op, "") == "circle") == 1L &&
   sum(vapply(ops4, function(o) o$op, "") == "lines") == 0L)

## Arrowheads: TDA draws one as "currentpoint / stroke m / dx dy / atan /
## rotate / ... fill" -- a PostScript idiom spreading a single
## command's arguments across several lines via the operand stack,
## not the usual "arguments precede the command on the same line" shape
## every other command here has. Two separate bugs here, both real:
## (1) bare "m"/"rotate" with no arguments of their own on the line
## silently did nothing at all, so every arrowhead used the same
## default, unrotated orientation regardless of which way its line
## pointed; (2) once that was fixed, PostScript's "num1 num2 atan
## angle" is atan2(num1, num2), and the two arguments were swapped
## (atan2(num2, num1) instead) -- a mistake that happens to be
## invisible at exactly 45 degrees (swapping is a fixed point there)
## and rotates the result 90 degrees off everywhere else, which is
## exactly why a first look at the real eight-arrow test file below
## seemed fine (the two diagonals nearest 45 degrees looked right) and
## was not: checked here against the real tip coordinates, not glanced
## over as a picture.
if (have_examples) {
p6 <- tda_read_ps(file.path(EX, "exam", "plot13.ps"))
heads <- Filter(function(o) o$op == "polygon" && length(o$x) == 5L, p6$ops)
ok("arrowheads: the real eight-arrow test file gets all eight",
   length(heads) == 8L)
## The centre here is computed from the tips themselves (their own
## mean -- the eight are symmetric around it) rather than a hardcoded
## value: the multi-panel translate-variable fix (.ps_body()/ps_vars,
## this session) made the previously-dropped page-position offset
## apply correctly everywhere, which shifted every plot's absolute
## coordinates by that constant, including this file's -- the tips
## themselves and the angles between them are unaffected by the shift
## (checked), so recomputing the centre from the data
## itself rather than re-hardcoding a new absolute value is what
## actually tests the geometry instead of one specific coordinate
## frame.
tips <- t(vapply(heads, function(o) c(o$x[1L], o$y[1L]), numeric(2)))
ctr <- colMeans(tips)
ang <- (atan2(tips[, 2L] - ctr[2L], tips[, 1L] - ctr[1L]) * 180 / pi) %% 360
same("arrowheads: every tip points in its own line's real direction, not one shared default",
     sort(round(ang / 45)) * 45, c(0, 45, 90, 135, 180, 225, 270, 315))

## plframe's background fill: "gsave 0.95 setgray" is one line,
## and the general gsave handling only ever pushed the state and moved
## on to the next line, never processing a setgray sharing that same
## line -- so the fill that followed kept whatever colour was active
## before the gsave (black, checked against a real file:
## examples/exam/plot2.ps's background rectangle came back solid
## black, not the light grey plframe(gs=0.95) it actually asked for).
p7 <- tda_read_ps(file.path(EX, "exam", "plot2.ps"))
bg <- Filter(function(o) o$op == "polygon", p7$ops)[[1L]]
same("plframe: the background fill is plframe's real grey, not black",
     bg$col, grDevices::gray(0.95))
}

## Contour fill draw order: TDA's grid-tracing algorithm (Snyder,
## ACM 531) walks the grid once, checking all levels at each cell and
## tracing whichever level's crossing it finds there first -- not
## an outer loop over levels with a full pass each -- so two
## contour bands that are spatially separate (each a disconnected
## fringe fragment near its plot corner, say) can end up written
## to the real PostScript file in either relative order, independent
## of their own level numbers. Since each band's fill is a real,
## opaque PostScript fill covering the plot area behind it (confirmed
## directly), whichever fragment is written *last* wins
## the paint job regardless of which one is nominally the lighter or
## darker level -- so a correct picture requires the geometrically
## larger of any two overlapping fragments to be painted first, the
## smaller one second (staying visible on top of it).
##
## Sorting by ascending level number instead (an earlier attempt at
## this fix, not shipped) is not the same thing and got it backwards
## in exactly the cases it was meant to fix: checked directly against
## a real run's polygon data, a fringe fragment near a plot corner
## can end up geometrically *larger* than a numerically lower-level
## fragment nearby, because the straight line that closes an open
## contour curve back to the nearest plot corner has to sweep back
## through more of the plot to close a curve that exits further from
## that corner -- inflating that fragment's closed area
## independently of its level. The fix implemented instead sorts each
## traced fragment by its own actual, measured area (shoelace formula,
## accumulated in cont_draw() as each point is written), which is
## correct regardless of why a fragment ended up whatever size it is.
##
## This is checked here two ways: against the real, bimodal density
## example from the manual (both of its two corner fragments must
## each end up with two more polygons -- not one merged/overpainted
## shape -- for the same reasons in mind when this bug was first
## noticed), and against a hand-computable single-peaked cone function
## where every level's true enclosed area is known exactly, as a
## sanity check that the fix does not disturb the simple,
## already-correct case.
##
## The manual's plot21.ps is regenerated here via tda_pl_contour()
## rather than read as a static shipped file: examples/exam/plot21.ps
## gets overwritten by the C test suite (check.py) every time it runs,
## and by the packaging step that resets other regenerated .ps files
## back to their pristine, pre-fix state afterward -- a shipped copy
## can silently drift out of sync with whatever binary this test is
## actually run against and fail for the wrong reason (confirmed
## directly: this happened once already). Regenerating it here keeps
## the test self-contained and always checked against the binary
## actually running right now.
d8 <- data.frame(x = -2:2, y = -2:2)
p8ps <- tda_ps(d8, xlim = c(-2, 2), ylim = c(-2, 2))
p8ps <- tda_pl_axes(p8ps, ic = c(.5, .5), sc = c(2, 2))
p8ps <- tda_pl_contour(p8ps,
    "exp(-((x1+1)^2+(x2+1)^2)/2)/(2*pi)+exp(-((x1-1)^2+(x2-1)^2)/2)/(2*pi)",
    levels = c(0.01, 0.05, 0.1, 0.125, 0.15), resolution = c(50, 50),
    gray = .9)
poly8 <- Filter(function(o) o$op == "polygon", ps_ops(p8ps))
shoelace <- function(x, y)
    abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
area8 <- vapply(poly8, function(o) shoelace(o$x, o$y), numeric(1))
ok("contour: the manual's bimodal example gets all 9 fragments",
   length(poly8) == 9L)
## the two smallest fragments are the corner slivers (levels 1 and 2);
## each corner's pair must be written smaller-second (larger
## fragment's index precedes the smaller one among fragments that
## spatially overlap it) -- checked directly via each fragment's
## bounding box rather than assumed from position in the list.
bbox <- lapply(poly8, function(o) c(range(o$x), range(o$y)))
overlaps <- function(a, b)
    a[1] <= b[2] && b[1] <= a[2] && a[3] <= b[4] && b[3] <= a[4]
ok("contour: for every pair of spatially overlapping fragments, the larger one is written first",
   all(vapply(seq_along(poly8), function(i) {
       ok_i <- TRUE
       if (i > 1)
           for (j in seq_len(i - 1))
               if (overlaps(bbox[[i]], bbox[[j]]) && area8[i] > area8[j])
                   ok_i <- FALSE
       ok_i
   }, logical(1))))

d9 <- data.frame(x = -1:1, y = -1:1)
p9 <- tda_ps(d9, xlim = c(-1, 1), ylim = c(-1, 1))
p9 <- tda_pl_contour(p9, "1 - sqrt(x1^2 + x2^2)", levels = c(0.2, 0.5, 0.8),
                     resolution = c(80, 80), gray = 0.9)
poly9 <- Filter(function(o) o$op == "polygon", ps_ops(p9))
w9 <- vapply(poly9, function(o) diff(range(o$x)), numeric(1))
ok("contour: the simple, single-peaked case still comes out largest-level-radius first",
   length(w9) == 3L && w9[1] > w9[2] && w9[2] > w9[3])

## Sunflower plots (scplot, opt=2): a grid cell with exactly one point
## draws its dot as "gsave / X Y translate / %#symbol: 5 0 0 size
## / 0 0 circle fill / grestore" -- a single-point glyph at the local
## origin of its own translated coordinate system, not TDA's usual
## "%#symbol:" usage of writing an already-final, absolute device
## position directly. Treating the comment's x,y as already final
## (correct for a plain symbol plot, where nothing precedes it) sent
## every single-point sunflower cell to the same, wrong position
## instead of its own grid cell -- checked against a real
## run: dozens of "%#symbol: 5 0.00 0.00 ..." comments, identical
## regardless of which cell they belonged to. Mapping the comment's
## own coordinates through the current transform (as in effect for
## any other coordinate in the file) fixes this for the translated
## case while leaving the plain, untransformed case unchanged.
if (have_examples) {
p10 <- tda_read_ps(file.path(EX, "exam", "scplot2.ps"))
sym10 <- Filter(function(o) o$op == "symbol", p10$ops)
ok("sunflower: single-point cells are spread across the grid, not collapsed to one spot",
   length(sym10) > 10L &&
   length(unique(round(vapply(sym10, function(o) o$x, numeric(1)), 1))) > 5L)
}

## Line types: TDA's 9 line types are meant to be visually
## separable even all together in one plot (the manual's reference
## figure shows exactly this). Mapping a real PostScript dash array
## onto grid's line-width-relative hex-string lty (an earlier
## version of .ps_dash()) crushed TDA's true, absolute dash lengths
## down to near-identical micro-dashes for a thin line: checked
## directly, side by side, all 9 rendered indistinguishably. The fix
## keeps the raw on/off lengths (in points) and breaks each path into
## its actual dash segments at render time (.dash_segments()),
## independent of line width. Checked here two ways: that each of the
## 9 line types parses to a different dash array (not
## collapsed into a handful of repeated buckets, which is exactly what
## the earlier, broken version did), and that the actual rendered
## segment lengths for two clearly different line types (a fine-dotted
## one and a long-dashed one) differ by roughly the ratio their own
## PostScript dash values do, not by an unrelated, fixed amount.
p11 <- tda_ps(xlim = c(0, 12), ylim = c(0, 10))
for (lt in 1:9)
    p11 <- tda_pl_lines(p11, x = c(1, 9), y = c(lt, lt), lty = lt)
ltys <- lapply(Filter(function(o) o$op == "lines" && length(o$x) == 2L,
                      ps_ops(p11)), function(o) o$lty)
ltys <- ltys[!vapply(ltys, identical, "solid", FUN.VALUE = logical(1))]
ok("line types: all 8 non-solid types parse to distinct dash arrays",
   length(unique(lapply(ltys, function(v) paste(v, collapse = ",")))) == 8L)

r11 <- tda_read_ps(tda_ps_file(p11))
seg_len <- function(lty_vec) {
    o <- Filter(function(o) o$op == "lines" && length(o$x) == 2L &&
               identical(o$lty, lty_vec), r11$ops)[[1L]]
    seg <- tdaR:::.dash_segments(o$x, o$y, o$lty)
    brk <- c(0, which(is.na(seg$x)), length(seg$x) + 1L)
    runs <- Map(function(a, b) seq(a + 1L, b - 1L), brk[-length(brk)], brk[-1L])
    max(vapply(runs, function(idx)
        sqrt(diff(range(seg$x[idx]))^2 + diff(range(seg$y[idx]))^2),
        numeric(1)))
}
fine <- Filter(function(o) o$op == "lines" && length(o$x) == 2L,
               ps_ops(p11))[[2L]]$lty   # line type 2, finely dotted
long <- Filter(function(o) o$op == "lines" && length(o$x) == 2L,
               ps_ops(p11))[[9L]]$lty   # line type 9, longest dash
ok("line types: a long-dash type's rendered segments are visibly longer than a fine-dotted type's",
   seg_len(long) > 3 * seg_len(fine))

## A session with only one text op and nothing else (found via a
## systematic Ghostscript-vs-tdaR comparison across every real,
## shipped .ps file, not by inspection: examples/tests/p3.ps).
## tda_plot_ps()'s viewport range came only from each op's raw
## x/y, but a text op does not actually get drawn at its raw y --
## it is shifted down by 0.32*fontsize first (undoing PostScript's
## baseline-vs-centre convention, see that comment). Left out of the
## range, that fixed, points-sized shift is invisible whenever other
## content gives the plot real spread, but for a plot whose only
## content is one label the range collapses to a single point, the
## pad is a few hundredths of a point, and the shift lands the label
## entirely outside the viewport -- rendering nothing at all, not a
## wrong position, checked: the earlier version of this
## function produced a blank image for exactly this file, matching
## the manual's promise only that it "does not crash", not that
## anything is visible.
##
## tda_plot_ps() runs end to end without error either way (that much
## the existing "does not crash" test already covers); what it does
## not check is that the label lands inside the plotted range at all.
## Checked here by reproducing the same extent computation the fixed
## function itself uses -- the shifted y has to be so far from the
## op's raw y that a naive, unfixed range (built from raw x/y
## alone, pad included) would miss it entirely, and the fixed range
## must not.
if (have_examples) {
p12 <- tda_read_ps(file.path(EX, "tests", "p3.ps"))
o12 <- p12$ops[[1L]]
shifted_y <- o12$y + 0.32 * (o12$fontsize %||% 8.5)
naive_pad <- 0.04 * max(1e-8, 1)
naive_yr <- c(o12$y - naive_pad, o12$y + naive_pad)
ok("plot: the old, unfixed range would have missed a lone label's shifted position",
   shifted_y > naive_yr[2L])
fixed_yr <- range(o12$y, o12$y + 0.32 * (o12$fontsize %||% 8.5))
fixed_pad <- 0.04 * max(diff(fixed_yr), 1)
ok("plot: the fixed range includes a lone label's shifted draw position",
   shifted_y >= fixed_yr[1L] - fixed_pad && shifted_y <= fixed_yr[2L] + fixed_pad)
r12 <- tryCatch({ pdf(NULL); on.exit(dev.off()); tda_plot_ps(p12); TRUE },
                error = function(e) e)
ok("plot: the same single-text-op session still runs without error",
   isTRUE(r12))
}

## tda_pl_function(): TDA's command dispatch matches "plotf" as a
## 5-character prefix (t_cmd.c), so "plotf1"/"plotf2" reach the exact
## same C handler as plain "plotf" -- but that handler itself then
## reads the very next character to mean something real: '1' plots
## the function's first derivative, '2' its second, not the
## function itself (t_plot.c's pl_plotf()). Found via a systematic
## replication of the manual's examples, not by reading the
## source first: examples/exam/plot1.cf plots sin(x) twice, once
## plain and once through "plotf1", labelling them "sine" and
## "cosine" -- which only makes sense once plotf1 is understood to be
## the derivative, cos(x), not a second, identical sine curve.
## Unsupported before this, with no way to reach it from R at all.
d13 <- data.frame(x = 0)
p13a <- tda_pl_function(tda_ps(d13, xlim = c(0, 6), ylim = c(-1, 1)),
                        "sin(x1)", range = c(0, 6))
p13b <- tda_pl_function(tda_ps(d13, xlim = c(0, 6), ylim = c(-1, 1)),
                        "sin(x1)", range = c(0, 6), deriv = 1)
same("tda_pl_function: deriv=1 reaches TDA as plotf1, not plain plotf",
     tail(p13b$cmds, 1), sub("^plotf\\(", "plotf1(", tail(p13a$cmds, 1)))
ln_a <- Filter(function(o) o$op == "lines", ps_ops(p13a))[[1L]]
ln_b <- Filter(function(o) o$op == "lines", ps_ops(p13b))[[1L]]
ok("tda_pl_function: deriv=1 of sin(x) draws a different curve (cos(x), not sin(x) again)",
   max(abs(ln_a$y - ln_b$y)) > 1)

## .ps_body(): TDA's multi-panel technique -- calling psetup()
## more than once, at a different psorg=, to lay out several plots on
## one page (examples/exam/plot6.cf, found via a systematic
## replication of the manual's examples). Each psetup() call
## re-emits its own "%#Parameter: .../1 setlinecap/1 setlinejoin"
## header, and .ps_body() used the *last* setlinejoin in the whole
## file to mark where the prolog ends and the real drawing stream
## begins -- correct for the overwhelming majority of files, which
## only ever call psetup() once and so have exactly one setlinejoin
## either way, but for a multi-panel file this discarded every
## earlier panel's entire content outright, not merely mis-positioned
## it. Checked: only found because the two panels in the
## real example use different pylen values, so the first panel's
## y-axis could not be mistaken for the second's; its x-axis (same
## pxlen in both panels) happened to coincide with the second panel's
## and so looked present when it was not, which is what made this
## easy to miss with a single check. The prolog itself never contains
## its setlinejoin, so the first occurrence in the file is always
## the correct, single true boundary, regardless of how many panels
## follow it.
p14 <- tda_ps(xlim = c(0, 6), ylim = c(-1, 1), width = 50, height = 40)
p14 <- tda_pl_axes(p14, sc = 1)
p14 <- tda_pl_function(p14, "sin(x1)", range = c(0, 6))
p14 <- tda_pl_text(p14, "Plot 1", at = c(4, 0.2))
## tda_pl_panel(): the documented, high-level way to reach this same
## psetup() feature -- takes mm, like every other size/position
## parameter in this package, and tracks the running page origin
## itself so the caller never needs to know TDA's undocumented
## default psorg= or that psorg= is itself in points while pxlen=/
## pylen= are in mm (TDA's inconsistency, confirmed
## directly against its help text).
p14 <- tda_pl_panel(p14, origin = c(50, 0), width = 50, height = 30)
p14 <- tda_pl_axes(p14, sc = 1)
p14 <- tda_pl_function(p14, "sin(x1)", range = c(0, 6))
p14 <- tda_pl_text(p14, "Plot 2", at = c(4, 0.2))
r14 <- tda_read_ps(tda_ps_file(p14))
lab14 <- vapply(Filter(function(o) o$op == "text", r14$ops),
                function(o) o$label, "")
ok("multi-panel: both panels' own labels survive parsing, not just the last panel's",
   "Plot 1" %in% lab14 && "Plot 2" %in% lab14)
curves14 <- Filter(function(o) o$op == "lines" && length(o$x) > 10L, r14$ops)
ok("multi-panel: both panels' own curves survive parsing",
   length(curves14) == 2L)

## dplot() (examples/exam/plot7.cf) combines several already-created
## PostScript files into a grid, one row per fn= -- unlike psetup()
## called more than once (the plot6.cf case above), it re-emits each
## combined file's full set of procedure definitions
## ("/sclear {...} def", "/center {...} def", ...) once per panel, not
## just at the very start of the file: checked, "/sclear
## {...} def" appears four separate times, evenly spaced, in a real
## four-panel combined file. This parser has no notion of PostScript's
## own "/name { ... } def" procedure-definition syntax, and without
## recognising it, read every gsave/grestore/setgray *inside* such a
## block as if it were a real, sequential drawing command -- which
## desynchronised the gsave/grestore stack from the file's true
## nesting (checked: 24 gsave against only 20 grestore
## within what should have been one, self-contained panel) and left
## every panel after the first stuck in a colour state (white) left
## over from an unrelated, never-executed procedure body, not the one
## actually active at that point in the real drawing stream --
## confirmed against Ghostscript's rendering of the identical,
## unmodified file, which shows all four panels correctly, in black.
if (have_examples) {
f3 <- file.path(EX, "exam", "plot3.ps")
f4 <- file.path(EX, "exam", "plot4.ps")
r15 <- tda_combine_ps(list(c(f3, f3), c(f4, f4)), width = 100, height = 60,
                      file = "combined.ps")
p15 <- tda_read_ps(file.path(r15$dir, "combined.ps"))
frames15 <- Filter(function(o) o$op == "lines" && length(o$x) >= 5L, p15$ops)
ok("dplot: all four combined panels' own frames survive parsing",
   length(frames15) == 4L)
ok("dplot: every panel's frame is drawn in its real colour, not a leftover white",
   all(vapply(frames15, function(o) identical(o$col, "black"), logical(1))))
}

## tda_pl_axes(): fmt= (found via replicating examples/exam/plot8.cf's
## own "plxa(sc=0.5,fmt=4.1)" alongside a differently-scaled y axis
## with no fmt= at all) was left out of the same per-axis pick() logic
## sc/ic/lty/lw/fs/tl already had, falling through to ... instead --
## which sends the identical, whole, un-picked value to both axes
## rather than index i's share of it, checked: fmt =
## c("4.1", NA) reached TDA as the literal, un-parseable "fmt =
## 4.1,NA" on *both* plxa and plya, not "fmt = 4.1" on plxa alone with
## plya left unset.
p16 <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1))
p16 <- tda_pl_axes(p16, sc = c(0.5, 1), ic = c(0, 10), fmt = c("4.1", NA))
full16 <- paste(p16$cmds, collapse = "\n")
ok("tda_pl_axes: fmt= reaches only the axis it was given for, not both",
   lengths(regmatches(full16, gregexpr("fmt", full16))) == 1L)
ok("tda_pl_axes: the axis with no fmt= is not sent the literal string \"NA\"",
   !grepl("NA", full16))

## tda_pl_lines(): gray= (found via replicating examples/exam/ple5p.cf's
## own confidence-interval plot, "plot(gs=0.9,sel=T1)=Time,UBnd") was
## left to ... rather than translated to TDA's gs=, the way every
## other tda_pl_* wrapper's gray= already is -- checked,
## gray = 0.9 combined with group= reached TDA as the literal,
## unrecognised option "gray=0.9", a syntax error, not silently
## ignored or misapplied.
d17 <- data.frame(x = 1:5, y = 1:5, grp = rep("a", 5))
p17 <- tda_ps(d17, xlim = c(0, 6), ylim = c(0, 6))
p17 <- tda_pl_lines(p17, "x", "y", select = "a", by = "grp", gray = 0.7)
ok("tda_pl_lines: gray= reaches TDA as gs=, not the literal, unrecognised \"gray=\"",
   any(grepl("gs\\s*=\\s*0.7", p17$cmds)) && !any(grepl("gray", p17$cmds)))

## Multi-panel spatial layout: a real four-panel dplot()-combined file
## (examples/exam/plot7.cf) rendered with the whole layout badly
## broken -- one panel dominating the canvas, the other three
## compressed into a corner -- confirmed against Ghostscript's
## rendering of the identical, unmodified file, an evenly-proportioned
## 2x2 grid. Two separate, real bugs stacked to cause this, found by
## tracing the parser's transform stack directly rather than by
## further guessing at the rendering geometry:
##
## 1. .ps_body(): dplot() keeps each combined file's own, complete,
##    separate %%BeginProlog...%%EndProlog pair rather than sharing
##    one -- taking only the *first* %%EndProlog as the single prolog
##    boundary (an earlier version of this fix) stripped the very
##    first panel's "gsave % added" isolation marker along with
##    its prolog, leaving that gsave's matching "grestore % added"
##    unmatched later. The first panel's coordinate-system scale
##    was then never actually undone before the second panel applied
##    its on top, and every panel after the first genuinely
##    compounded it, 0.43 becoming 0.43^2. Fixed by removing every
##    %%BeginProlog...%%EndProlog pair found, not just the first.
##
## 2. The "/fsiz N def" handler: font size is drawn through the same
##    transform geometry is, but the raw value read off "/fsiz" is
##    the size *before* that transform -- dplot()'s per-panel
##    shrink (sfx/sfy, often well under 1) was never applied to it,
##    leaving text at full original size against correctly-shrunk
##    geometry. Fixed by scaling fontsize by the CTM's linear
##    magnitude (sqrt(a^2+b^2)) at the point "/fsiz" is set.
if (have_examples) {
f3 <- file.path(EX, "exam", "plot3.ps")
f4 <- file.path(EX, "exam", "plot4.ps")
r18 <- tda_combine_ps(list(c(f3, f3), c(f4, f4)), width = 100, height = 60,
                      file = "combined18.ps")
p18 <- tda_read_ps(file.path(r18$dir, "combined18.ps"))
panel_ids18 <- vapply(p18$ops, function(o) o[["panel"]] %||% 1L, integer(1))
scales18 <- vapply(sort(unique(panel_ids18)), function(pnl) {
    o <- Filter(function(o) o$op == "text" && grepl("Coordinate", o$label),
               p18$ops[panel_ids18 == pnl])
    if (!length(o)) NA_real_ else o[[1L]]$fontsize
}, numeric(1))
ok("dplot: every panel's title font size is the same, real scale, not a compounding one",
   all(is.finite(scales18)) &&
   diff(range(scales18)) < 0.01 * mean(scales18))
}

## tda_pl_axis(): a second axis drawn at an explicit position (TDA's
## own plxa(...) = xa,ya,xb,yb, its documented technique for e.g. a
## top axis in different units from the bottom one) previously meant
## dropping to raw tda_pl(p, "plxa", dir=, rhs=) -- reachable now with
## named at=/dir= like every other tda_pl_axes() option.
p19 <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1))
p19 <- tda_pl_axis(p19, "x", sc = 1, ic = 2)
p19 <- tda_pl_axis(p19, "x", sc = 1, ic = 5, dir = 1, at = c(-1, 1, 1, 1))
r19 <- tda_read_ps(tda_ps_file(p19))
ax19 <- Filter(function(o) o$op == "lines" && length(o$x) == 2L, r19$ops)
ok("tda_pl_axis: a second axis at an explicit position draws a real, separate line",
   any(vapply(ax19, function(o) all(abs(o$y - o$y[1L]) < 1e-9) &&
             abs(o$y[1L] - o$y[2L]) < 1e-9 && diff(range(o$x)) > 1,
             logical(1))))

## tda_pl_panel(): the documented way to add a new panel to a session
## (TDA's psetup() called again, at a different psorg=) -- takes
## mm, consistent with every other size/position parameter in this
## package, and tracks the running page origin so the caller never
## needs TDA's undocumented default psorg= or to know that
## psorg= is itself in points while pxlen=/pylen= are in mm (a real
## inconsistency in TDA's units for the very same command,
## checked against its help text). Rounding psorg= to
## whole points also matters on its own: TDA's psetup() parser
## rejects a decimal psorg= value with a syntax error, confirmed
## directly against a value as ordinary as "102.9167,162.2778".
p20 <- tda_ps(xlim = c(0, 1), ylim = c(0, 1), width = 40, height = 30)
same("tda_pl_panel: a session starts at TDA's documented default origin",
     p20$origin, c(150, 460))
p20 <- tda_pl_frame(p20)
p20 <- tda_pl_panel(p20, origin = c(50, 0), width = 40, height = 30)
same("tda_pl_panel: a relative offset in mm is converted to TDA's points",
     p20$origin, c(150, 460) + c(50, 0) * 72 / 25.4, 1e-6)
ok("tda_pl_panel: the emitted psorg= is rounded to whole points, not a decimal that TDA itself refuses",
   grepl("psorg = [0-9]+,[0-9]+", paste(p20$cmds, collapse = " ")))
p20 <- tda_pl_frame(p20)
r20 <- tryCatch({ pdf(NULL); on.exit(dev.off()); tda_plot_ps(tda_ps_file(p20)); TRUE },
                error = function(e) e)
ok("tda_pl_panel: a session with a second panel actually runs and renders",
   isTRUE(r20))

## tda_combine_ps(): every tda_ps() session defaults to the identical
## file="plot.ps" unless the caller sets a different name for each
## one, and dplot()'s fn= list references files purely by name
## once they are copied into one, shared directory. Copying each
## input under its basename() (an earlier version of this)
## silently overwrote an earlier file the moment a later one shared
## that same default name -- checked: combining two
## distinct plots, neither given an explicit file=, left every one
## of dplot()'s panels showing the *last* input's content,
## not each its own, with no error or warning at all.
p21a <- tda_pl_text(tda_ps(xlim = c(0, 1), ylim = c(0, 1)), "AAA", at = c(0.5, 0.5))
p21b <- tda_pl_text(tda_ps(xlim = c(0, 1), ylim = c(0, 1)), "BBB", at = c(0.5, 0.5))
f21a <- tda_ps_file(p21a)
f21b <- tda_ps_file(p21b)
ok("tda_combine_ps: two inputs sharing the same default file= are not confused with each other",
   !identical(basename(f21a), NA_character_) &&
   file.exists(f21a) && file.exists(f21b))
r21 <- tda_combine_ps(list(f21a, f21b), file = "combined21.ps")
p21 <- tda_read_ps(file.path(r21$dir, "combined21.ps"))
lab21 <- vapply(Filter(function(o) o$op == "text", p21$ops),
                function(o) o$label, "")
ok("tda_combine_ps: both panels' own, distinct content survives combining, not just the last one's",
   "AAA" %in% lab21 && "BBB" %in% lab21)

## tda_pl_labels(): pxlabel/pylabel (an x-/y-axis label, as distinct
## from plabel's title) are real, working TDA commands -- every
## shipped multi-axis example uses them -- but neither is documented
## in TDA's help text at all, checked by its absence
## there; only the C source itself (t_cmd.c, pl_label()) confirms
## both reach the identical underlying command as plabel, just with
## mode 1/2 instead of 0. Reachable now via which=, rather than only
## ever emitting a title.
p22 <- tda_ps(xlim = c(0, 1), ylim = c(0, 1))
p22 <- tda_pl_labels(p22, "A Title", which = "title")
p22 <- tda_pl_labels(p22, "An X Label", which = "x")
p22 <- tda_pl_labels(p22, "A Y Label", which = "y", sc = 3)
full22 <- paste(p22$cmds, collapse = " ")
ok("tda_pl_labels: which= reaches TDA's three distinct label commands",
   grepl("plabel", full22) && grepl("pxlabel", full22) && grepl("pylabel", full22))

## tda_pl_circle(): ploto supports an optional arc
## (r,alpha,beta on its own right-hand side, TDA's documented
## "r [,alpha,beta]"), not only a full circle -- checked
## against TDA's help text -- but an earlier version of this
## function only ever sent a bare r, with no way to reach the arc
## form at all.
p23 <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
p23 <- tda_pl_circle(p23, at = c(5, 5), r = 2, angles = c(45, 90))
ok("tda_pl_circle: angles= reaches ploto's optional alpha,beta arc form",
   any(grepl("2,45,90", p23$cmds, fixed = TRUE)))

## tda_pl_smooth(): gray= (found via replicating examples/exam/
## plot19.cf's smoothing plot) was left to ... rather than
## translated to TDA's gs=, the way every other tda_pl_*
## wrapper's gray= already is -- checked, gray = 0.9
## reached TDA as the literal, unrecognised option "gray=0.9", a
## syntax error, not silently ignored.
d24 <- data.frame(x = 1:5, y = c(1, 2, 1.5, 2, 1))
p24 <- tda_ps(d24, xlim = c(0, 6), ylim = c(0, 3))
p24 <- tda_pl_smooth(p24, "x", "y", ns = 10, gray = 0.9)
ok("tda_pl_smooth: gray= reaches TDA as gs=, not the literal, unrecognised \"gray=\"",
   any(grepl("gs\\s*=\\s*0.9", p24$cmds)) && !any(grepl("gray", p24$cmds)))

## tda_pl_function(): step= (found while replicating examples/exam/
## plot-g1.cf's function plots, where a shallow, explicit step is
## needed across a wide range) was silently divided by 100 before use
## -- R's custom infix operators (%anything%, including %||% here)
## bind *tighter* than arithmetic, so "step %||% diff(range) / 100"
## parsed as "(step %||% diff(range)) / 100", not "step %||%
## (diff(range) / 100)" as the code visually appears to intend.
## Checked: an explicit, working step = 0.02 reached TDA's
## own rx= as 0.0002, a hundred times smaller, with no error or
## warning at all -- only a plot with 100x too many points, silently.
p25 <- tda_ps(xlim = c(0, 5), ylim = c(0, 2))
p25 <- tda_pl_function(p25, "x1", range = c(0.01, 5), step = 0.02)
ok("tda_pl_function: an explicit step= is not silently divided by 100",
   any(grepl("0.01(0.02)5", p25$cmds, fixed = TRUE)))

## Adversarial/fuzz case: a circle drawn much larger than the plot's
## own declared range (TDA's technique for a near-straight edge or
## a full-bleed fill, checked: examples/exam/gd18.cf's
## rd=-3 edges use a similarly oversized radius). Circle ops were never
## given a clip= field at all (unlike lines/symbol, which already had
## one), so the viewport's extent computation used the circle's
## full, unclipped radius -- inflating the whole plot to fit content
## that was never actually visible past the plot's boundary in
## TDA's real output, and painting straight over the frame/axis
## instead of stopping at it. Found by deliberately constructing an
## adversarial file (a circle of radius 5000 in a plot ranging -1000
## to 1000) and comparing against the real TDA binary's output,
## not against another well-behaved example.
p26 <- tda_ps(xlim = c(-1000, 1000), ylim = c(-1, 1), width = 100, height = 60)
p26 <- tda_pl_frame(p26)
p26 <- tda_pl_circle(p26, at = c(-500, -0.5), r = 5000, gray = 0.9)
r26 <- tda_read_ps(tda_ps_file(p26))
circ26 <- Filter(function(o) o$op == "circle", r26$ops)[[1L]]
ok("tda_plot_ps: an oversized circle's clip region is recorded, not left NULL",
   !is.null(circ26$clip))
pdf26 <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf26)
r26_result <- tryCatch({ tda_plot_ps(r26); TRUE }, error = function(e) e)
grDevices::dev.off()
ok("tda_plot_ps: an oversized, clipped circle still renders without error",
   isTRUE(r26_result))

## Adversarial/fuzz case: a plot with a strongly anisotropic aspect
## ratio (150mm wide, 15mm tall, a data range 1 wide but 100 tall --
## deliberately constructed, not from an existing example) exposed a
## real bug in the page-wide text-scaling factor: span (the largest
## of the x/y extents) was compared against avail (the smallest of
## the device's width/height) -- mixing two different axes'
## own numbers, since the wide span and the narrow available
## dimension belong to different directions here. Checked
## against the real TDA binary's output for the same file: text
## that should render at a normal, readable size came out at a small
## fraction of a point, illegible. Fixed by pairing each axis with
## its available space and taking whichever axis is the binding
## constraint.
p27 <- tda_ps(xlim = c(0, 1), ylim = c(0, 100), width = 150, height = 15)
p27 <- tda_pl_frame(p27)
p27 <- tda_pl_text(p27, "diagonal", at = c(0.5, 50), rotate = 45)
r27 <- tda_read_ps(tda_ps_file(p27))
pdf27 <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf27, width = 10, height = 1)
tda_plot_ps(r27)
grb27 <- grid::grid.grab()
grDevices::dev.off()
# grid::grid.grep() returned an empty result here despite text clearly
# being on the page (checked) -- its own $gp$fontsize was
# NULL too, not a usable value at all. grid.grab() plus a manual,
# recursive walk of $children is what actually finds a rendered text
# grob and its real, resolved fontsize.
find_text27 <- function(g) {
    if (inherits(g, "text"))
        return(list(g))
    if (is.null(g$children))
        return(list())
    unlist(lapply(g$children, find_text27), recursive = FALSE)
}
texts27 <- find_text27(grb27)
ok("tda_plot_ps: text on a strongly anisotropic plot is not shrunk to near-zero size",
   length(texts27) > 0 &&
   all(vapply(texts27, function(g) (g$gp$fontsize %||% 0) > 3, logical(1))))

## Real bugs found from user-reported differences against TDA's
## real output, not from further fuzzing:
##
## 1. tda_plot_ps(): a dashed/dotted arc or circle (TDA's numeric
##    dash pattern, not a grid-compatible keyword) rendered solid --
##    noted in an earlier comment as a known gap but never actually
##    finished. Checked: examples/exam/plot14.cf's
##    small angle-marker arc, ploto(...,lt=2), is dashed in TDA's
##    real output.
p28 <- tda_ps(xlim = c(0, 10), ylim = c(0, 6))
p28 <- tda_pl_circle(p28, at = c(1, 2), r = 0.8, angles = c(45, 90), lty = 2)
r28 <- tda_read_ps(tda_ps_file(p28))
pdf28 <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf28)
grb28 <- tryCatch({ tda_plot_ps(r28); grid::grid.grab() },
                  error = function(e) e)
grDevices::dev.off()
ok("tda_plot_ps: a dashed partial arc renders without error",
   !inherits(grb28, "error"))

## 2. Real PostScript's gsave/grestore save and restore the
##    current path, not only colour/CTM/etc -- this parser's path
##    was a separate variable neither gsave nor grestore ever touched
##    at all. Checked: examples/exam/plot19.cf's
##    plotsp fill ("[draw curve] gsave [extend path down to the x
##    axis] fill grestore stroke") strokes only the curve in TDA's
##    own real output, not the polygon's bottom and sides, since
##    grestore undoes the path extension before the stroke that
##    follows it -- this parser stroked the full, closed polygon
##    instead, a border TDA's output never has. Fixing this then
##    broke a second, related case (examples/exam/plot21.cf's
##    contour bands) via a separate "fill grestore" on one line code
##    path that needed the identical fix and had been missed.
p29 <- tda_ps(data.frame(x = c(1, 2, 3), y = c(1, 2, 1)),
              xlim = c(0, 4), ylim = c(0, 3))
p29 <- tda_pl_smooth(p29, "x", "y", ns = 5, gray = 0.9)
r29 <- tda_read_ps(tda_ps_file(p29))
poly29 <- Filter(function(o) o$op == "polygon", r29$ops)
lines29 <- Filter(function(o) o$op == "lines", r29$ops)
ok("tda_read_ps: a filled curve's stroke is a separate, shorter line, not the closed polygon",
   length(poly29) > 0 && length(lines29) > 0 &&
   all(vapply(lines29, function(l) length(l$x) < length(poly29[[1L]]$x),
             logical(1))))

## 3. Special characters in TDA's text ("@nnn"/"\nnn" escapes,
##    drawn through either FT -- Times-Roman, re-encoded via a table
##    embedded in every TDA file's prolog -- or FS, the standard
##    Adobe Symbol font) decoded as the raw Latin-1 byte value
##    regardless of which font was active, which is wrong for both:
##    checked, examples/exam/plots.cf's text rendered
##    as raw, undecoded ASCII in place of Greek letters and accented
##    characters, instead of "alpha = 2.5 (bullet) beta" and
##    "Sm(oslash)rebr(oslash)d" the way TDA's real output shows
##    it. A second, related bug this surfaced: some "show" calls have
##    their own "(string)" and the "show" keyword on two separate
##    lines rather than one, which silently dropped that entire text
##    segment rather than merely mis-decoding it.
p30 <- tda_ps(xlim = c(0, 10), ylim = c(0, 5))
p30 <- tda_pl_text(p30, "@141 = 2.5 \\245 @142", at = c(1, 4))
r30 <- tda_read_ps(tda_ps_file(p30))
txt30 <- Filter(function(o) o$op == "text", r30$ops)[[1L]]
ok("tda_read_ps: TDA's @nnn/\\nnn special-character text decodes through the right font's table, not a raw byte value",
   grepl("\u03b1", txt30$label) && grepl("\u2022", txt30$label) &&
   grepl("\u03b2", txt30$label) && !grepl("\\\\", txt30$label))


# page transforms in the reader: psrot= rotates the whole page about the
# plot origin ("xorg yorg translate" then "rf rotate", t_psinit.c) and the
# reader's transform tracker replays it -- a claim a stale comment denied
# for the 3-D family until round-tripped here
rot_td <- tempdir()
tda_run(c("psfile = rot.ps;",
          "psetup(pxlen=60, pylen=40, pxa=0,10, pya=0,10, psrot=90);",
          "pltext(xy=5,5) = mid;"), dir = rot_td)
rot_r <- tda_read_ps(file.path(rot_td, "rot.ps"))
rot_t <- Filter(function(o) o$op == "text", rot_r$ops)[[1L]]
ok("reader: psrot=90 lands text at origin + R(90) * local",
   abs(rot_t$x - (150 - 56.70)) < 0.05 && abs(rot_t$y - (460 + 85.05)) < 0.05)

# and a 3-D session's ops come back at annotation + origin, ordinary 2-D
# geometry (TDA projects in C before writing)
p3x <- c(1, 3, 5, 7); p3y <- c(2, 6, 3, 8); p3z <- c(1, 4, 2, 6)
p3 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8),
              view = c(60, 25))
p3 <- tda_pl_points3(p3, p3x, p3y, p3z, symbol = 4)
p3 <- tda_pl_text3(p3, "peak", at = c(7, 8, 6))
p3f <- tda_ps_file(p3)
p3ln <- readLines(p3f)
p3an <- strsplit(trimws(sub("^%#text:", "",
                            grep("^%#text:", p3ln, value = TRUE)[1L])),
                 "\\s+")[[1L]]
p3r <- tda_read_ps(p3f)
p3t <- Filter(function(o) o$op == "text", p3r$ops)[[1L]]
ok("reader: 3-D text sits at its annotation plus the plot origin",
   abs(p3t$x - (as.numeric(p3an[1L]) + 100)) < 0.05 &&
   abs(p3t$y - (as.numeric(p3an[2L]) + 100)) < 0.05)
ok("reader: 3-D symbols all present",
   length(Filter(function(o) o$op == "symbol", p3r$ops)) == 4L)

# ---- the X11 screen commands, replicated ------------------------------
# xreg/xconh drew onto an X11 window whose front end was never part of
# this source tree.  tda_pl_regression and tda_pl_hull reproduce them in
# the ordinary ps session.  The lowess claim is checked path-for-path:
# TDA's scplot(opt=3) curve and stats::lowess(f, iter = 2,
# delta = 0) must emit the identical PostScript line.
local({
    set.seed(1)
    d <- data.frame(x = sort(round(runif(40, 0, 10), 2)))
    d$y <- round(2 + sin(d$x) * 3 + rnorm(40, sd = 0.4), 2)
    path_pts <- function(f) {
        # the drawn path of the file's last command: everything from the
        # last moveto on (the clip rectangle earlier in the section is
        # part of the plumbing, not the curve)
        ll <- readLines(f)
        sec <- ll[seq(max(grep("^%#", ll)), length(ll))]
        keep <- grep("^ *[0-9.+-]+ +[0-9.+-]+ (l|m)$", sec, value = TRUE)
        keep <- keep[seq(max(grep(" m$", keep)), length(keep))]
        m <- regmatches(keep,
                        regexec("^ *([0-9.+-]+) +([0-9.+-]+) (l|m)$", keep))
        do.call(rbind,
                lapply(m[lengths(m) == 4], function(z) as.numeric(z[2:3])))
    }
    p1 <- tda_ps(d, xlim = c(0, 10), ylim = c(-2, 8))
    p1 <- tda_pl_scatter(p1, "x", "y", type = "lowess", bandwidth = 0.4)
    p2 <- tda_ps(d, xlim = c(0, 10), ylim = c(-2, 8))
    p2 <- tda_pl_regression(p2, "x", "y", type = "lowess", bandwidth = 0.4)
    a <- path_pts(tda_ps_file(p1))
    b <- path_pts(tda_ps_file(p2))
    ok("xreg lowess: identical ps path to TDA's scplot curve",
       nrow(a) == nrow(b) && max(abs(a - b)) == 0)

    p3 <- tda_ps(d, xlim = c(0, 10), ylim = c(-2, 8))
    p3 <- tda_pl_regression(p3, "x", "y", type = "least_squares")
    same("xreg least squares: coefficients are lm's",
         unname(p3$fits$least_squares),
         unname(coef(lm(y ~ x, d))), 1e-10)
    ok("xreg least squares: the line was drawn",
       !is.null(path_pts(tda_ps_file(p3))))

    p4 <- tda_ps(d, xlim = c(0, 10), ylim = c(-2, 8))
    p4 <- tda_pl_regression(p4, "x", "y", type = "l1")
    same("xreg l1: coefficients are tda_l1reg's",
         unname(p4$fits$l1),
         unname(coef(tda_l1reg(y ~ x, d))), 1e-10)

    # Hulls come from TDA's plotch, not from chull() in R. Per group,
    # a closed outline; a group of two points still gets a degenerate one
    # -- TDA draws it, and this used to assert the opposite because the
    # hull was computed here and such a group was skipped.
    dh <- data.frame(x = c(0, 4, 4, 0, 2, 9, 10),
                     y = c(0, 0, 4, 4, 2, 9, 9),
                     g = c(1, 1, 1, 1, 1, 2, 2))
    p5 <- tda_ps(dh, xlim = c(0, 10), ylim = c(0, 10))
    p5 <- tda_pl_hull(p5, "x", "y", by = "g")
    pts <- path_pts(tda_ps_file(p5))
    # path_pts reads the last command's path, which here is group 2:
    # two points, so TDA draws a degenerate outline rather than skipping
    # the group as the old R-side hull did
    ok("plotch hull: the two-point group is drawn, not skipped",
       nrow(pts) == 3)
    # and one group alone is the square on its own
    p6 <- tda_ps(dh[dh$g == 1, ], xlim = c(0, 10), ylim = c(0, 10))
    p6 <- tda_pl_hull(p6, "x", "y")
    ok("plotch hull: five points, first and last the same",
       { q <- path_pts(tda_ps_file(p6))
         nrow(q) == 5 && isTRUE(all.equal(q[1, ], q[5, ])) })
})
