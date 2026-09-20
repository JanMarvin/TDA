# Add a command to a plot

`tda_pl` reaches any of TDA's thirty-seven plot commands. The named
functions are the common ones and only spell out their arguments;
anything without a named wrapper goes through `tda_pl` with the command
name.

## Usage

``` r
tda_pl(p, cmd, ..., rhs = NULL)

tda_pl_axes(
  p,
  sc = NULL,
  ic = NULL,
  lty = NULL,
  lw = NULL,
  fs = NULL,
  tl = NULL,
  fmt = NULL,
  ...
)

tda_pl_frame(p, lty = NULL, lw = NULL, gray = NULL, ...)

tda_pl_grid(p, at_x = NULL, at_y = NULL, lty = NULL, lw = NULL, ...)

tda_pl_lines(
  p,
  x,
  y,
  select = NULL,
  by = NULL,
  lty = NULL,
  band = NULL,
  rows = NULL,
  gray = NULL,
  ...
)

TDA_LTY

tda_pl_points(
  p,
  x,
  y,
  symbol = 1,
  size = NULL,
  lty = 0,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_smooth(p, x, y, ns = 2, lty = NULL, gray = NULL, ...)

tda_pl_circle(
  p,
  at,
  r = 1,
  angles = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_text(
  p,
  label,
  at,
  fs = NULL,
  rotate = NULL,
  white = NULL,
  symbol = NULL,
  center = NULL,
  ...
)

tda_pl_rect(p, x, y, lty = NULL, lw = NULL, gray = NULL, rotate = NULL, ...)

tda_pl_polyline(p, x, y, lty = NULL, lw = NULL, gray = NULL, ...)

tda_pl_function(
  p,
  expr,
  range,
  step = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  deriv = NULL,
  ...
)

tda_ps_file(x)
```

## Arguments

- p:

  a `tda_ps` session.

- cmd:

  the TDA command name, e.g. `"plotd"`, `"plsurf3"`.

- ...:

  further options for the underlying TDA command, e.g. `lt=` (line type)
  or `lw=` (line width in mm) for `tda_pl_grid`. Each command's full
  option list is TDA's – look it up with `tda_help("plxa")`,
  `tda_help("plxgrid")`, and so on (the name after `"pl"` in each
  function above); this page only names the options with their own R
  argument.

- rhs:

  the command's right-hand side, when it takes one.

- sc:

  for `tda_pl_axes`, the distance between major tick marks (TDA's
  default: chosen automatically). For `tda_pl_labels`, the label's
  distance from the top of the plot, in mm (default 0). For
  `tda_pl_axes` specifically, this (and `ic`/`lty`/ `lw`/`fs`/`tl`
  below) can also be a pair, `c(x, y)`, to give the x and y axis their
  value instead of sharing one – `sc = c(1, 0.1)` for very different
  scales on the two axes, say. A single value (the usual case) is still
  used for both, unchanged.

- ic:

  for `tda_pl_axes`, the number of minor sub-intervals between major
  tick marks.

- lty:

  line type, either an R name – `"solid"`, `"dashed"`, `"dotted"`,
  `"dotdash"`, `"longdash"`, `"twodash"` – or a TDA number. With `by=`
  and no `lty`, the groups cycle through that list in order, as they
  would in R.

- lw:

  line width in mm, for `tda_pl_grid`, `tda_pl_axes`, `tda_pl_frame`,
  and `tda_pl_circle`.

- fs:

  font size in mm, for `tda_pl_axes` (axis labels) and `tda_pl_labels`
  (the label itself).

- tl:

  for `tda_pl_axes`, the length of the tick marks, in mm.

- fmt:

  for `tda_pl_axes`, the numeric format for tick labels (TDA's `fmt=` on
  `plxa`/`plya`, e.g. `"4.1"` for one decimal place in a field 4
  characters wide).

- gray:

  for `tda_pl_smooth`, grey scale value (TDA's `gs=`).

- at_x, at_y:

  for `tda_pl_grid`, where to draw the vertical and horizontal grid
  lines; default to the axis's tick positions.

- x, y:

  for `tda_pl_lines`/`tda_pl_points`, a column name in the session's
  data, or the numbers themselves – R's `lines(x = 1:10, y = 1:10)`. A
  vector the same length as the session's data is added as a new column;
  a different length gets its own, separate data matrix for just that
  one series (TDA's `clear;`+`nvar(dfile=...)` trick, restored
  afterward), useful for overlaying a fitted curve at a different
  resolution than the raw data it was fit to – see
  [`tda_npreg`](https://janmarvin.github.io/TDA/reference/tda_npreg.md)'s
  example. For everything else reachable through `tda_pl`, variables or
  coordinates, as names in the session's data or as numbers; for
  `tda_pl_rect`, two opposite corners.

- select:

  draw only this level of `by` – TDA's `sel=`, a row selection, said in
  R (the same name every other wrapper with a `sel=` uses).

- by:

  a grouping column in the plot's data. Without `select`, one line is
  drawn per level.

- band:

  a pair of column names holding the upper and lower bound of a
  confidence band, shaded behind the line the way
  `examples/ehhnew/ehc8.cf` does it.

- rows:

  an explicit row selection, for drawing two series of different lengths
  out of one frame.

- symbol:

  for `tda_pl_text`, draw a marker symbol at this position instead of
  (or as well as) the text (TDA's `s=` on `pltext`).

- size:

  symbol size.

- ns:

  degree of smoothing.

- at:

  where to draw it, as a pair of coordinates.

- r:

  radius of a circle, in the plot's units.

- angles:

  for `tda_pl_circle`, the start and end angle in degrees,
  counterclockwise, of a partial arc rather than a full circle (TDA's
  optional `alpha,beta` on `ploto`'s right-hand side, `r[,alpha,beta]`);
  omit for a full circle.

- label:

  text to draw.

- rotate:

  rotation in degrees, for `tda_pl_text` and `tda_pl_rect` (around their
  reference point) and `tda_pl_ellipse` (of its main axis).

- white:

  for `tda_pl_text`, draw the label on a white background rather than
  transparently.

- center:

  for `tda_pl_text`, centre the string on `at` rather than TDA's default
  of starting it there (TDA's `sc=1` on `pltext`).

- expr:

  an expression in TDA's language.

- range, step:

  for `tda_pl_function`, the range of the function's argument and the
  increment to evaluate it at – TDA's `rx = a(d)b`, built from these so
  the range does not need to be written in TDA's step notation by hand.

- deriv:

  for `tda_pl_function`, plot the function's first (`1`) or second (`2`)
  derivative instead of the function itself (TDA's `plotf1`/`plotf2`);
  omit to plot the function as written.

## Value

The session, with the command added.

## See also

Other plotting:
[`plot.tda_ple()`](https://janmarvin.github.io/TDA/reference/plot.tda_ple.md),
[`tda_check_ps()`](https://janmarvin.github.io/TDA/reference/tda_check_ps.md),
[`tda_pl_arc()`](https://janmarvin.github.io/TDA/reference/tda_pl_arc.md),
[`tda_pl_axis()`](https://janmarvin.github.io/TDA/reference/tda_pl_axis.md),
[`tda_pl_graph()`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md),
[`tda_pl_hist`](https://janmarvin.github.io/TDA/reference/tda_pl_hist.md),
[`tda_pl_panel()`](https://janmarvin.github.io/TDA/reference/tda_pl_panel.md),
[`tda_pl_regression()`](https://janmarvin.github.io/TDA/reference/tda_pl_regression.md),
[`tda_pl_scatter()`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md),
[`tda_plot_ps()`](https://janmarvin.github.io/TDA/reference/tda_plot_ps.md),
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md),
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
set.seed(1)
d <- data.frame(v = rnorm(60, 50, 10))
p <- tda_ps(d, xlim = c(0, 100), ylim = c(0, 0.05))
p <- tda_pl_axes(p)

# plotd (a density plot) has no named wrapper -- it needs an evaluation
# grid (x=) and exactly one variable on the right-hand side, which is not
# the R column name but the internal one tda_ps() assigned; look it up in
# the session's xlab/xname the same way the named wrappers do
vname <- p$xname[p$xlab == "v"]
p <- tda_pl(p, "plotd", x = "0(5)100", rhs = vname)
p
#> TDA plot session (plot.ps), 100 x 70 mm, x: [0,100], y: [0,0.05]
#> Data: 60 rows, v 
#> 3 drawing commands:
#>   plxa( );
#>   plya( );
#>   plotd( x = 0(5)100, ) = Vv;
#> (not yet rendered -- plot() or tda_ps_file() runs it)

# every named wrapper, on one plot, so each one is a working example
# rather than a bare signature -- x/y take a column name or a number
# everywhere, and every one of these actually draws something
d2 <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
p2 <- tda_ps(d2, xlim = c(0, 6), ylim = c(0, 12), width = 90, height = 70)
p2 <- tda_pl_axes(p2)              # tick marks and their numbers
p2 <- tda_pl_frame(p2)             # a box around the plotting area
p2 <- tda_pl_grid(p2, lw = 0.15)   # TDA's default (lw = 0.05mm) is
                                   # thin to the point of illegible on a
                                   # modern screen; heavier here on purpose
p2 <- tda_pl_lines(p2, "x", "y")   # the data, joined point to point
p2 <- tda_pl_points(p2, "x", "y", symbol = 1, size = 4)  # a marker per point
p2 <- tda_pl_smooth(p2, "x", "y", ns = 4, lty = "dashed")  # a smoothed fit
p2 <- tda_pl_circle(p2, at = c(3, 5), r = 0.3)    # mark a location
p2 <- tda_pl_text(p2, "peak", at = c(4.3, 9))     # a text label
p2 <- tda_pl_rect(p2, x = c(4, 5), y = c(11, 11.8))         # a rectangle
p2 <- tda_pl_polyline(p2, c(0.2, 1, 0.2), c(11, 11, 11.8))  # a raw shape
# its row along the top, so it doesn't cross the data line below
p2 <- tda_pl_function(p2, "sin(x)+11", range = c(0, 6), step = 0.2)
tda_ps_file(p2)   # the PostScript path, once the session has been run
#> [1] "/tmp/RtmpbyWI98/tda1f9b3f0dac8/plot.ps"
p2
#> TDA plot session (plot.ps), 90 x 70 mm, x: [0,6], y: [0,12]
#> Data: 5 rows, x, y 
#> 13 drawing commands:
#>   plxa( );
#>   plya( );
#>   plframe( );
#>   plxgrid( lw = 0.15, ) = 0,2,4,6,8,10,12;
#>   plygrid( lw = 0.15, ) = 0,1,2,3,4,5,6;
#>   plot( ) = Vx,Vy;
#>   plot( s = 1, lt = 0, fs = 4, ) = Vx,Vy;
#>   plots( ns = 4, lt = 5, ) = Vx,Vy;
#>   ploto( xy = 3,5, ) = 0.3;
#>   pltext( xy = 4.3,9, ) = "peak";
#>   plrec( ) = 4,11,1,0.800000000000001;
#>   plotp( ) = 0.2,11,1,11,0.2,11.8;
#>   plotf( rx = 0(0.2)6, ) = sin(x)+11;
#> (not yet rendered -- plot() or tda_ps_file() runs it)
```
