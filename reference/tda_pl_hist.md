# Further plot commands

Named wrappers for the plot commands that take arguments worth spelling
out. Everything else is reachable with
[`tda_pl`](https://janmarvin.github.io/TDA/reference/tda_pl.md).

## Usage

``` r
tda_pl_density(
  p,
  x,
  at,
  bandwidth = NULL,
  kernel = c("uniform", "triangle", "quartic", "epanechnikov"),
  scale = NULL,
  clip = TRUE,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_histogram(
  p,
  x,
  breaks,
  weights = NULL,
  scale = NULL,
  closed = c("left", "right"),
  vertical_lines = TRUE,
  clip = TRUE,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_contour(
  p,
  expr,
  levels,
  resolution = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_curve(
  p,
  x,
  y,
  symbol = NULL,
  size = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  smooth = NULL,
  ...
)

tda_pl_ellipse(
  p,
  at,
  axes = c(1, 0.5),
  rotate = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_labels(
  p,
  label,
  which = c("title", "x", "y"),
  sc = NULL,
  fs = NULL,
  ...
)

tda_pl_grid_lines(p, ...)
```

## Arguments

- p:

  a `tda_ps` session.

- x, y:

  variables in the session's data, or coordinates.

- at:

  evaluation points for a density, or a pair of coordinates.

- bandwidth:

  kernel bandwidth for `tda_pl_density`; TDA's default is 1.

- kernel:

  kernel shape for `tda_pl_density`: `"uniform"` (default),
  `"triangle"`, `"quartic"`, or `"epanechnikov"`.

- scale:

  for `tda_pl_histogram` or `tda_pl_density`, a scaling factor applied
  to the bar heights or curve (`ploth`/`plotd`'s `dscal=`, default 1).

- clip:

  for `tda_pl_histogram` or `tda_pl_density`, clip to the plot's
  declared range (the default); `FALSE` turns this off
  (`ploth`/`plotd`'s `nc=1`).

- lty, lw, gray:

  line type, line width in mm, and grey level from 0 (black) to 1
  (white), for `tda_pl_histogram`, `tda_pl_density`, `tda_pl_contour`,
  `tda_pl_curve`, and `tda_pl_ellipse`.

- ...:

  further options for the underlying TDA command – look them up with
  [`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md),
  e.g. `tda_help("ploth")` for `tda_pl_histogram`; see
  [`tda_pl`](https://janmarvin.github.io/TDA/reference/tda_pl.md).

- breaks:

  interval boundaries for a histogram.

- weights:

  for `tda_pl_histogram`, a variable in the session's data giving each
  case's weight.

- closed:

  for `tda_pl_histogram`, which side of each interval is closed:
  `"left"` (default) or `"right"` – `ploth`'s own `s=` (0/1).

- vertical_lines:

  for `tda_pl_histogram`, draw the vertical lines between bars (the
  default); `FALSE` omits them (`ploth`'s `ns=1`).

- expr:

  an expression in TDA's language.

- levels:

  contour levels for `tda_pl_contour`.

- resolution:

  for `tda_pl_contour`, the evaluation grid size, as `c(nx, ny)` or one
  number for both; TDA's default is `c(10, 10)`.

- symbol:

  marker symbol, for `tda_pl_curve`.

- size:

  marker size in mm, for `tda_pl_curve`.

- smooth:

  the number of intervals for Akima smoothing of a convex hull's
  outline; TDA's `ns=`. `plotch` smooths only when this is 2 or more AND
  the hull has more than two vertices – below either it draws the plain
  outline (t_plot.c, pl_plotch).

- axes:

  for `tda_pl_ellipse`, the half-lengths of the main and second axis,
  `c(a, b)`, in the plot's x and y units.

- rotate:

  rotation in degrees, for `tda_pl_ellipse`, of its main axis.

- label:

  the text `tda_pl_labels` places at the top of the plot.

- which:

  which label to set: `"title"` (TDA's `plabel`, at the top of the
  plot), `"x"` (`pxlabel`, below the x axis) or `"y"` (`pylabel`, left
  of the y axis). All three reach the same underlying TDA command
  (`pl_label()`; none of the three appear in TDA's help text, despite
  being real, working commands every shipped multi-axis example relies
  on), so they share every other option here too.

- sc:

  for `tda_pl_labels`, the label's distance from the top of the plot, in
  mm (default 0).

- fs:

  font size in mm, for `tda_pl_labels`.

## Value

The session, with the command added.

## See also

Other plotting:
[`plot.tda_ple()`](https://janmarvin.github.io/TDA/reference/plot.tda_ple.md),
[`tda_check_ps()`](https://janmarvin.github.io/TDA/reference/tda_check_ps.md),
[`tda_pl()`](https://janmarvin.github.io/TDA/reference/tda_pl.md),
[`tda_pl_arc()`](https://janmarvin.github.io/TDA/reference/tda_pl_arc.md),
[`tda_pl_axis()`](https://janmarvin.github.io/TDA/reference/tda_pl_axis.md),
[`tda_pl_graph()`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md),
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
d <- data.frame(v = round(rnorm(60) * 10 + 50, 1))
# ploth's bars are density-scaled (each bar's area, not its height, is
# its share of the cases), not raw counts -- ylim has to fit that
# small scale or the bars are real but too short to see
p <- tda_ps(d, xlim = c(0, 100), ylim = c(0, 0.05))
p <- tda_pl(p, "plxa", sc = 20, ic = 1)

# ploth takes the interval boundaries with breaks=, not tp=: passing tp=
# is silently ignored and every bar comes out flat.
tda_pl_histogram(p, "v", breaks = seq(0, 100, 10))
#> TDA plot session (plot.ps), 100 x 70 mm, x: [0,100], y: [0,0.05]
#> Data: 60 rows, v 
#> 2 drawing commands:
#>   plxa( sc = 20, ic = 1, );
#>   ploth( x = 0 (10) 100, ) = Vv;
#> (not yet rendered -- plot() or tda_ps_file() runs it)

# plotd is a kernel density estimate, and at= are the points it is
# evaluated at rather than the intervals of a histogram.
tda_pl_density(p, "v", at = seq(0, 100, 5))
#> TDA plot session (plot.ps), 100 x 70 mm, x: [0,100], y: [0,0.05]
#> Data: 60 rows, v 
#> 2 drawing commands:
#>   plxa( sc = 20, ic = 1, );
#>   plotd( x = 0 (5) 100, ) = Vv;
#> (not yet rendered -- plot() or tda_ps_file() runs it)

# a curve through data, a contour of an expression in two variables, an
# ellipse, a grid of light lines at the axis ticks, and a title
d2 <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
p2 <- tda_ps(d2, xlim = c(0, 6), ylim = c(0, 10))
p2 <- tda_pl_axes(p2)
p2 <- tda_pl_curve(p2, "x", "y")
p2 <- tda_pl_contour(p2, "x^2+y^2", levels = c(10, 30, 50))
p2 <- tda_pl_ellipse(p2, at = c(3, 5), axes = c(1.5, 0.8))
p2 <- tda_pl_grid(p2, lw = 0.15)
p2 <- tda_pl_labels(p2, "my plot title")
p2
#> TDA plot session (plot.ps), 100 x 70 mm, x: [0,6], y: [0,10]
#> Data: 5 rows, x, y 
#> 8 drawing commands:
#>   plxa( );
#>   plya( );
#>   plotch( ) = Vx,Vy;
#>   plotc( x = 10,30,50, ) = x^2+y^2;
#>   plote( xy = 3,5, ) = 1.5,0.8;
#>   plxgrid( lw = 0.15, ) = 0,2,4,6,8,10;
#>   plygrid( lw = 0.15, ) = 0,1,2,3,4,5,6;
#>   plabel( ) = "my plot title";
#> (not yet rendered -- plot() or tda_ps_file() runs it)
```
