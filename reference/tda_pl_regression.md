# Add a fitted regression curve or convex hulls to a plot

`tda_pl_regression` draws the fitted curve that TDA's screen-only `xreg`
command used to add to an X11 scatterplot – the X11 front end was never
part of this source tree, so the three fits it offered are reproduced
here and drawn through the ordinary PostScript session: a lowess curve,
a least-squares line, or an L1-norm (least-absolute-deviations) line,
one per group. The lowess curve is TDA's, drawn by `scplot` – the same
one
[`tda_pl_scatter`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md)
draws, and the one `xreg` drew. (`npreg(opt = 4)` runs the same C
routine but cannot give the same curve: it resets its band width to 1
when it is below epsilon, before the lowess branch, so the delta
shortcut is never 0, while `scplot` passes 0.) The least-squares line is
ordinary regression of y on x; the L1 line is TDA's `l1reg` via
[`tda_l1reg`](https://janmarvin.github.io/TDA/reference/tda_l1reg.md).
Like `xreg`, the straight-line fits are drawn from the group's smallest
x to its largest.

## Usage

``` r
tda_pl_regression(
  p,
  x,
  y,
  type = c("lowess", "least_squares", "l1"),
  by = NULL,
  bandwidth = 0.5,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_hull(
  p,
  x,
  y,
  by = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  symbol = NULL,
  size = NULL,
  smooth = NULL,
  expand = NULL,
  ...
)
```

## Arguments

- p:

  A plot session started with
  [`tda_ps`](https://janmarvin.github.io/TDA/reference/tda_ps.md).

- x, y:

  Column names in the session's data.

- type:

  The fit: `"lowess"` (the default, as in `xreg`), `"least_squares"`, or
  `"l1"`.

- by:

  Optional name of a grouping column; one fit (or hull) per group.
  Without it, all points are one group.

- bandwidth:

  The lowess span, `xreg`'s `sig=` (default 0.5, clamped there to (0,
  1)).

- lty, lw, gray:

  Line type, width (mm), and gray level for the drawn curve.

- ...:

  Further options passed to the drawing command.

- symbol, size:

  for `tda_pl_hull`, a marker symbol and its size in millimetres drawn
  at the points; `plotch`'s `s=` and `fs=`.

- smooth:

  the number of intervals for Akima smoothing of the hull's outline;
  `plotch`'s `ns=`.

- expand:

  a margin added to a smoothed hull, in millimetres; TDA's `ic=`. Each
  point of the smoothed outline is pushed away from the hull's centroid
  in proportion to its distance from it along each axis, so the outline
  grows without changing shape (t_plot.c, pl_plotch).

## Value

The session, with the curve(s) added. The fitted coefficients (for the
straight-line types) or curves are kept in `p$fits`, one entry per call
and group.

## Details

`tda_pl_hull` is `xconh`: the convex hull of each group's points, drawn
as a closed outline.

The third screen command, `xplotf`, plotted values straight from a file;
[`tda_pl_scatter`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md)
and
[`tda_pl_lines`](https://janmarvin.github.io/TDA/reference/tda_pl.md)
with raw vectors already do that.

## See also

Other plotting:
[`plot.tda_ple()`](https://janmarvin.github.io/TDA/reference/plot.tda_ple.md),
[`tda_check_ps()`](https://janmarvin.github.io/TDA/reference/tda_check_ps.md),
[`tda_pl()`](https://janmarvin.github.io/TDA/reference/tda_pl.md),
[`tda_pl_arc()`](https://janmarvin.github.io/TDA/reference/tda_pl_arc.md),
[`tda_pl_axis()`](https://janmarvin.github.io/TDA/reference/tda_pl_axis.md),
[`tda_pl_graph()`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md),
[`tda_pl_hist`](https://janmarvin.github.io/TDA/reference/tda_pl_hist.md),
[`tda_pl_panel()`](https://janmarvin.github.io/TDA/reference/tda_pl_panel.md),
[`tda_pl_scatter()`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md),
[`tda_plot_ps()`](https://janmarvin.github.io/TDA/reference/tda_plot_ps.md),
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md),
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = round(runif(40, 0, 10), 2))
d$y <- round(1 + 0.8 * d$x + rnorm(40), 2)
p <- tda_ps(d, xlim = c(0, 10), ylim = c(0, 12))
p <- tda_pl_frame(p)
p <- tda_pl_scatter(p, "x", "y", symbol = 1)
p <- tda_pl_regression(p, "x", "y", type = "least_squares",
                       lty = "dashed")
p$fits
#> $least_squares
#>     alpha      beta 
#> 1.3000324 0.7544393 
#> 
```
