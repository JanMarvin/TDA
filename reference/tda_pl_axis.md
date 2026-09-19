# Draw a single axis, optionally at an explicit position

`tda_pl_axes` draws both axes together, at the plot's boundary.
`tda_pl_axis` draws one axis at a time and, via `at`, can place it
anywhere in the plot rather than only at the boundary – TDA's technique
for a second x or y axis (e.g. a top axis in different units from the
bottom one), which `tda_pl_axes` alone cannot reach.

## Usage

``` r
tda_pl_axis(
  p,
  which = c("x", "y"),
  sc = NULL,
  ic = NULL,
  lty = NULL,
  lw = NULL,
  fs = NULL,
  tl = NULL,
  fmt = NULL,
  at = NULL,
  dir = NULL,
  ...
)
```

## Arguments

- p:

  a [`tda_ps`](tda_ps.md) session.

- which:

  `"x"` or `"y"`, which axis to draw.

- sc, ic, lty, lw, fs, tl, fmt:

  as in [`tda_pl_axes`](tda_pl.md), for this one axis.

- at:

  the axis line's endpoints, `c(xa, ya, xb, yb)` in data units, when it
  is not simply the plot's boundary (TDA's `plxa(...) = xa,ya,xb,yb`
  right-hand side). Omit for an ordinary axis at the plot boundary.

- dir:

  which side of the axis line the tick labels are written on (TDA's
  `dir=`); only meaningful together with `at`.

- ...:

  further options for `plxa`/`plya`.

## Value

The session, with the command added.

## See also

Other plotting: [`plot.tda_ple()`](plot.tda_ple.md),
[`tda_check_ps()`](tda_check_ps.md), [`tda_pl()`](tda_pl.md),
[`tda_pl_arc()`](tda_pl_arc.md), [`tda_pl_graph()`](tda_pl_graph.md),
[`tda_pl_hist`](tda_pl_hist.md), [`tda_pl_panel()`](tda_pl_panel.md),
[`tda_pl_regression()`](tda_pl_regression.md),
[`tda_pl_scatter()`](tda_pl_scatter.md),
[`tda_plot_ps()`](tda_plot_ps.md), [`tda_ps()`](tda_ps.md),
[`tda_ps3()`](tda_ps3.md), [`tda_read_ps()`](tda_read_ps.md)

## Examples

``` r
p <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1))
p <- tda_pl_axis(p, "x", sc = 1, ic = 2)
p <- tda_pl_axis(p, "y", sc = 0.5, ic = 0)
# a second x axis, drawn at y = 1 instead of the plot's bottom edge
p <- tda_pl_axis(p, "x", sc = 1, ic = 5, dir = 1, at = c(-1, 1, 1, 1))
```
