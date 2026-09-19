# Curved arcs

`plotk`: draws an arc, optionally ending in an arrow, between each
consecutive pair of points – a curved alternative to
[`tda_pl_lines`](tda_pl.md) for showing a directed relationship without
the line hiding what is under it. Fixed along the way: a genuine
rendering bug where a shallow arc (TDA's technique for one – a small
slice of a large-radius circle) was measured by the full circle's radius
when sizing the plot, inflating the bounding box to fit an almost
entirely off-screen circle and shrinking everything else into a corner;
the fix (using the arc segment's actual extent) is general, not specific
to this function.

## Usage

``` r
tda_pl_arc(p, x, y, curvature = 0, arrow = NULL, lty = NULL, lw = NULL, ...)
```

## Arguments

- p:

  a `tda_ps` session.

- x, y:

  coordinates the arcs connect, in order; at least two points.

- curvature:

  how much the arc bows out, as a fraction of the distance between the
  two points; `0` is a straight line. TDA's own `sc=`.

- arrow:

  size of an arrowhead at the end of each arc, as `c(length, width)` in
  mm; no arrowhead if not given.

- lty, lw:

  line type and line width in mm for the arc.

- ...:

  further options for the command.

## Value

The session, with the command added.

## See also

Other plotting: [`plot.tda_ple()`](plot.tda_ple.md),
[`tda_check_ps()`](tda_check_ps.md), [`tda_pl()`](tda_pl.md),
[`tda_pl_axis()`](tda_pl_axis.md), [`tda_pl_graph()`](tda_pl_graph.md),
[`tda_pl_hist`](tda_pl_hist.md), [`tda_pl_panel()`](tda_pl_panel.md),
[`tda_pl_regression()`](tda_pl_regression.md),
[`tda_pl_scatter()`](tda_pl_scatter.md),
[`tda_plot_ps()`](tda_plot_ps.md), [`tda_ps()`](tda_ps.md),
[`tda_ps3()`](tda_ps3.md), [`tda_read_ps()`](tda_read_ps.md)

## Examples

``` r
p <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
p <- tda_pl_frame(p)
p <- tda_pl_arc(p, c(2, 8), c(2, 8), curvature = 3, arrow = c(1.5, 1))
plot(p)
```
