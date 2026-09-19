# A TDA plot

Opens a plot session. Drawing commands are added to it and it is run
when the plot is drawn or its file is asked for, which mirrors how TDA
works: a PostScript file is opened, a coordinate system set up, commands
issued, and the file closed.

## Usage

``` r
tda_ps(
  data = NULL,
  width = 100,
  height = 70,
  xlim = NULL,
  ylim = NULL,
  file = "plot.ps",
  origin = NULL,
  log = NULL,
  ...
)
```

## Arguments

- data:

  optional data frame the drawing commands may refer to.

- width, height:

  size of the plotting area in millimetres.

- xlim, ylim:

  ranges of the coordinate system.

- file:

  name for the PostScript file.

- origin:

  the session's physical position on the page, `c(x, y)` in mm from the
  page origin (TDA's `psorg=`). Only useful together with
  [`tda_pl_panel`](tda_pl_panel.md), since it fixes a known starting
  point that panel's `origin=` (relative to the previous panel) is added
  onto; omit it for an ordinary, single-panel session, where TDA's
  default page position is used.

- log:

  which axes are logarithmic (TDA's `pxa(log=1)`/ `pya(log=1)`, a
  per-axis option nested inside `pxa=`/ `pya=` rather than an ordinary
  `psetup` option): any of `"x"`, `"y"`, `"xy"`, or `""`/`NULL` (the
  default) for neither.

- ...:

  further options for `psetup`.

## Value

A `tda_ps` session.

## See also

Other plotting: [`plot.tda_ple()`](plot.tda_ple.md),
[`tda_check_ps()`](tda_check_ps.md), [`tda_pl()`](tda_pl.md),
[`tda_pl_arc()`](tda_pl_arc.md), [`tda_pl_axis()`](tda_pl_axis.md),
[`tda_pl_graph()`](tda_pl_graph.md), [`tda_pl_hist`](tda_pl_hist.md),
[`tda_pl_panel()`](tda_pl_panel.md),
[`tda_pl_regression()`](tda_pl_regression.md),
[`tda_pl_scatter()`](tda_pl_scatter.md),
[`tda_plot_ps()`](tda_plot_ps.md), [`tda_ps3()`](tda_ps3.md),
[`tda_read_ps()`](tda_read_ps.md)

## Examples

``` r
d <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
p <- tda_ps(d, xlim = c(0, 6), ylim = c(0, 10))
p <- tda_pl_axes(p)
p <- tda_pl_lines(p, "x", "y")
p
#> TDA plot session (plot.ps), 100 x 70 mm, x: [0,6], y: [0,10]
#> Data: 5 rows, x, y 
#> 3 drawing commands:
#>   plxa( );
#>   plya( );
#>   plot( ) = Vx,Vy;
#> (not yet rendered -- plot() or tda_ps_file() runs it)
```
