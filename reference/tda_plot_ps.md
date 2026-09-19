# Draw TDA's PostScript output with grid

Replays what
[`tda_read_ps`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)
parsed onto the current grid device, so a plot TDA produced appears in
an R graphics device rather than only in a `.ps` file.

## Usage

``` r
tda_plot_ps(x, newpage = TRUE, cex = NULL, ...)
```

## Arguments

- x:

  a path, a TDA run, or the result of
  [`tda_read_ps`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md).

- newpage:

  start a new page before drawing.

- cex:

  character expansion for the text, relative to the size TDA asked for.
  The sizes in the file are PostScript points against a plot of perhaps
  80mm; drawn into a device several times that, they come out small, so
  the default scales them with the drawing, times
  `getOption("tdaR.ps.cex", 1)`.

- ...:

  passed to
  [`tda_read_ps`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md).

## Value

The parsed operations, invisibly.

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
[`tda_pl_regression()`](https://janmarvin.github.io/TDA/reference/tda_pl_regression.md),
[`tda_pl_scatter()`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md),
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md),
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
# a PostScript file with something in it: five points on a quadratic.
# psetup needs its ranges -- axes alone draw nothing, and tda_plot_ps
# then has nothing to render.
r <- tda_run(c("nvar(noc = 5, X = case, Y = case * case);",
               "psfile = p.ps;",
               "psetup(pxlen = 90, pylen = 50, pxa = 0,6, pya = 0,30);",
               "plxa(sc = 1); plya(sc = 10); plframe;",
               "plot = X,Y;", "psclose;"))
p <- tda_plot_ps(file.path(r$dir, "p.ps"))
```
