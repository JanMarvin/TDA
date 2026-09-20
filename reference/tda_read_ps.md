# Read TDA's PostScript output

Parses the PostScript a plot command writes into a table of drawing
operations: line segments, polygons, symbols and text, with the graphics
state that applied to each.

## Usage

``` r
tda_read_ps(file, which = NULL)
```

## Arguments

- file:

  a `.ps` file written by TDA, or a run whose directory contains one.

- which:

  when `file` is a run, the name of the PostScript file.

## Value

A list with `ops`, the drawing operations, and `bbox`, the bounding box
TDA declared. Each operation carries the plot command that produced it,
taken from the `%#` annotations TDA writes into the stream.

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
[`tda_plot_ps()`](https://janmarvin.github.io/TDA/reference/tda_plot_ps.md),
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md),
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md)

## Examples

``` r
d <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
p <- tda_ps(d, xlim = c(0, 6), ylim = c(0, 10))
p <- tda_pl_axes(p)
p <- tda_pl_lines(p, "x", "y")
pdf(NULL)  # plot() runs the session and draws it; capture the run
p <- plot(p)
dev.off()
#> agg_record_1f87755ee3b6 
#>                       2 
ps <- tda_read_ps(p$run, which = p$file)
ps$bbox
#> [1] 145 455 438 663
ps$ops[[3]]$op
#> [1] "lines"
```
