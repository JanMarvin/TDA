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

Other plotting: [`plot.tda_ple()`](plot.tda_ple.md),
[`tda_check_ps()`](tda_check_ps.md), [`tda_pl()`](tda_pl.md),
[`tda_pl_arc()`](tda_pl_arc.md), [`tda_pl_axis()`](tda_pl_axis.md),
[`tda_pl_graph()`](tda_pl_graph.md), [`tda_pl_hist`](tda_pl_hist.md),
[`tda_pl_panel()`](tda_pl_panel.md),
[`tda_pl_regression()`](tda_pl_regression.md),
[`tda_pl_scatter()`](tda_pl_scatter.md),
[`tda_plot_ps()`](tda_plot_ps.md), [`tda_ps()`](tda_ps.md),
[`tda_ps3()`](tda_ps3.md)

## Examples

``` r
d <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
p <- tda_ps(d, xlim = c(0, 6), ylim = c(0, 10))
p <- tda_pl_axes(p)
p <- tda_pl_lines(p, "x", "y")
pdf(NULL)  # plot() runs the session and draws it; capture the run
p <- plot(p)
dev.off()
#> agg_record_210543aae572 
#>                       2 
ps <- tda_read_ps(p$run, which = p$file)
ps$bbox
#> [1] 145 455 438 663
ps$ops[[3]]$op
#> [1] "lines"
```
