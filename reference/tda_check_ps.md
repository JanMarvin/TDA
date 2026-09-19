# Audit a parse against the PostScript it came from

Counts the drawing operators in the file and the operations the parser
produced, and reports anything unaccounted for. Reading a plot back is
easy to get half right – an operator that is silently ignored looks
exactly like one that is not there – so this compares the two directly
rather than leaving it to the eye.

## Usage

``` r
tda_check_ps(file, ...)
```

## Arguments

- file:

  a `.ps` file, a TDA run, or a parse from
  [`tda_read_ps`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md).

- ...:

  passed to
  [`tda_read_ps`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md).

## Value

A data frame with one row per operator: how often it appears in the
file, whether the parser handles it, and how many operations resulted.

## See also

Other plotting:
[`plot.tda_ple()`](https://janmarvin.github.io/TDA/reference/plot.tda_ple.md),
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
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
p <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
p <- tda_pl_circle(p, at = c(5, 5), r = 2)
tda_check_ps(tda_ps_file(p))
#>        operator in_file handled ignorable unaccounted
#> 1             m       2    TRUE     FALSE       FALSE
#> 2             l       3    TRUE     FALSE       FALSE
#> 5        stroke       2    TRUE     FALSE       FALSE
#> 7     closepath       1    TRUE     FALSE       FALSE
#> 8       newpath       1   FALSE      TRUE       FALSE
#> 9          clip       1   FALSE      TRUE       FALSE
#> 10        gsave       1    TRUE     FALSE       FALSE
#> 11     grestore       1    TRUE     FALSE       FALSE
#> 12          arc       1    TRUE     FALSE       FALSE
#> 18 setlinewidth       1    TRUE     FALSE       FALSE
#> 19      setdash       1    TRUE     FALSE       FALSE
#> 20    translate       1    TRUE     FALSE       FALSE
#> 21       rotate       1    TRUE     FALSE       FALSE
#> 22        scale       1    TRUE     FALSE       FALSE
#> 25     showpage       1   FALSE      TRUE       FALSE
```
