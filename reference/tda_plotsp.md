# Akima-smoothed curve through points

TDA's `plotsp`.

## Usage

``` r
tda_plotsp(x, y, ...)
```

## Arguments

- x, y:

  coordinates.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
# a smooth curve through five points that zigzag between 0 and 1
f <- tda_plotsp(x = c(0, 1, 2, 3, 4), y = c(0, 1, 0, 1, 0))
tda_plot_ps(tda_read_ps(f))
```
