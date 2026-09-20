# Contour plot of a matrix

TDA's `plotcm` over a value matrix at the given levels.

## Usage

``` r
tda_plotcm(z, levels, ...)
```

## Arguments

- z:

  value matrix.

- levels:

  contour levels.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
# a bowl-shaped surface over an 8 x 8 grid: contours at three heights
z <- outer(1:8, 1:8, function(i, j) (i - 4)^2 + (j - 4)^2)
f <- tda_plotcm(z, levels = c(2, 6, 12))
tda_plot_ps(tda_read_ps(f))
```
