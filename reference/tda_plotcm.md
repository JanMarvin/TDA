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
z <- outer(1:8, 1:8, function(i, j) (i - 4)^2 + (j - 4)^2)
r <- tda_plotcm(z, levels = c(2, 6, 12))
file.exists(grep("[.]ps$", r, value = TRUE))   # the finished plot
#> [1] TRUE
```
