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

  passed to [`tda_run`](tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
r <- tda_plotsp(x = c(0, 1, 2, 3, 4), y = c(0, 1, 0, 1, 0))
file.exists(grep("[.]ps$", r, value = TRUE))   # the finished plot
#> [1] TRUE
```
