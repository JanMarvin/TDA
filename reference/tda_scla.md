# Separable clusters of a valued graph

TDA's `scla` on a distance matrix.

## Usage

``` r
tda_scla(d, ...)
```

## Arguments

- d:

  symmetric distance matrix.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

`clusters`, the rows TDA wrote (one numeric vector per line of the
separable-cluster listing), and `output`.

## Examples

``` r
d <- as.matrix(dist(c(0, 0.1, 5, 5.1)))
r <- tda_scla(d)
r$clusters
#> NULL
```
