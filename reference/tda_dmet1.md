# Distances made metric

TDA's `dmet1`: adjusts a distance matrix (given as an undirected valued
graph) toward metricity.

## Usage

``` r
tda_dmet1(d, tolerance = 1e-04, max_iterations = 20, ...)
```

## Arguments

- d:

  symmetric distance matrix.

- tolerance, max_iterations:

  deviation tolerance and cap.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the output file's lines.

## Examples

``` r
# 5 breaks the triangle inequality (1 + 2 < 5), so the
# metricized matrix shortens it to 3
d <- rbind(c(0, 1, 5), c(1, 0, 2), c(5, 2, 0))
r <- tda_dmet1(d)
r$modified
#>      [,1] [,2] [,3]
#> [1,]    0    1    3
#> [2,]    1    0    2
#> [3,]    3    2    0
```
