# Projection of proximities

TDA's `gpro` on a proximity matrix given as an undirected valued graph.

## Usage

``` r
tda_gpro(d, max_iterations = 50, ...)
```

## Arguments

- d:

  symmetric proximity matrix.

- max_iterations:

  iteration cap.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

`projection`, the projected coordinates, and `output`.

## Examples

``` r
r <- tda_gpro(as.matrix(dist(cbind(c(0, 0, 3, 3), c(0, 2, 0, 2)))))
r$projection
#>      [,1]      [,2]       [,3]
#> [1,]    1 0.0000000  0.0000000
#> [2,]    2 1.9688721 -0.3517456
#> [3,]    3 0.5275889  2.9531913
#> [4,]    4 2.4963804  2.6015133
```
