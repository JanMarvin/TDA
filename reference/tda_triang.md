# Triangulation of points

TDA's `triang` over an x,y point set.

## Usage

``` r
tda_triang(x, y, ...)
```

## Arguments

- x, y:

  coordinate vectors.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the output file's lines (the triangle list).

## Examples

``` r
# a Delaunay triangulation of five points: the edges of the triangles
tda_triang(x = c(0, 2, 1, 3, 1), y = c(0, 0, 1, 1, 2))
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    1    5
#> [4,]    2    4
#> [5,]    2    3
#> [6,]    3    4
#> [7,]    3    5
#> [8,]    4    5
```
