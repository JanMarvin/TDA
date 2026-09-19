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

  passed to [`tda_run`](tda_run.md).

## Value

the output file's lines (the triangle list).

## Examples

``` r
tda_triang(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1))
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    2    4
#> [4,]    2    3
#> [5,]    3    4
```
