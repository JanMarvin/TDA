# The assignment problem

TDA's `gap` on a cost matrix: the permutation that assigns each row to a
column at minimal total cost (Carpaneto and Toth's algorithm), the
manual's 7.3.1.1. The matrix goes to TDA as a full matrix graph
(`gdd(opt=7)`), so it need not be symmetric.

## Usage

``` r
tda_gap(cost, ...)
```

## Arguments

- cost:

  a square numeric cost matrix.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

`assignment`, the permutation as a two-column matrix (i, p(i));
`permuted`, the cost matrix with its rows permuted accordingly; and
`output`.

## Examples

``` r
cost <- rbind(c(60, 0, 0, 76, 0, 0), c(0, 40, 18, 0, 60, 24),
              c(60, 16, 2, 4, 0, 40), c(0, 27, 18, 3, 55, 75),
              c(0, 40, 62, 16, 11, 53), c(28, 4, 10, 84, 0, 16))
tda_gap(cost)$assignment
#>      [,1] [,2]
#> [1,]    1    6
#> [2,]    2    4
#> [3,]    3    3
#> [4,]    4    1
#> [5,]    5    5
#> [6,]    6    2
```
