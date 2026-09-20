# Hasse diagram of multivariate patterns

TDA's `ghd`: the rows of `x` are (assumed unique) patterns; the command
derives the dominance order's Hasse diagram.

## Usage

``` r
tda_ghd(x, ...)
```

## Arguments

- x:

  matrix or data frame of patterns.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the edge list lines of the diagram.

## Examples

``` r
# five patterns over two criteria; a row dominates another when it
# is at least as large in both.  The Hasse diagram keeps only the
# covering relations: (0,0) < (1,0), (0,1) < (1,1) < (2,1)
tda_ghd(rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(2, 1)))
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    1    3
#> [3,]    2    4
#> [4,]    3    4
#> [5,]    4    5
```
