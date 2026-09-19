# Column permutations toward graph agreement

TDA's `gap` (CACM algorithm 548) on two integer-valued graphs, given as
matrices the way
[`tda_gqap`](https://janmarvin.github.io/TDA/reference/tda_gqap.md)
takes them.

## Usage

``` r
tda_gap_permute(g1, g2, ...)
```

## Arguments

- g1, g2:

  symmetric integer matrices of equal dimension.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

`permuted`, the permuted matrix, and `output`.

## Examples

``` r
f <- rbind(c(0, 3, 1), c(3, 0, 2), c(1, 2, 0))
d <- rbind(c(0, 1, 4), c(1, 0, 2), c(4, 2, 0))
r <- tda_gap_permute(f, d)
r$permuted
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    2    2
#> [3,]    3    3
```
