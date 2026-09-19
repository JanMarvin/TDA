# Unfolding of preference data

TDA's `unf`: each row ranks the alternatives (larger is preferred).

## Usage

``` r
tda_unf(x, ...)
```

## Arguments

- x:

  matrix of preference values.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the printed output.

## Examples

``` r
r <- tda_unf(rbind(c(3, 2, 1), c(1, 3, 2), c(2, 3, 1)))
cat(tail(tda_payload(r), 4), sep = "\n")   # the best permutation
#> rr: 1 3 2 
#> UNFMax=9
#> Best value (sum of taus): 3
#> Best permutation: 1 2 3 
```
