# All rank orders of a given size

`tda_rank_orders` enumerates every rank order of `size` objects,
including all ways of tying them. For three objects there are 13: six
strict orders, six with one pair tied, and one with all three tied.

## Usage

``` r
tda_rank_orders(size = 2, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- size:

  number of objects to rank, 2 or more.

- options:

  a named list of further TDA options, passed through. `m=` selects the
  enumeration variant.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame with one row per rank order: `ties`, the number of tie
groups, then one column per object (`r1`, `r2`, ...) giving its rank.
`attr(x, "run")` carries the run.

## Details

Useful as the sample space for a rank-order model – the set of outcomes
a distribution over rankings has to sum over – which is otherwise fiddly
to generate because of the ties.

## Examples

``` r
r <- tda_rank_orders(3)
nrow(r)          # 13
#> [1] 13
head(r)
#>   ties r1 r2 r3
#> 1    1  1  1  1
#> 2    2  1  2  2
#> 3    2  2  1  1
#> 4    2  2  1  2
#> 5    2  1  2  1
#> 6    2  2  2  1
table(r$ties)    # 1 all tied, 6 one pair tied, 6 strict
#> 
#> 1 2 3 
#> 1 6 6 
```
