# Rank-order data utilities

TDA's `rod`: rank orders per record, with options for standardization,
distances, medians, and tabulation.

## Usage

``` r
tda_rod(
  x,
  statistic = c("standard", "distances", "central", "graph", "graph5", "graph6",
    "between", "table"),
  ...
)
```

## Arguments

- x:

  matrix of rank values, one order per row.

- statistic:

  one of "standard", "distances", "central", "graph", "graph5",
  "graph6", "between", "table".

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the rank orders as a data frame in `$orders` (plus a count column for
the table statistic), with the full protocol in `$output`.

## Examples

``` r
r <- tda_rod(rbind(c(1, 2, 3), c(2, 3, 1), c(3, 1, 2)))
r$orders
#>   r1 r2 r3
#> 1  1  2  3
#> 2  2  3  1
#> 3  3  1  2
```
