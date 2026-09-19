# Split episodes at a time-varying covariate

A covariate that changes during an episode is handled by cutting the
episode in two at the moment it changes, so that each piece has one
value. TDA does this with the `split` option of `edef`, and
[`tda_episodes`](tda_episodes.md) and [`tda_rate`](tda_rate.md) take a
`split` argument that reaches it.

## Usage

``` r
tda_split(
  formula,
  data,
  at = NULL,
  grid = NULL,
  vars = NULL,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  the usual `Surv()` formula.

- data:

  a data frame.

- at:

  one or more variables holding a split time per episode, as names or as
  a vector, matrix or data frame – one column per cut.

- grid:

  fixed times at which every episode is cut, as a numeric vector.
  Combines with `at`.

- vars:

  further variables to carry into the split data.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame of split episodes, one row per piece.

## Details

`at` names a variable holding, for each episode, the time at which the
change happens – measured on the same axis as the episode's start and
end. An episode whose covariate never changes is given a value beyond
its end, conventionally a large number, and is not cut.

`edef` cuts once per split variable, so several are needed to cut an
episode more than once. `grid` is the common case said directly: a set
of fixed times every episode is cut at, one split variable per time.
`grid = seq(60, 420, by = 60)` is how TDA's `rrdat.d60` was built, 600
job episodes becoming 1021 pieces of at most five years.

## See also

Other rate models: [`tda_constrain()`](tda_constrain.md),
[`tda_control()`](tda_control.md), [`tda_dple()`](tda_ltb.md),
[`tda_rates()`](tda_rates.md), [`tda_survivor()`](tda_survivor.md),
[`tda_transitions()`](tda_transitions.md),
[`vcov.tda_fit()`](tda_rate.md)

## Examples

``` r
d <- tda_rrdat()
# marriage during the job episode, or beyond its end if never
d$MarrDate <- ifelse(d$TMAR <= 0, 10000, d$TMAR - d$TStart)
head(tda_split(Surv(TFP, DES) ~ EDU, d, at = "MarrDate"))
#>   episode case subsample transition org des  ts  tf EDU MarrDate
#> 1       1    1         2          1   0   0   0 124  17      124
#> 2       1    1         2          2   0   0 124 428  17      124
#> 3       2    1         1          1   0   1   0  46  10      169
#> 4       3    1         1          1   0   1   0  34  10      123
#> 5       4    1         2          1   0   0   0  89  10       89
#> 6       4    1         2          2   0   1  89 220  10       89
# every five years, however long the episode
nrow(tda_split(Surv(TFP, DES) ~ EDU, d, grid = seq(60, 420, by = 60)))
#> [1] 1021
```
