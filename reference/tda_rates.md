# The estimated rate as a single table

The companion to
[`tda_survivor`](https://janmarvin.github.io/TDA/reference/tda_survivor.md):
the rate a model implies, with the same column names whichever model
produced it. `prate =` at fit time is what asks TDA to compute it.

## Usage

``` r
tda_rates(x)
```

## Arguments

- x:

  a fit from
  [`tda_rate`](https://janmarvin.github.io/TDA/reference/tda_rate.md) or
  [`tda_coxph`](https://janmarvin.github.io/TDA/reference/tda_rate.md).

## Value

A data frame with `group`, `time`, `rate`, `survivor` and `density`, or
`NULL` if the fit carries no rate table. For a Cox fit `rate` is the
baseline rate at each event time and `cumrate` the cumulative baseline
rate beside it (already cumulative: TDA's CumRate column).

## See also

Other rate models:
[`tda_constrain()`](https://janmarvin.github.io/TDA/reference/tda_constrain.md),
[`tda_control()`](https://janmarvin.github.io/TDA/reference/tda_control.md),
[`tda_dple()`](https://janmarvin.github.io/TDA/reference/tda_ltb.md),
[`tda_split()`](https://janmarvin.github.io/TDA/reference/tda_split.md),
[`tda_survivor()`](https://janmarvin.github.io/TDA/reference/tda_survivor.md),
[`tda_transitions()`](https://janmarvin.github.io/TDA/reference/tda_transitions.md),
[`vcov.tda_fit()`](https://janmarvin.github.io/TDA/reference/tda_rate.md)

## Examples

``` r
d <- tda_rrdat()
f <- tda_rate(Surv(TFP, DES) ~ 1, d, prate = "0(10)300")
head(tda_rates(f))
#>   group time       rate  survivor     density
#> 1     0    0 0.01123044 1.0000000 0.011230445
#> 2     0   10 0.01123044 0.8937721 0.010037458
#> 3     0   20 0.01123044 0.7988286 0.008971200
#> 4     0   30 0.01123044 0.7139707 0.008018209
#> 5     0   40 0.01123044 0.6381271 0.007166451
#> 6     0   50 0.01123044 0.5703402 0.006405174
```
