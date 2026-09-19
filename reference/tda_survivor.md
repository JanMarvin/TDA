# The survivor function as a single table

TDA writes one block per group and its examples select a row range out
of the flat file with `tsel = Case[16,,30]`. This returns the whole
thing as one long data frame instead, with a `group` column, so it can
be plotted, subset or joined without knowing how many blocks there were.

## Usage

``` r
tda_survivor(x, conf.int = 0.95)
```

## Arguments

- x:

  a fit from [`tda_ltb`](tda_ltb.md), [`tda_ple`](tda_ltb.md) or
  [`tda_km`](tda_ltb.md).

- conf.int:

  width of the confidence band to add as `lower` and `upper`, or `NULL`
  for none. TDA has no band of its own – it reports the standard error
  and leaves the band to the caller – so these are the usual normal
  limits, clamped to \[0, 1\].

## Value

A data frame with `group`, `time`, `survivor`, `std.err`, and where
available `density`, `rate`, `lower` and `upper`.

## Details

The same call works for [`tda_ltb`](tda_ltb.md) and
[`tda_ple`](tda_ltb.md): both name their time column `time` here,
whatever the underlying command called it, and both carry `survivor` and
`std.err`. Where the estimator provides them, `density` and `rate` come
too.

## See also

Other rate models: [`tda_constrain()`](tda_constrain.md),
[`tda_control()`](tda_control.md), [`tda_dple()`](tda_ltb.md),
[`tda_rates()`](tda_rates.md), [`tda_split()`](tda_split.md),
[`tda_transitions()`](tda_transitions.md),
[`vcov.tda_fit()`](tda_rate.md)

## Examples

``` r
d <- data.frame(t = c(4, 3, 1, 5, 8, 2), s = c(1, 1, 0, 1, 1, 1),
                g = c(1, 1, 1, 2, 2, 2))
f <- tda_ltb(Surv(t, s) ~ as.factor(g), d, tp = seq(0, 10, 2))
tda_survivor(f)
#>   group time  survivor   std.err   density      rate    lower    upper
#> 1     1    0 1.0000000 0.0000000 0.0000000 0.0000000 1.000000 1.000000
#> 2     1    2 1.0000000 0.0000000 0.2500000 0.3333333 1.000000 1.000000
#> 3     1    4 0.5000000 0.3535534 0.2500000 1.0000000 0.000000 1.000000
#> 4     2    0 1.0000000 0.0000000 0.0000000 0.0000000 1.000000 1.000000
#> 5     2    2 1.0000000 0.0000000 0.1666667 0.2000000 1.000000 1.000000
#> 6     2    4 0.6666667 0.2721655 0.1666667 0.3333333 0.133232 1.000000
#> 7     2    6 0.3333333 0.2721655 0.0000000 0.0000000 0.000000 0.866768
#> 8     2    8 0.3333333 0.2721655 0.1666667 1.0000000 0.000000 0.866768
```
