# Build a linear parameter constraint

[`tda_rate`](tda_rate.md) takes constraints in TDA's notation, where
parameters are numbered in the order they are estimated:
`"b3 - b10 = 0"` makes the third and tenth equal. This builds such a
string from a fitted model's coefficient names, so the numbering does
not have to be counted out by hand.

## Usage

``` r
tda_constrain(fit, ...)
```

## Arguments

- fit:

  a fit from [`tda_rate`](tda_rate.md), used for its coefficient names.

- ...:

  constraints, each a character vector of two coefficient names to be
  held equal, or a formula like `EDU ~ PRES`.

## Value

A character vector of constraints for `tda_rate`.

## See also

Other rate models: [`tda_control()`](tda_control.md),
[`tda_dple()`](tda_ltb.md), [`tda_rates()`](tda_rates.md),
[`tda_split()`](tda_split.md), [`tda_survivor()`](tda_survivor.md),
[`tda_transitions()`](tda_transitions.md),
[`vcov.tda_fit()`](tda_rate.md)

## Examples

``` r
d <- tda_rrdat(states = 4)
f <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, d)
tda_constrain(f, c("0->1: EDU", "0->2: EDU"))
#> [1] "b2 - b5 = 0"
```
