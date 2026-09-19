# A variable's value by category, for tda_qreg

Marks a set of columns, used inside a [`tda_qreg`](tda_qreg.md) formula,
as the same underlying variable's value at each category of the response
– price for choice 1, choice 2, choice 3, say, one coefficient shared
across all of them – rather than the formula's ordinary
one-column-per-predictor, shared-across-categories shape. TDA's
`"(Z1,Z2,...)"` grouped variable-list syntax (what its source calls "Z
variables" or "generic variables"), reached through R's `Surv()`-style
formula convention (from the `survival` package – not a dependency of
this one, just the same idea) instead of a separate argument.

## Usage

``` r
lvl(...)
```

## Arguments

- ...:

  one column per category of the response, in category order.

## Value

A matrix, for
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)'s benefit –
not meaningful on its own.

## Details

`lvl()` is only meaningful inside a `tda_qreg` formula – it is not a
real transformation, just a marker `tda_qreg` looks for and removes
before fitting, the same way `Surv()` is not itself computed by
`coxph()`, only recognised by it.

## Examples

``` r
# Daganzo's classic mode-choice example (examples/exam/qr4.cf in TDA's
# own sources): commute time by car (Z1), bus (Z2), and train (Z3), one
# coefficient for time itself, shared across all three modes
# (seeded synthetic version: choices generated from utility
# -0.3 * time + Gumbel noise, so the estimate recovers about -0.3)
set.seed(21)
n <- 80
d <- data.frame(Z1 = runif(n, 5, 30), Z2 = runif(n, 5, 30),
                Z3 = runif(n, 5, 30))
u <- -0.3 * as.matrix(d) + matrix(-log(-log(runif(3 * n))), n)
d$Y <- max.col(u)
tda_qreg(Y ~ lvl(Z1, Z2, Z3), d, model = "multinomial_logit",
         nq = 3, intercept = FALSE)$estimates
#>   Idx Cat Term Variable      Coeff      Error   C/Error Signif
#> 1   1  NA   Z1   Lvl1_1 -0.2783749 0.04944514 -5.629974      1
```
