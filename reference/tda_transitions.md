# Transitions covered by a fit

The origin and destination state pairs a fitted model covers. A
multi-state fit returns more than one row; a single-transition fit
returns one.

## Usage

``` r
tda_transitions(x)
```

## Arguments

- x:

  a fitted model.

## Value

A data frame with `Org` and `Des`, or `NULL`.

## See also

Other rate models: [`tda_constrain()`](tda_constrain.md),
[`tda_control()`](tda_control.md), [`tda_dple()`](tda_ltb.md),
[`tda_rates()`](tda_rates.md), [`tda_split()`](tda_split.md),
[`tda_survivor()`](tda_survivor.md), [`vcov.tda_fit()`](tda_rate.md)

## Examples

``` r
set.seed(9)
mn <- 60
m1 <- data.frame(ts = 0, tf = round(rexp(mn, 0.1), 1) + 0.5, org = 0,
                 des = sample(c(0, 1, 2), mn, TRUE, c(.25, .45, .30)),
                 x = round(rnorm(mn), 2))
mi <- which(m1$des == 1)
m2 <- data.frame(ts = m1$tf[mi],
                 tf = m1$tf[mi] + round(rexp(length(mi), 0.08), 1) + 0.5,
                 org = 1, des = sample(c(1, 2), length(mi), TRUE, c(.4, .6)),
                 x = m1$x[mi])
md <- rbind(m1, m2)
# three transitions: 0->1, 0->2 (censored at 2), 1->2
mf <- tda_rate(Surv(ts, tf, org, des) ~ x, md, model = "exponential")
tda_transitions(mf)
#>   Org Des
#> 1   0   1
#> 2   0   2
#> 3   1   2
```
