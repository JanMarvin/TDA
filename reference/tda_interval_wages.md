# Artificial interval-valued wage data

A small, reproducible dataset built for trying the whole interval family
on one realistic shape: n = 40 people, schooling asked and answered in
whole years (a point value), wages reported in brackets that widen with
income, with a top-coded highest bracket – the classic survey mix of
exact and interval-valued variables. Because the regressor is exact,
`tda_ivreg`'s exact method certifies sharp slope bounds instantly here;
widen the schooling values into intervals yourself to watch
certification get hard.

## Usage

``` r
tda_interval_wages()
```

## Value

A data frame with columns `wage_lo`, `wage_hi`, `school_lo`,
`school_hi`.

## Details

The generating truth is `log(wage) = 0.8 + 0.11 school + e`, which
produces a wage-level slope around 1.6 over this range; any interval
method's answer should be judged against the data's coarseness, not
against that number.

Works across the family: [`tda_imean`](tda_imean.md),
[`tda_ivariance`](tda_imean.md) (with `$sd`), [`tda_icov`](tda_icov.md),
[`tda_icorr`](tda_icov.md) (which honestly warns that it cannot certify
at this size), [`tda_ivreg`](tda_ivreg.md) and
[`tda_ilsreg`](tda_ilsreg.md).

## Examples

``` r
d <- tda_interval_wages()
tda_imean(~ iv(wage_lo, wage_hi), d)
#> Call: tda_imean(~iv(wage_lo, wage_hi), d)
#> 
#> Cases: 40 
#> Bounds: [7.45, 11.7]
# \donttest{
tda_ivreg(iv(wage_lo, wage_hi) ~ iv(school_lo, school_hi), d,
          method = "exact")
#> Call: tda_ivreg(iv(wage_lo, wage_hi) ~ iv(school_lo, school_hi), d, 
#>     method = "exact")
#> 
#> Cases: 40 
#> Mean of the response: [9.425, 9.425] 
#> Mean of the regressor: [12.725, 12.725] 
#> Variance of the regressor: [7.649375, 7.649375] 
#> 
#> Bounds on the slope: [0.319225, 1.754964]  (exact, certified)
# }
```
