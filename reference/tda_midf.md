# Distribution function bounds and estimate for interval-censored data

Given a sample of intervals `[lower_i, upper_i]` (each interval
containing an unobserved exact value), estimates the distribution
function of that unobserved value – TDA's `midf` (undocumented beyond
"calculates distribution function DM, and lower/upper bounds" in its own
header, which also understates the command's arity: it lists three
outputs, `midf(XL,XU,DL,DU,DM)`, but `m_midf` actually returns four,
`midf(XL,XU,PT,DL,DU,DM)` – the breakpoints themselves, `PT`, are a
result too, confirmed by counting `m_getmat` calls). Evaluated at every
distinct interval endpoint (the “induced partition” `breakpoints`):
`lower_bound` counts only intervals entirely below that point (a sure
lower bound on the true CDF); `upper_bound` counts every interval that
could be below it (a sure upper bound); `cdf` is a point estimate in
between, assuming each interval's unobserved value is uniformly
distributed within it. Every interval must have positive width
(`lower < upper`, swapped automatically if given in the other order).
Confirmed cell-for-cell by hand on a 3-interval example.

## Usage

``` r
tda_midf(lower, upper, dir = tempfile("tda"))
```

## Arguments

- lower, upper:

  numeric vectors, the same length, with `lower[i] != upper[i]` for
  every `i`.

- dir:

  working directory.

## Value

A list: `breakpoints`, the sorted distinct interval endpoints;
`lower_bound`, `upper_bound`, `cdf`, each the same length as
`breakpoints`.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf1()`](tda_midf1.md), [`tda_midf2()`](tda_midf2.md),
[`tda_midf3()`](tda_midf3.md), [`tda_minvs()`](tda_minvs.md),
[`tda_mkmet()`](tda_mkmet.md), [`tda_mkp()`](tda_mkp.md),
[`tda_mmul()`](tda_mmul.md), [`tda_mnc()`](tda_mnc.md),
[`tda_mnrow()`](tda_mnrow.md), [`tda_mpbl()`](tda_mpbl.md),
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_midf(c(1, 2, 0), c(3, 4, 2))
#> $breakpoints
#> [1] 0 1 2 3 4
#> 
#> $lower_bound
#> [1] 0.0000000 0.0000000 0.3333333 0.6666667 1.0000000
#> 
#> $upper_bound
#> [1] 0.3333333 0.6666667 1.0000000 1.0000000 1.0000000
#> 
#> $cdf
#> [1] 0.0000000 0.1666667 0.5000000 0.8333333 1.0000000
#> 
```
