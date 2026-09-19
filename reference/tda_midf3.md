# Average nearby endpoints for interval-censored data

For each interval `i`, averages the lower endpoints of every interval
(itself included) whose lower endpoint falls within
`[lower[i], upper[i]]`, and separately averages the upper endpoints
falling in that same range – TDA's `midf3`. A local smoothing of nearby
endpoints, not a distribution estimate like
[`tda_midf`](tda_midf.md)/[`tda_midf1`](tda_midf1.md)/
[`tda_midf2`](tda_midf2.md). Confirmed by hand on a 3-interval example.

## Usage

``` r
tda_midf3(lower, upper, dir = tempfile("tda"))
```

## Arguments

- lower, upper:

  numeric vectors, the same length.

- dir:

  working directory.

## Value

A list, both the same length as `lower`, in the original order:
`lower_avg`, the average nearby lower endpoint; `upper_avg`, the average
nearby upper endpoint.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_minvs()`](tda_minvs.md),
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
tda_midf3(c(1, 2, 0), c(3, 4, 2))
#> $lower_avg
#> [1] 1.5 2.0 1.0
#> 
#> $upper_avg
#> [1] 2.5 3.0 2.0
#> 
```
