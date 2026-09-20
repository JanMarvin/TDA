# Average nearby endpoints for interval-censored data

For each interval `i`, averages the lower endpoints of every interval
(itself included) whose lower endpoint falls within
`[lower[i], upper[i]]`, and separately averages the upper endpoints
falling in that same range – TDA's `midf3`. A local smoothing of nearby
endpoints, not a distribution estimate like
[`tda_midf`](https://janmarvin.github.io/TDA/reference/tda_midf.md)/[`tda_midf1`](https://janmarvin.github.io/TDA/reference/tda_midf1.md)/
[`tda_midf2`](https://janmarvin.github.io/TDA/reference/tda_midf2.md).

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

Other matrix algebra:
[`tda_mcel()`](https://janmarvin.github.io/TDA/reference/tda_mcel.md),
[`tda_mcent()`](https://janmarvin.github.io/TDA/reference/tda_mcent.md),
[`tda_mch()`](https://janmarvin.github.io/TDA/reference/tda_mch.md),
[`tda_mcross()`](https://janmarvin.github.io/TDA/reference/tda_mcross.md),
[`tda_mdiag()`](https://janmarvin.github.io/TDA/reference/tda_mdiag.md),
[`tda_mev()`](https://janmarvin.github.io/TDA/reference/tda_mev.md),
[`tda_mevs()`](https://janmarvin.github.io/TDA/reference/tda_mevs.md),
[`tda_midf()`](https://janmarvin.github.io/TDA/reference/tda_midf.md),
[`tda_midf1()`](https://janmarvin.github.io/TDA/reference/tda_midf1.md),
[`tda_midf2()`](https://janmarvin.github.io/TDA/reference/tda_midf2.md),
[`tda_minvs()`](https://janmarvin.github.io/TDA/reference/tda_minvs.md),
[`tda_mkmet()`](https://janmarvin.github.io/TDA/reference/tda_mkmet.md),
[`tda_mkp()`](https://janmarvin.github.io/TDA/reference/tda_mkp.md),
[`tda_mmul()`](https://janmarvin.github.io/TDA/reference/tda_mmul.md),
[`tda_mnc()`](https://janmarvin.github.io/TDA/reference/tda_mnc.md),
[`tda_mnrow()`](https://janmarvin.github.io/TDA/reference/tda_mnrow.md),
[`tda_mpbl()`](https://janmarvin.github.io/TDA/reference/tda_mpbl.md),
[`tda_mpfit()`](https://janmarvin.github.io/TDA/reference/tda_mpfit.md),
[`tda_mpinv()`](https://janmarvin.github.io/TDA/reference/tda_mpinv.md),
[`tda_mpit()`](https://janmarvin.github.io/TDA/reference/tda_mpit.md),
[`tda_mple()`](https://janmarvin.github.io/TDA/reference/tda_mple.md),
[`tda_mpz()`](https://janmarvin.github.io/TDA/reference/tda_mpz.md),
[`tda_mscal1()`](https://janmarvin.github.io/TDA/reference/tda_mscal1.md),
[`tda_msqrtd()`](https://janmarvin.github.io/TDA/reference/tda_msqrtd.md),
[`tda_msvd()`](https://janmarvin.github.io/TDA/reference/tda_msvd.md),
[`tda_mwvec()`](https://janmarvin.github.io/TDA/reference/tda_mwvec.md)

## Examples

``` r
lo <- c(0, 1, 2, 3, 1)
up <- c(2, 3, 4, 5, 2)
# for each interval, the mean lower and upper endpoint of the
# intervals overlapping it
tda_midf3(lo, up)
#> $lower_avg
#> [1] 1.000000 1.750000 2.500000 3.000000 1.333333
#> 
#> $upper_avg
#> [1] 2.000000 2.333333 2.750000 4.000000 2.000000
#> 
```
