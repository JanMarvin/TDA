# Distribution function bounds and estimate for interval-censored data

Given a sample of intervals `[lower_i, upper_i]` (each interval
containing an unobserved exact value), estimates the distribution
function of that unobserved value – TDA's `midf` (undocumented beyond
"calculates distribution function DM, and lower/upper bounds" in its own
header, which also understates the command's arity: it lists three
outputs, `midf(XL,XU,DL,DU,DM)`, but `m_midf` actually returns four,
`midf(XL,XU,PT,DL,DU,DM)` – the breakpoints themselves, `PT`, are a
result too). Evaluated at every distinct interval endpoint (the “induced
partition” `breakpoints`): `lower_bound` counts only intervals entirely
below that point (a sure lower bound on the true CDF); `upper_bound`
counts every interval that could be below it (a sure upper bound); `cdf`
is a point estimate in between, assuming each interval's unobserved
value is uniformly distributed within it. Every interval must have
positive width (`lower < upper`, swapped automatically if given in the
other order). Confirmed cell-for-cell by hand on a 3-interval example.

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

Other matrix algebra:
[`tda_mcel()`](https://janmarvin.github.io/TDA/reference/tda_mcel.md),
[`tda_mcent()`](https://janmarvin.github.io/TDA/reference/tda_mcent.md),
[`tda_mch()`](https://janmarvin.github.io/TDA/reference/tda_mch.md),
[`tda_mcross()`](https://janmarvin.github.io/TDA/reference/tda_mcross.md),
[`tda_mdiag()`](https://janmarvin.github.io/TDA/reference/tda_mdiag.md),
[`tda_mev()`](https://janmarvin.github.io/TDA/reference/tda_mev.md),
[`tda_mevs()`](https://janmarvin.github.io/TDA/reference/tda_mevs.md),
[`tda_midf1()`](https://janmarvin.github.io/TDA/reference/tda_midf1.md),
[`tda_midf2()`](https://janmarvin.github.io/TDA/reference/tda_midf2.md),
[`tda_midf3()`](https://janmarvin.github.io/TDA/reference/tda_midf3.md),
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
# five observations known only to intervals: [0,2], [1,3], [2,4], [3,5], [1,2]
lo <- c(0, 1, 2, 3, 1)
up <- c(2, 3, 4, 5, 2)
# at each endpoint: the share of intervals surely below it (lower
# bound), possibly below it (upper bound), and the estimate between
tda_midf(lo, up)
#> $breakpoints
#> [1] 0 1 2 3 4 5
#> 
#> $lower_bound
#> [1] 0.0 0.0 0.4 0.6 0.8 1.0
#> 
#> $upper_bound
#> [1] 0.2 0.6 0.8 1.0 1.0 1.0
#> 
#> $cdf
#> [1] 0.0 0.1 0.5 0.7 0.9 1.0
#> 
```
