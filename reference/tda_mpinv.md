# Invert a permutation

TDA's `mpinv`: given a permutation `p` of `1:n`, returns `q` such that
`q[p[i]] == i` – matching base R's `order(p)`.

## Usage

``` r
tda_mpinv(p, dir = tempfile("tda"))
```

## Arguments

- p:

  an integer permutation of `1:length(p)`.

- dir:

  working directory.

## Value

An integer vector, the inverse permutation.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_minvs()`](tda_minvs.md), [`tda_mkmet()`](tda_mkmet.md),
[`tda_mkp()`](tda_mkp.md), [`tda_mmul()`](tda_mmul.md),
[`tda_mnc()`](tda_mnc.md), [`tda_mnrow()`](tda_mnrow.md),
[`tda_mpbl()`](tda_mpbl.md), [`tda_mpfit()`](tda_mpfit.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_mpinv(c(3, 1, 4, 2))   # == order(c(3, 1, 4, 2))
#> [1] 2 4 1 3
```
