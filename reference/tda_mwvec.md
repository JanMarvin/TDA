# Backward weighted running mean

For each position `i`, averages `a[j]` weighted by `weights[j]` over
every later position (`tda_mwvec`: every `j > i`; `tda_mwvec1`: every
`j` with `t[j] > t[i]`, letting an arbitrary order vector `t` stand in
for position – `tda_mwvec` is the special case `t = seq_along(a)`) –
TDA's `mwvec`/`mwvec1`. Falls back to `a[i]` unchanged wherever the
weights beyond `i` sum to zero (including the last position, which
always has nothing after it). Confirmed against a hand-written R
equivalent on a random instance; no base-R builtin matches this
directly.

## Usage

``` r
tda_mwvec(a, weights, dir = tempfile("tda"))

tda_mwvec1(a, weights, t, dir = tempfile("tda"))
```

## Arguments

- a, weights:

  numeric vectors of the same length.

- dir:

  working directory.

- t:

  for `tda_mwvec1`, a numeric order vector the same length as `a`; ties
  are included on neither side (only strictly later positions count).

## Value

A numeric vector, the same length as `a`.

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
[`tda_mpinv()`](tda_mpinv.md), [`tda_mpit()`](tda_mpit.md),
[`tda_mple()`](tda_mple.md), [`tda_mpz()`](tda_mpz.md),
[`tda_mscal1()`](tda_mscal1.md), [`tda_msqrtd()`](tda_msqrtd.md),
[`tda_msvd()`](tda_msvd.md)

## Examples

``` r
tda_mwvec(c(10, 20, 30, 40), c(1, 2, 1, 3))
#> [1] 31.66667 37.50000 40.00000 40.00000
```
