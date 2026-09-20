# Weighted mean of the later elements

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
[`tda_msvd()`](https://janmarvin.github.io/TDA/reference/tda_msvd.md)

## Examples

``` r
tda_mwvec(c(10, 20, 30, 40), c(1, 2, 1, 3))
#> [1] 31.66667 37.50000 40.00000 40.00000
```
