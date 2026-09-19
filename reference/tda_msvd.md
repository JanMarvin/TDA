# Singular value decomposition

`tda_msvd` returns just the singular values – TDA's `msvd`. `tda_msvd1`
returns the full decomposition `A = U %*% diag(d) %*% t(V)` – TDA's
`msvd1`, matching base R's [`svd()`](https://rdrr.io/r/base/svd.html)
convention exactly (singular vectors as columns of `u`/`v`), confirmed
by reconstructing `A` from the result. `A` must have at least as many
rows as columns for either. As with any singular vector, each column's
sign is arbitrary and need not match
[`svd()`](https://rdrr.io/r/base/svd.html)'s choice.

## Usage

``` r
tda_msvd(A, dir = tempfile("tda"))

tda_msvd1(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix with at least as many rows as columns.

- dir:

  working directory.

## Value

`tda_msvd` returns a plain numeric vector of singular values.
`tda_msvd1` returns a list: `d`, the singular values; `u`, `v`, the
singular vectors as columns.

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
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
A <- matrix(c(1, 3, 5, 2, 4, 7), 3, 2)
tda_msvd(A)     # == svd(A)$d
#> [1] 10.1914283  0.3671377
tda_msvd1(A)    # == svd(A), up to per-column sign
#> $d
#> [1] 10.1914283  0.3671377
#> 
#> $u
#>           [,1]        [,2]
#> [1,] 0.2167839  0.93892288
#> [2,] 0.4904541 -0.34146390
#> [3,] 0.8440732 -0.04273442
#> 
#> $v
#>           [,1]       [,2]
#> [1,] 0.5797531 -0.8147922
#> [2,] 0.8147922  0.5797531
#> 
```
