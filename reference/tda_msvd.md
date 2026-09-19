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
[`tda_mwvec()`](https://janmarvin.github.io/TDA/reference/tda_mwvec.md)

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
