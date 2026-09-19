# Eigen decomposition of a general (possibly non-symmetric) matrix

Eigenvalues and eigenvectors of a square matrix that need not be
symmetric, so eigenvalues and eigenvectors may be complex – TDA's `mev`.
`t_mat.c`'s dispatch comment understates this command's arity
(`mev(A,ER,EI,EV)`, four arguments); what `m_mev` actually parses is
five (`mev(A,ER,EI,EVR,EVI)`), splitting the eigenvector matrix into
separate real and imaginary parts – confirmed by counting `m_getmat`
calls and matched against `eigen1`'s header, which does list five.
Column `i` of the result, `EVR[, i] + 1i * EVI[, i]`, is the right
eigenvector for eigenvalue `ER[i] + 1i * EI[i]`, i.e.
`A %*% v == lambda * v` – the same convention and ordering as base R's
`eigen(A)`, confirmed on two matrices, one with real eigenvalues and one
with an asymmetric complex pair (not simply plus/minus of each other,
which the first, simpler test matrix happened to have and which turned
out to mask a construction mistake in an earlier verification attempt –
recorded here because it is exactly the kind of false confirmation a
too-simple or too-symmetric test case can produce, the same lesson
[`tda_mevs`](https://janmarvin.github.io/TDA/reference/tda_mevs.md)
needed).

## Usage

``` r
tda_mev(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a square numeric matrix.

- dir:

  working directory.

## Value

A list: `values`, a complex vector of eigenvalues; `vectors`, a complex
matrix whose columns are the corresponding right eigenvectors.

## Details

TDA's `mev(A, ER, EI, EVR, EVI)` requires all four output operands (none
are optional); the wrapper supplies them and folds the real and
imaginary parts into complex `$values` and `$vectors`.

## See also

Other matrix algebra:
[`tda_mcel()`](https://janmarvin.github.io/TDA/reference/tda_mcel.md),
[`tda_mcent()`](https://janmarvin.github.io/TDA/reference/tda_mcent.md),
[`tda_mch()`](https://janmarvin.github.io/TDA/reference/tda_mch.md),
[`tda_mcross()`](https://janmarvin.github.io/TDA/reference/tda_mcross.md),
[`tda_mdiag()`](https://janmarvin.github.io/TDA/reference/tda_mdiag.md),
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
[`tda_msvd()`](https://janmarvin.github.io/TDA/reference/tda_msvd.md),
[`tda_mwvec()`](https://janmarvin.github.io/TDA/reference/tda_mwvec.md)

## Examples

``` r
tda_mev(matrix(c(0, -1, 1, 0), 2, 2))   # eigenvalues +-i
#> $values
#> [1] 0+1i 0-1i
#> 
#> $vectors
#>                 [,1]            [,2]
#> [1,] 7.395571e-32-1i 7.395571e-32+1i
#> [2,] 1.000000e+00+0i 1.000000e+00+0i
#> 
#> $er
#> [1] 0 0
#> 
#> $ei
#> [1]  1 -1
#> 
#> $evr
#>              [,1]         [,2]
#> [1,] 7.395571e-32 7.395571e-32
#> [2,] 1.000000e+00 1.000000e+00
#> 
#> $evi
#>      [,1] [,2]
#> [1,]   -1    1
#> [2,]    0    0
#> 
```
