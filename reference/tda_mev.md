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
with a asymmetric complex pair (not simply plus/minus of each other,
which the first, simpler test matrix happened to have and which turned
out to mask a construction mistake in an earlier verification attempt –
recorded here because it is exactly the kind of false confirmation a
too-simple or too-symmetric test case can produce, the same lesson
[`tda_mevs`](tda_mevs.md) needed).

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

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mevs()`](tda_mevs.md), [`tda_midf()`](tda_midf.md),
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
