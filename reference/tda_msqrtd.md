# Diagonal square root and its reciprocal

Given a square matrix `A`, returns a diagonal matrix (zero off-diagonal)
built from its diagonal – TDA's `msqrtd` (`diag(sqrt(diag(A)))`; refuses
a negative diagonal element) and `msqrti` (`diag(1/sqrt(diag(A)))`;
refuses an (almost) zero diagonal element). Off-diagonal entries of `A`
are ignored entirely, same as [`tda_minvd`](tda_minvs.md).

## Usage

``` r
tda_msqrtd(A, dir = tempfile("tda"))

tda_msqrti(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a square numeric matrix.

- dir:

  working directory.

## Value

A diagonal numeric matrix.

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
[`tda_mscal1()`](tda_mscal1.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_msqrtd(matrix(c(4, 1, 1, 9), 2))    # diag(2, 3)
#>      [,1] [,2]
#> [1,]    2    0
#> [2,]    0    3
tda_msqrti(matrix(c(4, 1, 1, 9), 2))    # diag(0.5, 1/3)
#>      [,1]      [,2]
#> [1,]  0.5 0.0000000
#> [2,]  0.0 0.3333333
```
