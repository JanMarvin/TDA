# Matrix size and norms

Scalar summaries of a matrix – TDA's `mnrow`, `mncol`, `mnorm` (max
absolute element), `mnorm1` (sum of absolute elements), `mnorm2`
(Frobenius/Euclidean norm, `sqrt(sum(A^2))`), and `mtrace`
(`sum(diag(A))`, using `min(nrow(A), ncol(A))` diagonal elements for a
non-square `A`, same as base R's `sum(diag(A))`). Each returns a plain
number, not a 1x1 matrix.

## Usage

``` r
tda_mnrow(A, dir = tempfile("tda"))

tda_mncol(A, dir = tempfile("tda"))

tda_mnorm(A, dir = tempfile("tda"))

tda_mnorm1(A, dir = tempfile("tda"))

tda_mnorm2(A, dir = tempfile("tda"))

tda_mtrace(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix.

- dir:

  working directory.

## Value

A single number.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_minvs()`](tda_minvs.md), [`tda_mkmet()`](tda_mkmet.md),
[`tda_mkp()`](tda_mkp.md), [`tda_mmul()`](tda_mmul.md),
[`tda_mnc()`](tda_mnc.md), [`tda_mpbl()`](tda_mpbl.md),
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
A <- matrix(c(2, 1, 4, 5, 3, 8, 9, 7, 1), 3)
tda_mnrow(A); tda_mncol(A)
#> [1] 3
#> [1] 3
tda_mnorm(A); tda_mnorm1(A); tda_mnorm2(A); tda_mtrace(A)
#> [1] 9
#> [1] 40
#> [1] 15.81139
#> [1] 6
```
