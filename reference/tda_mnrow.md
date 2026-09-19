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
