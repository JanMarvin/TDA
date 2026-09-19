# Zero out matrix elements by a threshold test

Sets `A[i,j]` to 0 wherever the chosen comparison against `x` holds,
leaving every other element unchanged – TDA's `mnc` (documented in
`tda.hlp`). `test` chooses the comparison: `"<="`, `"<"`, `">="`, `">"`
(against the raw value), or the same four prefixed `"abs"` (against
`abs(A[i,j])`).

## Usage

``` r
tda_mnc(
  A,
  x,
  test = c("<=", "<", ">=", ">", "abs<=", "abs<", "abs>=", "abs>"),
  dir = tempfile("tda")
)
```

## Arguments

- A:

  a numeric matrix.

- x:

  the threshold.

- test:

  which comparison zeroes an element; see above.

- dir:

  working directory.

## Value

A numeric matrix.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_minvs()`](tda_minvs.md), [`tda_mkmet()`](tda_mkmet.md),
[`tda_mkp()`](tda_mkp.md), [`tda_mmul()`](tda_mmul.md),
[`tda_mnrow()`](tda_mnrow.md), [`tda_mpbl()`](tda_mpbl.md),
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
A <- matrix(c(1, 5, 3, -2, 4, 6), 2, 3)
tda_mnc(A, 3, "<=")   # zero out every element <= 3
#>      [,1] [,2] [,3]
#> [1,]    0    0    4
#> [2,]    5    0    6
```
