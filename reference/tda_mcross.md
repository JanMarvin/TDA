# Cross-product and an arithmetic sequence as a column

`tda_mcross` is TDA's `mcross`, matching base R's `crossprod(A)`
(`t(A) %*% A`). `tda_mnum` builds the length-`n` arithmetic sequence
starting at `x` with step `d` – TDA's `mnum`, matching
`seq(x, by = d, length.out = n)`.

## Usage

``` r
tda_mcross(A, dir = tempfile("tda"))

tda_mnum(x, d, n, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix.

- dir:

  working directory.

- x:

  the starting value, for `tda_mnum`.

- d:

  the step, for `tda_mnum`.

- n:

  the length, for `tda_mnum`.

## Value

`tda_mcross` returns a matrix; `tda_mnum` returns a plain numeric
vector.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mdiag()`](tda_mdiag.md), [`tda_mev()`](tda_mev.md),
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
tda_mcross(matrix(1:6, 3, 2))    # == crossprod(matrix(1:6, 3, 2))
#>      [,1] [,2]
#> [1,]   14   32
#> [2,]   32   77
tda_mnum(5, 2, 4)                # == seq(5, by = 2, length.out = 4)
#> [1]  5  7  9 11
```
