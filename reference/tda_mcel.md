# Build an edge list from an adjacency-style matrix

Lists every `(i, j)` with `A[i, j] >= x` as one row `(i, j, A[i, j])` –
TDA's `mcel`, scanning row by row. Refuses if no element of `A` reaches
`x` (there would be nothing to return).

## Usage

``` r
tda_mcel(A, x, dir = tempfile("tda"))
```

## Arguments

- A:

  a square numeric matrix.

- x:

  the threshold; an entry is an edge when it is at least this large.

- dir:

  working directory.

## Value

A 3-column numeric matrix: `i`, `j`, and the value `A[i, j]`, one row
per edge.

## See also

Other matrix algebra: [`tda_mcent()`](tda_mcent.md),
[`tda_mch()`](tda_mch.md), [`tda_mcross()`](tda_mcross.md),
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
tda_mcel(matrix(c(5, 0, 1, 0, 3, 0, 2, 0, 4), 3, 3), 3)
#>      [,1] [,2] [,3]
#> [1,]    1    1    5
#> [2,]    2    2    3
#> [3,]    3    3    4
```
