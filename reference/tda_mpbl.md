# Permute to block triangular form

Finds a symmetric permutation that puts a square matrix into block
lower-triangular (`tda_mpbl`) or block upper-triangular (`tda_mpbu`)
form (Duff & Reid, CACM 529) – TDA's `mpbl`/`mpbu`.
`B[i, j] == A[p[i], p[j]]` for every `i, j` (confirmed against raw TDA
output, not merely assumed from `m_perm`'s header comment); `block`
gives, for each row/column of `B`, which block it belongs to (blocks
numbered in the order `B` is arranged, so entries between an earlier and
a later block are always zero on the side `tda_mpbl`/`tda_mpbu` promises
– confirmed on an instance with cross-block coupling, not one where
every block trivially has size 1).

## Usage

``` r
tda_mpbl(A, dir = tempfile("tda"))

tda_mpbu(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a square numeric matrix.

- dir:

  working directory.

## Value

A list: `B`, the permuted, block-triangular matrix; `p`, the permutation
(`B == A[p, p]`); `block`, an integer vector the same length as `p`
giving each position's block number.

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
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
A <- matrix(c(1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1), 4, 4)
r <- tda_mpbl(A)
r$B
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    2    1    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    3    1
```
