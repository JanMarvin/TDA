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
A <- matrix(c(1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1), 4, 4)
r <- tda_mpbl(A)
r$B
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    2    1    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    3    1
```
