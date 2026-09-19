# Check an ordering against row/column dominance

Flags each row/column index of a square matrix `A` where the current
index order fails a running dominance test against its mirror image –
TDA's `mch` (undocumented beyond "check row and column sum conditions"
in its own header; the exact rule below was read from the loop body and
confirmed against a from-scratch re-implementation of the same loops,
not derived from a named textbook procedure). For each index `i`,
scanning outward (`j` from `i+1` to `n`, then separately from `i-1` down
to `1`), `A[i, j]` accumulates into a running total `z` and its mirror
`A[j, i]` accumulates into a running total `s`; index `i` is flagged the
first time `z` falls behind `s` by more than a small tolerance. In
effect: is `A[i, ]`'s share of an outward pair always at least as large,
cumulatively, as the corresponding mirrored entries in `A[, i]` – a
diagnostic for whether the current `1:n` ordering already respects that
dominance, without finding a better ordering the way
[`tda_mpz`](tda_mpz.md)/ [`tda_mpbl`](tda_mpbl.md) do.

## Usage

``` r
tda_mch(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a square numeric matrix.

- dir:

  working directory.

## Value

An integer vector, length `nrow(A)`: 1 where that index fails the test,
0 where it passes.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mcross()`](tda_mcross.md),
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
tda_mch(matrix(c(0, 1, 0, 5, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE))
#> [1] 1 1 0
```
