# Check an ordering against row/column dominance

Flags each row/column index of a square matrix `A` where the current
index order fails a running dominance test against its mirror image –
TDA's `mch` (undocumented beyond "check row and column sum conditions"
in its own header; the rule below is the loop body's, not a named
textbook procedure). For each index `i`, scanning outward (`j` from
`i+1` to `n`, then separately from `i-1` down to `1`), `A[i, j]`
accumulates into a running total `z` and its mirror `A[j, i]`
accumulates into a running total `s`; index `i` is flagged the first
time `z` falls behind `s` by more than a small tolerance. In effect: is
`A[i, ]`'s share of an outward pair always at least as large,
cumulatively, as the corresponding mirrored entries in `A[, i]` – a
diagnostic for whether the current `1:n` ordering already respects that
dominance, without finding a better ordering the way
[`tda_mpz`](https://janmarvin.github.io/TDA/reference/tda_mpz.md)/
[`tda_mpbl`](https://janmarvin.github.io/TDA/reference/tda_mpbl.md) do.

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

Other matrix algebra:
[`tda_mcel()`](https://janmarvin.github.io/TDA/reference/tda_mcel.md),
[`tda_mcent()`](https://janmarvin.github.io/TDA/reference/tda_mcent.md),
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
# A[i, j] = 5 means i strongly precedes j.  Ordered 1 < 2 < 3 the
# rows dominate their columns and every index passes (0); the same
# matrix with the order reversed fails everywhere (1).
A <- matrix(c(0, 5, 5,
              1, 0, 5,
              1, 1, 0), 3, byrow = TRUE)
tda_mch(A)
#> [1] 0 0 0
tda_mch(A[3:1, 3:1])
#> [1] 1 1 1
```
