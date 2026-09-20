# Kemeny distance between rankings

Pairwise Kemeny (Kendall-tau-with-ties) distance between every pair of
rows of `A`, each row a ranking (or score) over the same set of items –
TDA's `mkmet`. For two rankings, each pair of items contributes 0 if
both rankings order (or tie) the pair the same way, 1 if one ranking
ties the pair and the other doesn't, or 2 if the rankings strictly
disagree on the pair's order.

## Usage

``` r
tda_mkmet(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix, one ranking per row.

- dir:

  working directory.

## Value

A symmetric numeric matrix with zero diagonal, the Kemeny distance
between every pair of rows.

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
# four rankings of the same four items: b swaps the last two of a,
# c reverses a, d ties items in pairs
R <- rbind(a = c(1, 2, 3, 4),
           b = c(1, 2, 4, 3),
           c = c(4, 3, 2, 1),
           d = c(1, 1, 2, 2))
tda_mkmet(R)   # a-b: 2 (one pair reversed); a-c: 12 (all six pairs)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    2   12    2
#> [2,]    2    0   10    2
#> [3,]   12   10    0   10
#> [4,]    2    2   10    0
```
