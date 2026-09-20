# Iterate a Leslie matrix

Projects an age-structured population forward using a Leslie matrix –
TDA's `mpit`/`mpit1`. Row `i` of `fertility_survival` holds
`(fertility[i], survival[i])`: `fertility[i]` is the birth rate
contributed by age class `i` (used only for class 1 of the next
generation); `survival[i]` is the fraction of age class `i` that
survives into age class `i + 1` (`survival[n]`, the last row, is
unused). `tda_mpit1` adds a constant vector `immigration` to every age
class after each projection step.

## Usage

``` r
tda_mpit(fertility_survival, population, iterations, dir = tempfile("tda"))

tda_mpit1(
  fertility_survival,
  population,
  immigration,
  iterations,
  dir = tempfile("tda")
)
```

## Arguments

- fertility_survival:

  an n x 2 numeric matrix, columns `(fertility, survival)` as above.

- population:

  the initial population vector, length n.

- iterations:

  how many projection steps to take.

- dir:

  working directory.

- immigration:

  for `tda_mpit1`, a constant vector added to the population after each
  step, length n.

## Value

A matrix with `iterations + 1` rows (the initial population, then one
row per step) and `n` columns.

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
[`tda_mpbl()`](https://janmarvin.github.io/TDA/reference/tda_mpbl.md),
[`tda_mpfit()`](https://janmarvin.github.io/TDA/reference/tda_mpfit.md),
[`tda_mpinv()`](https://janmarvin.github.io/TDA/reference/tda_mpinv.md),
[`tda_mple()`](https://janmarvin.github.io/TDA/reference/tda_mple.md),
[`tda_mpz()`](https://janmarvin.github.io/TDA/reference/tda_mpz.md),
[`tda_mscal1()`](https://janmarvin.github.io/TDA/reference/tda_mscal1.md),
[`tda_msqrtd()`](https://janmarvin.github.io/TDA/reference/tda_msqrtd.md),
[`tda_msvd()`](https://janmarvin.github.io/TDA/reference/tda_msvd.md),
[`tda_mwvec()`](https://janmarvin.github.io/TDA/reference/tda_mwvec.md)

## Examples

``` r
F <- rbind(c(0, 0), c(2, 0.5), c(1, 0.3))
tda_mpit(F, c(100, 50, 20), iterations = 3)
#>      [,1] [,2] [,3]
#> [1,]  100   50   20
#> [2,]  120    0   25
#> [3,]   25    0    0
#> [4,]    0    0    0
```
