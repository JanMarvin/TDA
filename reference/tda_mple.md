# Product-limit (Kaplan-Meier) distribution and jumps

The product-limit (Kaplan-Meier) estimate of the CDF for a vector of
possibly right-censored event times – TDA's `mple`, via Efron's
redistribute-to-the-right algorithm, applied one observation at a time
in sorted order (events before censoring at the same time). When no two
events fall at exactly the same time, this matches the standard
product-limit CDF exactly – confirmed against
[`survival::survfit`](https://rdrr.io/pkg/survival/man/survfit.html)
across random instances (events may still tie with a censored
observation; only ties \*among events\* are excluded). When two or more
events tie, processing them one at a time rather than as a single
simultaneous risk-set reduction means only the *last* observation
processed within that tied group reaches the value standard Kaplan-Meier
software reports for that time – earlier ones in the same group show
smaller, order-dependent intermediate values. Confirmed by hand: two
events tied at `time = 2` out of 3 total observations give
`F = (1/3, 2/3, 1)` here, where standard Kaplan-Meier reports
`F(2) = 2/3` for both tied observations. This is a property of the
algorithm as TDA implements it, not a bug and not this wrapper's choice
to make. Separately: `F` always reaches exactly 1 at the highest-time
observation, *even when that observation is censored* – redistribution
only moves probability mass to later observations, never destroys it,
and the highest-time observation has nothing later to redistribute its
share to, so all of the mass ends up accounted for there by
construction. Standard Kaplan-Meier software instead leaves the survival
curve at whatever value it last reached, undefined beyond a censored
tail. Confirmed by hand: `time = c(1, 2, 3)`, `censored = c(0, 0, 1)`
gives `F = (1/3, 2/3, 1)` here. This is the low-level building block
behind TDA's
[`tda_ple`](https://janmarvin.github.io/TDA/reference/tda_ltb.md), which
most users want instead – this one works directly on plain vectors, with
no grouping, formula interface, or standard errors.

## Usage

``` r
tda_mple(time, censored, dir = tempfile("tda"))
```

## Arguments

- time:

  a numeric vector of event/censoring times.

- censored:

  a numeric vector the same length as `time`: 1 where the observation is
  censored, 0 where it is an observed event.

- dir:

  working directory.

## Value

A list, both the same length as `time`, in the original (unsorted)
order: `F`, the estimated CDF at each observation; `jumps`, the increase
in `F` at each observation (0 for a censored observation).

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
[`tda_mpit()`](https://janmarvin.github.io/TDA/reference/tda_mpit.md),
[`tda_mpz()`](https://janmarvin.github.io/TDA/reference/tda_mpz.md),
[`tda_mscal1()`](https://janmarvin.github.io/TDA/reference/tda_mscal1.md),
[`tda_msqrtd()`](https://janmarvin.github.io/TDA/reference/tda_msqrtd.md),
[`tda_msvd()`](https://janmarvin.github.io/TDA/reference/tda_msvd.md),
[`tda_mwvec()`](https://janmarvin.github.io/TDA/reference/tda_mwvec.md)

## Examples

``` r
tda_mple(c(2, 3, 3, 5, 7), c(0, 0, 1, 0, 0))
#> $F
#> [1] 0.2 0.4 0.4 0.7 1.0
#> 
#> $jumps
#> [1] 0.2 0.2 0.0 0.3 0.3
#> 
```
