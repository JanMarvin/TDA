# Balanced repeated replication design

`brr`: builds the replicate-weight matrix balanced repeated replication
needs – a variance-estimation technique from survey statistics that
works by refitting an estimator on `nr` half-sample replicates of the
data (`nr`, the number of replications required, is determined from
`strata` and `secu` and reported back) and taking the spread across them
as the variance, rather than a model-based formula. This returns the
design itself – TDA's orthogonal (Hadamard-based) construction, not
practical to reproduce independently in R – so it can be applied to
whatever estimator is being replicated; actually refitting a model `nr`
times and computing the BRR variance from the results is ordinary R code
once the design is in hand, not shown here since it depends entirely on
what is being estimated.

## Usage

``` r
tda_brr(strata, secu, dir = tempfile("tda"), ...)
```

## Arguments

- strata:

  number of strata.

- secu:

  number of sampling units per stratum; `2` in the overwhelming majority
  of real BRR designs, and the only value tested here.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A list: `replications` (`nr`, how many are needed), and `design`, the
`strata` by `replications` matrix itself (TDA's description:
`design[i, j]` is which sampling unit of stratum `i` to use in
replication `j`).

## Details

Not every `(strata, secu)` combination succeeds: `secu` must be at least
2 (one stratum, one unit gives nothing to replicate), and at least one
combination that satisfies every other documented constraint
(`brr = 4,2`) can still fail with “no success in calculating orthogonal
weights”: 2, 3, and 8 strata at `secu = 2` all succeed, 4 does not. This
is a limitation in TDA's Galois-field construction for specific sizes,
not something this wrapper can lift.

## See also

Other descriptive statistics:
[`tda_cov()`](https://janmarvin.github.io/TDA/reference/tda_cov.md),
[`tda_dstat()`](https://janmarvin.github.io/TDA/reference/tda_dstat.md),
[`tda_freq()`](https://janmarvin.github.io/TDA/reference/tda_freq1.md),
[`tda_independence()`](https://janmarvin.github.io/TDA/reference/tda_independence.md),
[`tda_ineq()`](https://janmarvin.github.io/TDA/reference/tda_ineq.md),
[`tda_loglin()`](https://janmarvin.github.io/TDA/reference/tda_loglin.md),
[`tda_quant()`](https://janmarvin.github.io/TDA/reference/tda_quant.md),
[`tda_rcorr()`](https://janmarvin.github.io/TDA/reference/tda_rcorr.md),
[`tda_segr()`](https://janmarvin.github.io/TDA/reference/tda_segr.md),
[`tda_subm()`](https://janmarvin.github.io/TDA/reference/tda_subm.md)

## Examples

``` r
b <- tda_brr(3, 2)
b$replications
#> [1] 4
b$design
#>          rep1 rep2 rep3 rep4
#> stratum1    0    0    1    0
#> stratum2    0    1    0    0
#> stratum3    1    0    0    0
```
