# Simple Rasch model

TDA's `rmod` on binary item response data.

## Usage

``` r
tda_rmod(x, ...)
```

## Arguments

- x:

  matrix of 0/1 responses, persons in rows.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

centered item parameters in `$items` (joint ML is identified only up to
a constant, and its item estimates carry the known upward bias of about
k/(k-1) for k items), the pattern parameters in `$patterns`, and the
maximised log likelihood in `$logLik`.

## Examples

``` r
set.seed(1)
ab <- rnorm(80)
x <- sapply(c(-1, 1), function(d) rbinom(80, 1, plogis(ab - d)))
tda_rmod(x)$items   # easier item first, harder second
#>     item1     item2 
#>  2.590267 -2.590267 
```
