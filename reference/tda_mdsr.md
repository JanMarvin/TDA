# Axes in an MDS configuration

TDA's `mdsr`: regresses external variables onto a two-dimensional
configuration to draw interpretable axes.

## Usage

``` r
tda_mdsr(config, x, ...)
```

## Arguments

- config:

  two-column matrix, the configuration.

- x:

  matrix of external variables (same row count).

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the printed output.

## Examples

``` r
cfg <- cbind(c(0, 1, 2, 3), c(0, 1, 0, 1))
r <- tda_mdsr(cfg, x = cbind(v = c(1, 2, 3, 4)))
cat(head(tda_payload(r), 12), sep = "\n")
#> mdsr(...)=...
#> Find axes in configuration. Current memory: 402084 bytes.
#> Configuration: X, Y
#> X minimum:     0.0000   maximum:     3.0000   mean:     1.5000 
#> Y minimum:     0.0000   maximum:     1.0000   mean:     0.7500 
#> Variable: V1
#> Optimal phi: 0
#> Rank correlation: 1
#> Axis: (0,0.75) to (3,0.75)
```
