# Nonmetric multidimensional scaling

TDA's `mdsn1` (Kruskal) on a distance matrix.

## Usage

``` r
tda_mdsn(
  d,
  ties = c("ignore", "primary", "secondary"),
  max_iterations = 100,
  restarts = 5,
  ...
)
```

## Arguments

- d:

  symmetric distance matrix.

- ties:

  tie handling. Only `"ignore"` is implemented in `mdsn1` (its primary
  and secondary projections are commented out in the C, and the option
  would silently do nothing); the other two values error and point to
  [`tda_mds`](tda_mds.md), whose `mdsn` has them.

- max_iterations:

  iteration cap.

- restarts:

  random starting configurations, TDA's `ns=`; the one with the lowest
  stress is returned. Nonmetric MDS on few points has degenerate local
  minima (clusters of near-coincident points with small stress), and
  restarts are the guard against them.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

with direct exports enabled, a list with `stress`, `configuration` (n x
2) and `output`; otherwise the printed output lines.

## Examples

``` r
d <- as.matrix(dist(cbind(c(0, 0, 3, 3), c(0, 2, 0, 2))))
tda_mdsn(d)
#> <tda_mds_nonmetric>
#> stress: 0 
#> configuration: 4 x 2 matrix
#> (full TDA output in $output, 408 lines)
```
