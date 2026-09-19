# Quadratic assignment (GRASP)

TDA's `gqap`: given a flow graph and a distance graph (both undirected,
integer-valued), search for the assignment minimizing total cost, with
the GRASP heuristic of CACM algorithm 754.

## Usage

``` r
tda_gqap(flow, dist, max_iterations = 100, ...)
```

## Arguments

- flow, dist:

  symmetric integer matrices of equal dimension.

- max_iterations:

  iteration cap.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

with direct exports enabled, a list with `cost`, `assignment` and
`output`; otherwise the printed output lines.

## Examples

``` r
f <- rbind(c(0, 3, 1), c(3, 0, 2), c(1, 2, 0))
d <- rbind(c(0, 1, 4), c(1, 0, 2), c(4, 2, 0))
tda_gqap(f, d)
#> <tda_assignment>
#> cost: 22 
#> assignment: 1 2 3 
#> permuted: 9 x 6 table
#> (full TDA output in $output, 58 lines)
```
