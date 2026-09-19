# Forward control flows in a directed graph

TDA's `gfcf` on a directed graph loaded from an edge list.

## Usage

``` r
tda_gfcf(edges, threshold = 0.5, variant = 1, ...)
```

## Arguments

- edges:

  data frame: from, to.

- threshold:

  control threshold, default 0.5.

- variant:

  1 or 2 (the two forms the command documents).

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

one row per node: the node, how many nodes it controls, and the
controlled nodes as a comma-joined string; the pair table TDA writes
(node, controlled node, path length L, share S, M, R – the manual's Box
2 of 7.6.1.1) is the `"pairs"` attribute.

## Examples

``` r
tda_gfcf(data.frame(i = c(1, 2, 3), j = c(2, 3, 4)))
#>   node n_controlled controlled
#> 1    1            3      2,3,4
#> 2    2            2        3,4
#> 3    3            1          4
#> 4    4            0           
```
