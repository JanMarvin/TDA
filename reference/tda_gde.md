# Undirected graph from node lists

TDA's `gde`: builds an undirected graph from a data matrix where `node`
numbers each unit and `id` identifies the objects they share; units
sharing an id become adjacent.

## Usage

``` r
tda_gde(node, id, output = c("edges", "matrix"), loops = TRUE, ...)
```

## Arguments

- node, id:

  integer vectors of equal length.

- output:

  "edges" for an edge list, "matrix" for the square adjacency matrix.

- loops:

  FALSE removes loops.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the values TDA wrote, as a matrix (one row per line of the output file),
from the export channel.

## Examples

``` r
tda_gde(node = c(1, 2, 2, 3), id = c(10, 10, 20, 20))
#>      [,1] [,2] [,3]
#> [1,]    1    2    1
#> [2,]    2    3    1
```
