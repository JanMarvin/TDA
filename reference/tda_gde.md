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

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the values TDA wrote, as a matrix (one row per line of the output file),
from the export channel.

## Examples

``` r
# three persons (node) and the organisations they belong to (id):
# person 1 in 10 and 20, person 2 in 20 and 30, person 3 in 30.
# Sharing an organisation makes two persons adjacent: 1-2 and 2-3.
tda_gde(node = c(1, 1, 2, 2, 3), id = c(10, 20, 20, 30, 30))
#>      [,1] [,2] [,3]
#> [1,]    1    2    1
#> [2,]    2    3    1
```
