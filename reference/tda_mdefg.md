# Adjacency matrix of a graph

TDA's `mdefg`: builds the adjacency matrix of the current graph, with
`missing` substituted for absent edges (TDA's `sc=`, default -1 there; 0
here, which gives the usual weighted adjacency matrix).

## Usage

``` r
tda_mdefg(edges, missing = 0, directed = FALSE, valued = TRUE, ...)
```

## Arguments

- edges:

  edge list (data frame with from, to and optionally a value column).

- missing:

  value written where no edge exists.

- directed, valued:

  graph type flags, as in [`tda_g`](tda_g.md).

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the numeric adjacency matrix.

## Examples

``` r
e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3),
                value = c(5, 2, 7))
tda_mdefg(e)
#>      [,1] [,2] [,3]
#> [1,]    0    5    2
#> [2,]    5    0    7
#> [3,]    2    7    0
```
