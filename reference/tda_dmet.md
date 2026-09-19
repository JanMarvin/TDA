# Convert a distance matrix into a metric

`dmet`: an undirected valued graph's edge values, read as a distance
matrix, do not always satisfy the triangle inequality (a direct edge can
be "longer" than going via a third node). This finds the smallest
constant that, added to every distance, makes the whole matrix satisfy
it – \\c = \max(0, \max\_{i,j,k} d\_{ij} - d\_{ik} - d\_{kj})\\ – and
returns both the original and the corrected matrix.

## Usage

``` r
tda_dmet(g, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- g:

  a [`tda_graph`](tda_graph.md), undirected and valued.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A list: `original` and `modified`, both distance matrices, and
`constant`, the value added to get from one to the other.

## See also

Other graph analysis: [`plot.tda_graph()`](plot.tda_graph.md),
[`tda_g()`](tda_g.md), [`tda_g_analyses`](tda_g_analyses.md),
[`tda_graph()`](tda_graph.md),
[`tda_locate_line()`](tda_locate_line.md), [`tda_ptree()`](tda_ptree.md)

## Examples

``` r
e <- data.frame(from = c(1, 2, 3), to = c(2, 3, 1), value = c(1, 2, 5))
g <- tda_graph(e, directed = FALSE)   # 1-3 direct (5) exceeds
                                      # 1-2-3 via node 2 (1+2=3)
tda_dmet(g)
#> $original
#>      [,1] [,2] [,3]
#> [1,]    0    1    5
#> [2,]    1    0    2
#> [3,]    5    2    0
#> 
#> $modified
#>      [,1] [,2] [,3]
#> [1,]    0    3    7
#> [2,]    3    0    4
#> [3,]    7    4    0
#> 
#> $constant
#> [1] 2
#> 
```
