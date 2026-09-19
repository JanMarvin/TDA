# A graph for TDA's graph commands

A graph for TDA's graph commands

## Usage

``` r
tda_graph(
  edges,
  directed = TRUE,
  form = c("edges", "dissimilarity"),
  options = list()
)
```

## Arguments

- edges:

  an edge list – from, to, and one value, or two values when the graph
  carries a different one in each direction, as the shipped `gd1.dat`
  does – or a square adjacency matrix, or a `dist`.

- directed:

  whether the graph is directed. `FALSE` emits `gt=2`. Several commands
  –
  [`tda_g_components`](https://janmarvin.github.io/TDA/reference/tda_g_analyses.md)
  among them, and
  [`plot.tda_graph`](https://janmarvin.github.io/TDA/reference/plot.tda_graph.md)'s
  `layout = "tree"` – only work on an undirected graph, and TDA itself
  refuses a directed one for them; `directed = TRUE` (the default here)
  is right for ordinary, arrow-drawn graphs, but a tree needs `FALSE`,
  since a tree's layout comes from which nodes are connected, not which
  way an edge happens to point.

- form:

  how TDA is to read it: `"edges"` for an edge list, which is `gdd`
  option 1, or `"dissimilarity"` for the upper triangle of a
  dissimilarity matrix, option 4.

- options:

  a named list of further `gdd` options.

## Value

A `tda_graph`.

## See also

Other graph analysis:
[`plot.tda_graph()`](https://janmarvin.github.io/TDA/reference/plot.tda_graph.md),
[`tda_dmet()`](https://janmarvin.github.io/TDA/reference/tda_dmet.md),
[`tda_g()`](https://janmarvin.github.io/TDA/reference/tda_g.md),
[`tda_g_analyses`](https://janmarvin.github.io/TDA/reference/tda_g_analyses.md),
[`tda_locate_line()`](https://janmarvin.github.io/TDA/reference/tda_locate_line.md),
[`tda_ptree()`](https://janmarvin.github.io/TDA/reference/tda_ptree.md)

## Examples

``` r
e <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4), value = 1)
# gcon works on an undirected graph only, so `directed = FALSE`
tda_g_components(tda_graph(e, directed = FALSE))
#> Call: tda_g_components(tda_graph(e, directed = FALSE))
#> 
#> Cases: 4 
#> 
#>  component n i node_i
#>          1 4 1      1
#>          1 4 2      2
#>          1 4 3      3
#>          1 4 4      4
```
