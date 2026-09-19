# Print a tree's data structure

`ptree`: the same tree
[`plot.tda_graph`](https://janmarvin.github.io/TDA/reference/plot.tda_graph.md)
draws (via `pltree`), but as data rather than a picture – one row per
node, its parent, its first child, and its next sibling, the linked-list
representation TDA itself builds internally to draw one.

## Usage

``` r
tda_ptree(g, root = NULL, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- g:

  a
  [`tda_graph`](https://janmarvin.github.io/TDA/reference/tda_graph.md),
  undirected and shaped like a tree.

- root:

  which node to treat as the root; the first node in `g`'s data by
  default.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame: `node`, `edge_from`, `edge_to`, `value`, `parent`,
`first_child`, `next_sibling` (`0` where a node has none).

## Details

TDA's `gn=` (which graph, when more than one is open) is not optional
here in practice: `ptree`, unlike `pltree`, never calls the setup step
that gives it a safe default, so leaving it out reliably fails with
“current graph is not a tree” even for a graph that unambiguously is one
(`pltree` accepts the identical graph) – `gn = 1` is always sent
explicitly to work around it.

## See also

Other graph analysis:
[`plot.tda_graph()`](https://janmarvin.github.io/TDA/reference/plot.tda_graph.md),
[`tda_dmet()`](https://janmarvin.github.io/TDA/reference/tda_dmet.md),
[`tda_g()`](https://janmarvin.github.io/TDA/reference/tda_g.md),
[`tda_g_analyses`](https://janmarvin.github.io/TDA/reference/tda_g_analyses.md),
[`tda_graph()`](https://janmarvin.github.io/TDA/reference/tda_graph.md),
[`tda_locate_line()`](https://janmarvin.github.io/TDA/reference/tda_locate_line.md)

## Examples

``` r
e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 4), value = 1)
g <- tda_graph(e, directed = FALSE)   # 4 nodes, 3 edges: a tree
tda_ptree(g, root = 1)
#>   node edge_from edge_to value parent first_child next_sibling
#> 1    1         2       1     0      0           3            0
#> 2    2         3       1     1      1           4            0
#> 3    3         4       2     1      1           0            2
#> 4    4         0       0     1      2           0            0
```
