# Draw a graph

TDA's graph-drawing command, `pltree`, only draws trees. For anything
else – a cycle, more than one component, or a directed graph – this lays
the nodes out on a circle and draws it with
[`tda_pl_graph`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md)
instead, which takes explicit coordinates and so can draw any graph.

## Usage

``` r
# S3 method for class 'tda_graph'
plot(
  x,
  width = 110,
  height = 80,
  file = "graph.ps",
  newpage = TRUE,
  layout = c("auto", "circle", "tree"),
  ...
)
```

## Arguments

- x:

  a
  [`tda_graph`](https://janmarvin.github.io/TDA/reference/tda_graph.md).

- width, height:

  size of the plotting area in millimetres.

- file:

  name for the PostScript file.

- newpage:

  start a new page before drawing.

- layout:

  `"auto"` (the default) uses `pltree` when the graph is actually a tree
  (connected, one fewer edge than nodes) and a circular layout
  otherwise; `"circle"` always uses the circular layout, `"tree"` always
  tries `pltree` and lets TDA refuse a non-tree graph with its own
  error.

- ...:

  further options for `pltree` (`layout = "tree"`) or
  [`tda_pl_graph`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md)
  (`layout = "circle"`).

## Value

The parsed drawing operations, invisibly.

## See also

Other graph analysis:
[`tda_dmet()`](https://janmarvin.github.io/TDA/reference/tda_dmet.md),
[`tda_g()`](https://janmarvin.github.io/TDA/reference/tda_g.md),
[`tda_g_analyses`](https://janmarvin.github.io/TDA/reference/tda_g_analyses.md),
[`tda_graph()`](https://janmarvin.github.io/TDA/reference/tda_graph.md),
[`tda_locate_line()`](https://janmarvin.github.io/TDA/reference/tda_locate_line.md),
[`tda_ptree()`](https://janmarvin.github.io/TDA/reference/tda_ptree.md)

## Examples

``` r
e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 4), value = 1)
g <- tda_graph(e, directed = FALSE)  # 4 nodes, 3 edges: already a tree
pdf(NULL); plot(g); dev.off()
#> agg_record_1f865f44eb69 
#>                       2 

# a graph with a cycle is not a tree, and this one is directed besides --
# both are fine with the default circular layout
e2 <- data.frame(from = c(1, 1, 2, 2, 3, 4), to = c(2, 3, 3, 4, 4, 5))
g2 <- tda_graph(e2, directed = TRUE)
pdf(NULL); plot(g2); dev.off()
#> agg_record_1f865f44eb69 
#>                       2 
```
