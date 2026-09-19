# Draw a graph with explicit node positions

`plg` is TDA's graph drawing command, and a different thing from
`pltree`: it takes the node coordinates and the edges and draws them,
with control over node shape and shading, edge curvature and arrowheads.
`pltree` by contrast computes a tree layout and draws edges only.

## Usage

``` r
tda_pl_graph(p, nodes, edges = NULL, ...)
```

## Arguments

- p:

  a [`tda_ps`](tda_ps.md) session.

- nodes:

  a data frame with the node positions. It needs `x` and `y`; `id`
  numbers the nodes, defaulting to the row number, `grey` and `shape`
  set `gs=` and `gt=`, `label` a text string instead of the node number
  (`str=`), `size` the node's display size in mm (`rd=`), and `lty`,
  `lw`, `fontsize` its line type, line width, and font size
  (`lt=`/`lw=`/`fs=`).

- edges:

  a data frame of edges with `from` and `to`, holding node ids. `curve`
  sets `rd=`, the radius that bends an edge, `arrow` a pair of arrowhead
  dimensions, `label` the value written on the edge (`ic=`), `lty` its
  line type, which is how the two directions of a reciprocal pair are
  told apart, `lw` its line width, `fontsize` the edge label's font size
  (`fs=`), and `loop`, for a self-edge, which side of the node it loops
  from (`dir=`).

- ...:

  further options for `plg`, e.g. `nmax=` (the maximum number of nodes,
  rarely needed since TDA's default of 100 is usually enough).

## Value

The session, with the command added.

## See also

Other plotting: [`plot.tda_ple()`](plot.tda_ple.md),
[`tda_check_ps()`](tda_check_ps.md), [`tda_pl()`](tda_pl.md),
[`tda_pl_arc()`](tda_pl_arc.md), [`tda_pl_axis()`](tda_pl_axis.md),
[`tda_pl_hist`](tda_pl_hist.md), [`tda_pl_panel()`](tda_pl_panel.md),
[`tda_pl_regression()`](tda_pl_regression.md),
[`tda_pl_scatter()`](tda_pl_scatter.md),
[`tda_plot_ps()`](tda_plot_ps.md), [`tda_ps()`](tda_ps.md),
[`tda_ps3()`](tda_ps3.md), [`tda_read_ps()`](tda_read_ps.md)

## Examples

``` r
nodes <- data.frame(x = c(1, 3, 1, 3), y = c(2, 2, 1, 1))
edges <- data.frame(from = c(1, 1), to = c(2, 3))
p <- tda_ps(xlim = c(0, 4), ylim = c(0, 3), width = 40, height = 30)
p <- tda_pl_graph(p, nodes, edges)
```
