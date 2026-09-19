# Graph analyses

One function per graph command, each a thin wrapper: the command name,
the options it takes and the table TDA writes. Everything TDA offers
that has no named wrapper is still reachable with
[`tda_g`](https://janmarvin.github.io/TDA/reference/tda_g.md).

## Usage

``` r
tda_g_components(g, edges = FALSE, ...)

tda_g_shortest(g, format = c("reachable", "matrix", "matrix_no_labels"), ...)

tda_g_mst(g, algorithm = c("kruskal", "prim"), max = FALSE, sort = FALSE, ...)

tda_g_cliques(
  g,
  algorithm = c("bron_kerbosch", "harary_ross"),
  min_size = NULL,
  sort = FALSE,
  ...
)

tda_g_cycles(
  g,
  format = c("fundamental_nodes", "fundamental_edges", "all_v1", "all_v2"),
  ...
)

tda_g_cutpoints(g, per_node = FALSE, ...)

tda_g_flow(g, matrix = FALSE, ...)

tda_g_degrees(g)

tda_g_edges(
  g,
  format = c("edges", "lower_triangle", "lower_triangle_diag", "full", "full_square"),
  ...
)

tda_g_transitive(g, ...)

tda_g_aggregate(g, ...)

tda_g_symmetric(g, ...)

tda_g_toposort(g, relabel_graph = FALSE, ...)

tda_g_nodes(g, ...)

tda_g_dot(g, ...)

tda_g_links(g, backward = FALSE, ...)

tda_g_dcycles(
  g,
  format = c("per_cycle", "per_node", "progress"),
  stop_at_first = FALSE,
  ...
)

tda_g_dblocks(g, ...)

tda_g_compact(g, per_node = FALSE, ...)

tda_g_gcliques(g, ...)

tda_g_random(g, ...)

tda_g_ownership(g, ...)

tda_g_backward(g, ...)

tda_g_paths(g, ...)

tda_g_flowcontrol(g, ...)

tda_g_centred(g, ...)

tda_g_spantrees(g, edges = FALSE, ...)

tda_g_reachable(
  g,
  format = c("count", "list", "strong_nodes", "strong_edges"),
  ...
)

tda_g_neighbourhoods(g, ...)

tda_g_independent(g, min_size = NULL, use_complement = FALSE, ...)

tda_g_spanning(g)

tda_g_eigen(g, n = NULL, tol = NULL, control = NULL, ...)
```

## Arguments

- g:

  a
  [`tda_graph`](https://janmarvin.github.io/TDA/reference/tda_graph.md).

- edges:

  for `tda_g_components`, report each component as an edge list instead
  of a node list (`opt=2`). For `tda_g_spantrees`, an edge list instead
  of one record per tree (`gnst`'s `opt=2`).

- ...:

  options for the command.

- format:

  for `tda_g_shortest`: `"reachable"` (default, only reachable pairs),
  `"matrix"` (a square matrix with row/ column labels), or
  `"matrix_no_labels"` – `gsp`'s `opt=`. For `tda_g_cycles`:
  `"fundamental_nodes"` (default), `"fundamental_edges"`, `"all_v1"`, or
  `"all_v2"` – `gcyc`'s `opt=`. For `tda_g_edges`: `"edges"` (default),
  or the adjacency matrix as `"lower_triangle"`,
  `"lower_triangle_diag"`, `"full"`, or `"full_square"` – `gdp`'s
  `opt=`. For `tda_g_dcycles`: `"per_cycle"` (default), `"per_node"`, or
  `"progress"` (report progress to the console instead of writing
  anything) – `gdcyc`'s `opt=`. For `tda_g_reachable`: `"count"`
  (default, number of reachable nodes), `"list"` (plus the node
  numbers), `"strong_nodes"`, or `"strong_edges"` (strongly connected
  components) – `gdcon`'s `opt=`. Five different commands' own `opt=`,
  each with its own meaning.

- algorithm:

  for `tda_g_mst`: `"kruskal"` (default) or `"prim"`. For
  `tda_g_cliques`: `"bron_kerbosch"` (default) or `"harary_ross"`.
  Different commands' own `alg=`.

- max:

  for `tda_g_mst`, find a maximum spanning tree instead of a minimum one
  (`max=1`).

- sort:

  for `tda_g_mst` and `tda_g_cliques`, sort the output by node number.

- min_size:

  for `tda_g_cliques`, the minimum clique size (`gcliq`'s `min=`,
  default 3). For `tda_g_independent`, the minimum size of the
  independent set to search for (`giset`'s `min=`).

- per_node:

  for `tda_g_cutpoints`, one record per node instead of one per
  component (`gcut`'s `opt=2`). For `tda_g_compact`, the same idea for
  `gcset`'s `opt=2`.

- matrix:

  for `tda_g_flow`, report the result as a square flow matrix instead of
  one record per flow (`opt=2`).

- relabel_graph:

  for `tda_g_toposort`, return the graph with its new labels instead of
  the list of old/new labels (`gsort`'s own `opt=2`).

- backward:

  for `tda_g_links`, backward links instead of forward ones (`gdln`'s
  `opt=2`).

- stop_at_first:

  for `tda_g_dcycles`, stop after the first cycle is found (`gdcyc`'s
  `dopt=1`).

- use_complement:

  for `tda_g_independent`, search the complementary graph instead
  (`giset`'s `cg=1`).

- n:

  for `tda_g_eigen`, the number of eigenvalues/vectors to compute
  (`gev`'s `n=`, default 1).

- tol:

  for `tda_g_eigen`, the required accuracy (`gev`'s `eps=`).

- control:

  for `tda_g_eigen`, convergence settings from
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md)
  (only its `maxit` field applies, as `gev`'s `mxit=`).

## Value

An object carrying `table` – the records the command writes, in that
command's layout (edges for the tree and path commands, node rows for
the centrality and component ones; see the command's
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md)
entry) – or, for the few commands whose output is ragged or non-tabular,
`table` as a ragged record list or `text` as plain lines.

`tda_g_backward` also carries `summary`, one row per node, from the
second file `gbcf` writes; TDA does not produce it unless a file is
named for it, so the wrapper asks by default.

## Details

The descriptions below are Rohwer's, from TDA's manual; the exact
wording and the full option list for any command are available with
`tda_help("gmst")` and so on.

## Connectivity and structure

- `tda_g_components`:

  connected components (`gcon`).

- `tda_g_reachable`:

  for each node of a directed graph, the set of nodes reachable from it
  by a directed path (`gdcon`).

- `tda_g_cutpoints`:

  cut nodes of an undirected graph (`gcut`).

- `tda_g_transitive`:

  transitive closure (`gtcl`).

- `tda_g_symmetric`:

  whether every edge from i to j is matched by one from j to i (`gsym`).

- `tda_g_toposort`:

  topological sort of a directed graph (`gsort`).

## Subgraphs

- `tda_g_cliques`:

  all cliques – maximally connected subgraphs – of an undirected graph
  (`gcliq`).

- `tda_g_independent`:

  a maximal independent set, found by a greedy randomized search
  procedure (GRASP), so the answer is a good one rather than a proven
  optimum (`giset`).

- `tda_g_compact`:

  compact sets in an undirected valued graph (`gcset`).

- `tda_g_neighbourhoods`:

  compares the neighbourhoods of the nodes and groups those that agree.
  Needs an unvalued graph; the wrapper arranges that (`gcni`).

## Paths, trees and cycles

- `tda_g_shortest`:

  a shortest path for each pair of nodes. Directed or undirected, valued
  or unvalued; in an unvalued graph every edge counts as 1 (`gsp`).

- `tda_g_paths`:

  all paths between any two nodes, with their minimal or maximal length
  (`gep`).

- `tda_g_mst`:

  a minimum, or maximum, spanning tree for each connected component of
  an undirected valued graph (`gmst`).

- `tda_g_spantrees`:

  enumerates every spanning tree of an undirected graph (`gnst`).

- `tda_g_cycles`:

  a fundamental set of cycles, and optionally all cycles, of an
  undirected graph (`gcyc`).

- `tda_g_dcycles`:

  cycles in a directed graph (`gdcyc`).

- `tda_g_links`:

  forward and backward links in a directed graph (`gdln`).

## Flow, control and ownership

- `tda_g_flow`:

  a maximal flow between each pair of nodes joined by a directed path
  (`gflow`).

- `tda_g_flowcontrol`:

  the quantities from which measures of flow control are built, for a
  directed valued graph (`gfc`).

- `tda_g_ownership`:

  integrated ownership in a directed valued graph (`gio`).

- `tda_g_backward`:

  direct and indirect backward control (`gbcf`).

- `tda_g_centred`:

  node-centred networks (`gnc`).

## Description and output

- `tda_g_nodes`:

  information about the nodes of the current graph (`gni`).

- `tda_g_eigen`:

  eigenvalues and eigenvectors of the adjacency matrix of an undirected
  graph (`gev`).

- `tda_g_aggregate`:

  aggregates the nodes of a valued graph (`gda`).

- `tda_g_dblocks`:

  blocks of a directed graph (`gdcset`).

- `tda_g_random`:

  generates an unvalued graph with a given number of nodes and edges
  (`gcd`).

- `tda_g_dot`, `tda_g_edges`:

  write the graph out as a Graphviz dot file, or an edge list/adjacency
  matrix, respectively (`gdot`, `gdp`).

- `tda_g_degrees`:

  node degrees from `gni`, TDA's node-information command: in-degree,
  out-degree, loops, and a degree column (for an undirected graph each
  incident edge is reported on both sides, so the degree is the
  in-degree, not the sum). Takes no further arguments.

## What each command requires

TDA refuses a command whose graph is of the wrong kind rather than
answering wrongly, so the error message is worth reading. Most of the
`gd*` commands want a directed graph, `gcliq`, `gcut`, `gcyc`, `gmst`,
`gnst`, `gcset` and `gnc` want an undirected one, and `gcni` wants an
unvalued one.

## See also

Other graph analysis:
[`plot.tda_graph()`](https://janmarvin.github.io/TDA/reference/plot.tda_graph.md),
[`tda_dmet()`](https://janmarvin.github.io/TDA/reference/tda_dmet.md),
[`tda_g()`](https://janmarvin.github.io/TDA/reference/tda_g.md),
[`tda_graph()`](https://janmarvin.github.io/TDA/reference/tda_graph.md),
[`tda_locate_line()`](https://janmarvin.github.io/TDA/reference/tda_locate_line.md),
[`tda_ptree()`](https://janmarvin.github.io/TDA/reference/tda_ptree.md),
[`tda_ucl()`](https://janmarvin.github.io/TDA/reference/tda_ucl.md)

## Examples

``` r
# a small directed graph: 1 and 2 both feed 3, which feeds 4, which
# feeds 5
e <- data.frame(from = c(1, 1, 2, 2, 3, 4),
                to   = c(2, 3, 3, 4, 4, 5),
                value = c(1, 2, 1, 3, 2, 1))

# undirected questions want an undirected graph
g <- tda_graph(e, directed = FALSE)
tda_g_components(g)$table
#>   component n i node_i
#> 1         1 5 1      1
#> 2         1 5 2      2
#> 3         1 5 3      3
#> 4         1 5 4      4
#> 5         1 5 5      5
tda_g_cliques(g)$table
#>   V1 V2 V3 V4 V5 V6
#> 1  1  1  3  2  3  1
#> 2  1  2  3  2  3  4
tda_g_mst(g)$table
#>   tree node_i node_j value
#> 1    1      1      2     1
#> 2    1      2      3     1
#> 3    1      4      5     1
#> 4    1      3      4     2
tda_g_cycles(g)$table
#> 2 records, lengths 5-6
#> [1]  1 1 3 1 2 
#> [2]  1 2 4 3 1 2 
tda_g_cutpoints(g)$table
#>   V1 V2 V3 V4
#> 1  1  5  1  4
tda_g_compact(g)$table
#>   V1 V2 V3 V4 V5 V6 V7
#> 1  1  1  2  1  2  4  5
tda_g_spantrees(g)$table
#>   V1 V2 V3 V4 V5 V6 V7
#> 1  1  1  0  1  1  2  4
#> 2  1  2  0  1  1  3  4
#> 3  1  3  0  1  2  2  4
#> 4  1  4  0  1  2  3  4
#> 5  1  5  0  1  4  2  4
#> 6  1  6  0  3  1  2  4
#> 7  1  7  0  3  1  3  4
#> 8  1  8  0  4  1  3  4
tda_g_spanning(g)$table
#>   tree node_i node_j value
#> 1    1      1      2     1
#> 2    1      2      3     1
#> 3    1      3      4     2
#> 4    1      4      5     1
tda_g_independent(g)$table
#>   V1 V2 V3
#> 1  1  2  2
#> 2  2  5  5
tda_g_centred(g)$table
#>   V1 V2 V3 V4        V5
#> 1  1  1  3  3 1.0000000
#> 2  2  2  4  5 0.8333333
#> 3  3  3  4  5 0.8333333
#> 4  4  4  4  4 0.6666667
#> 5  5  5  2  1 1.0000000

# directed ones want a directed graph
d <- tda_graph(e, directed = TRUE)
tda_g_toposort(d)$table
#>   i node_i label_i
#> 1 1      1       1
#> 2 2      2       2
#> 3 3      3       3
#> 4 4      4       4
#> 5 5      5       5
tda_g_reachable(d)$table
#>   V1 V2 V3
#> 1  1  1  4
#> 2  2  2  3
#> 3  3  3  2
#> 4  4  4  1
tda_g_transitive(d)$table
#>   V1 V2 V3 V4 V5 V6 V7
#> 1  1  1  0  1  2  4  5
#> 2  2  2 -1  0  1  3  4
#> 3  3  3 -1 -1  0  2  3
#> 4  4  4 -1 -1 -1  0  1
#> 5  5  5 -1 -1 -1 -1  0
tda_g_nodes(d)$table
#>   V1 V2 V3 V4 V5
#> 1  1  1  0  2  0
#> 2  2  2  1  2  0
#> 3  3  3  2  1  0
#> 4  4  4  2  1  0
#> 5  5  5  1  0  0
tda_g_links(d)$table
#>   index node n_links links
#> 1     1    1       2  2, 3
#> 2     2    2       2  3, 4
#> 3     3    3       1     4
#> 4     4    4       1     5
tda_g_paths(d)$table
#>    V1 V2 V3 V4 V5 V6 V7 V8 V9 V10
#> 1   1  2  1  2  1  1  1  1  1   1
#> 2   1  3  1  3  2  2  1  2  2   2
#> 3   1  4  1  4  3  3  2  4  3   4
#> 4   1  5  1  5  3  3  3  5  4   5
#> 5   2  3  2  3  1  1  1  1  1   1
#> 6   2  4  2  4  2  2  1  3  2   3
#> 7   2  5  2  5  2  2  2  4  3   4
#> 8   3  4  3  4  1  1  1  2  1   2
#> 9   3  5  3  5  1  1  2  3  2   3
#> 10  4  5  4  5  1  1  1  1  1   1
tda_g_flow(d)$table
#>    index i j node_i node_j flow
#> 1      1 1 2      1      2    1
#> 2      2 1 3      1      3    3
#> 3      3 1 4      1      4    3
#> 4      4 1 5      1      5    1
#> 5      5 2 3      2      3    1
#> 6      6 2 4      2      4    4
#> 7      7 2 5      2      5    1
#> 8      8 3 4      3      4    2
#> 9      9 3 5      3      5    1
#> 10    10 4 5      4      5    1
tda_g_flowcontrol(d)$table
#>    i j node_i node_j f1 f2 f3 f4 f5 f6
#> 1  1 2      1      2  1 -1 -1 -1 -1 -1
#> 2  1 3      1      3  3 -1  2 -1 -1 -1
#> 3  1 4      1      4  3 -1  2  1 -1 -1
#> 4  1 5      1      5  1 -1  1  1  0 -1
#> 5  2 3      2      3  1 -1 -1 -1 -1 -1
#> 6  2 4      2      4  4 -1 -1  3 -1 -1
#> 7  2 5      2      5  1 -1 -1  1  0 -1
#> 8  3 4      3      4  2 -1 -1 -1 -1 -1
#> 9  3 5      3      5  1 -1 -1 -1  0 -1
#> 10 4 5      4      5  1 -1 -1 -1 -1 -1
tda_g_backward(d)$table
#>   i node_i n_a node_j a_ji n_y node_k y_ki
#> 1 2      2   1      1    1   1      1    1
#> 2 3      3   2      1    2   1      1    3
#> 3 3      3   2      2    1   1     -1    0
#> 4 4      4   2      2    3   1      1    5
#> 5 4      4   2      3    2   1     -1    0
#> 6 5      5   1      4    1   1      1    1
tda_g_degrees(d)[, c("node", "in_degree", "out_degree")]
#>   node in_degree out_degree
#> 1    1         0          2
#> 2    2         1          2
#> 3    3         2          1
#> 4    4         2          1
#> 5    5         1          0
tda_g_edges(d)$table         # the raw edge dump (gdp)
#>   V1 V2 V3 V4 V5
#> 1  1  2  1  2  1
#> 2  1  3  1  3  2
#> 3  2  3  2  3  1
#> 4  2  4  2  4  3
#> 5  3  4  3  4  2
#> 6  4  5  4  5  1
tda_g_random(d)$table        # a fresh random graph, ignores e's edges
#>    V1 V2 V3
#> 1   1  2  1
#> 2   1  3  1
#> 3   1  4  1
#> 4   1  5  1
#> 5   1  6  1
#> 6   1  7  1
#> 7   1  8  1
#> 8   1  9  1
#> 9   1 10  1
#> 10  2  1  1
#> 11  2  3  1
#> 12  2  4  1
#> 13  2  5  1
#> 14  2  6  1
#> 15  2  7  1
#> 16  2  8  1
#> 17  2  9  1
#> 18  2 10  1
#> 19  3  1  1
#> 20  3  2  1
#> 21  3  4  1
#> 22  3  5  1
#> 23  3  6  1
#> 24  3  7  1
#> 25  3  8  1
#> 26  3  9  1
#> 27  3 10  1
#> 28  4  1  1
#> 29  4  2  1
#> 30  4  3  1
#> 31  4  5  1
#> 32  4  6  1
#> 33  4  7  1
#> 34  4  8  1
#> 35  4  9  1
#> 36  4 10  1
#> 37  5  1  1
#> 38  5  2  1
#> 39  5  3  1
#> 40  5  4  1
#> 41  5  6  1
#> 42  5  7  1
#> 43  5  8  1
#> 44  5  9  1
#> 45  5 10  1
#> 46  6  1  1
#> 47  6  2  1
#> 48  6  3  1
#> 49  6  4  1
#> 50  6  5  1
#> 51  6  7  1
#> 52  6  8  1
#> 53  6  9  1
#> 54  6 10  1
#> 55  7  1  1
#> 56  7  2  1
#> 57  7  3  1
#> 58  7  4  1
#> 59  7  5  1
#> 60  7  6  1
#> 61  7  8  1
#> 62  7  9  1
#> 63  7 10  1
#> 64  8  1  1
#> 65  8  2  1
#> 66  8  3  1
#> 67  8  4  1
#> 68  8  5  1
#> 69  8  6  1
#> 70  8  7  1
#> 71  8  9  1
#> 72  8 10  1
#> 73  9  1  1
#> 74  9  2  1
#> 75  9  3  1
#> 76  9  4  1
#> 77  9  5  1
#> 78  9  6  1
#> 79  9  7  1
#> 80  9  8  1
#> 81  9 10  1
#> 82 10  1  1
#> 83 10  2  1
#> 84 10  3  1
#> 85 10  4  1
#> 86 10  5  1
#> 87 10  6  1
#> 88 10  7  1
#> 89 10  8  1
#> 90 10  9  1
tda_g_dot(d)$text            # Graphviz source, not a table -- see tda_g()
#>  [1] "digraph g {"   "n1 [label=1];" "n2 [label=2];" "n3 [label=3];"
#>  [5] "n4 [label=4];" "n5 [label=5];" "n1 -> n2;"     "n1 -> n3;"    
#>  [9] "n2 -> n3;"     "n2 -> n4;"     "n3 -> n4;"     "n4 -> n5;"    
#> [13] "}"            
tda_g_shortest(g)$table      # all-pairs shortest paths
#>   index node to1 dist1 to2 dist2 to3 dist3 to4 dist4
#> 1     1    1   2     1   3     2   4     4   5     5
#> 2     2    2   1     1   3     1   4     3   5     4
#> 3     3    3   1     2   2     1   4     2   5     3
#> 4     4    4   1     4   2     3   3     2   5     1
#> 5     5    5   1     5   2     4   3     3   4     1

# a denser graph, with an actual triangle and edges valued in [0, 1] --
# what gcliques (a real clique, not just a connected pair) and gio (an
# ownership fraction) both need to find something on this small example
e2 <- data.frame(from = c(1, 2, 3, 1, 4), to = c(2, 3, 1, 4, 5),
                 value = c(0.3, 0.5, 0.8, 0.6, 0.9))
g2 <- tda_graph(e2, directed = FALSE)
d2 <- tda_graph(e2, directed = TRUE)
tda_g_gcliques(g2)$table
#>   V1 V2 V3 V4 V5
#> 1  1  3  1  2  3
#> 2  2  3  1  2  4
tda_g_ownership(d2)$table
#>    i node_i node_j a_ij  y_ij c iterations
#> 1  1      1      1  0.0 0.120 5          1
#> 2  1      1      2  0.3 0.300 5          1
#> 3  1      1      3  0.0 0.150 5          1
#> 4  1      1      4  0.6 0.600 5          1
#> 5  1      1      5  0.0 0.540 5          1
#> 6  2      2      1  0.0 0.400 5          2
#> 7  2      2      2  0.0 0.120 5          2
#> 8  2      2      3  0.5 0.500 5          2
#> 9  2      2      4  0.0 0.240 5          2
#> 10 2      2      5  0.0 0.216 5          2
#> 11 3      3      1  0.8 0.800 5          1
#> 12 3      3      2  0.0 0.240 5          1
#> 13 3      3      3  0.0 0.120 5          1
#> 14 3      3      4  0.0 0.480 5          1
#> 15 3      3      5  0.0 0.432 5          1
#> 16 4      4      4  0.0 0.000 2          1
#> 17 4      4      5  0.9 0.900 2          1
#> 18 5      5      5  0.0 0.000 1          1

# aggregate collapses named nodes into new ones -- rcn=, TDA's
# recode syntax, is required; without it there is nothing to do
tda_g_aggregate(d, rcn = "1[1,2],2[3,4,5]")$table
#>   V1 V2 V3 V4 V5
#> 1  1  1  1  1  1
#> 2  1  2  1  2  6
#> 3  2  2  2  2  3

# neighbourhoods wants an unvalued graph (gcni's restriction)
e3 <- data.frame(from = c(1, 1, 2, 2, 3, 4), to = c(2, 3, 3, 4, 4, 5))
tda_g_neighbourhoods(tda_graph(e3, directed = FALSE))$table
#> 3 records, lengths 5-7
#> [1]  1 3 3 1 1 
#> [2]  2 4 5 3 2 3 4 
#> [3]  3 2 1 1 5 

# this graph has no cycles (it is a DAG) and every edge already agrees
# with itself in one direction, so these two are empty rather than
# broken -- not every command has something to report
tda_g_dcycles(d)$table
#> NULL
tda_g_symmetric(d)$run$output[grepl("Nodes|symmetr", d$run$output)]
#> character(0)
tda_g_dblocks(d)$table  # 0 compact blocks in this graph, also real
#> NULL
tda_g_eigen(g)          # small graphs like these often do not converge
#> Call: tda_g_eigen(g)
#> 
#> Cases: 6 
#> 

# anything without a named wrapper goes through tda_g()
tda_g(d, "gdp")$table
#>   V1 V2 V3 V4 V5
#> 1  1  2  1  2  1
#> 2  1  3  1  3  2
#> 3  2  3  2  3  1
#> 4  2  4  2  4  3
#> 5  3  4  3  4  2
#> 6  4  5  4  5  1
```
