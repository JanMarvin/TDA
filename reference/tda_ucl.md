# Split an ordering into contiguous clusters

`tda_ucl` splits an ordering of a graph's nodes into `clusters`
contiguous groups so that the largest cluster diameter is as small as
possible, by the dynamic program of Alpert and Kahng (1997). The
diameter of a cluster is the largest edge value between any two of its
nodes.

## Usage

``` r
tda_ucl(
  g,
  clusters = 2,
  min_size = NULL,
  max_size = NULL,
  order = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- g:

  a
  [`tda_graph`](https://janmarvin.github.io/TDA/reference/tda_graph.md),
  undirected and valued.

- clusters:

  number of clusters, at least 2 (default 2).

- min_size, max_size:

  smallest and largest number of nodes a cluster may have. Defaults: 1,
  and the number of nodes.

- order:

  the ordering to split, as node indices `1..n` in the sequence they
  should be considered. `NULL`, the default, uses the graph's node
  order.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame with one row per node: `position` in the ordering, `node`
(internal index), `cluster`, and `label` (the node number the graph was
built with). `attr(x, "diameter")` is the largest cluster diameter, the
quantity minimised, and `attr(x, "run")` carries the run.

## Details

The clusters are runs of the *ordering*, not arbitrary subsets: this
answers "where do I cut this sequence" rather than "which things group
together". Give the ordering with `order`; the default is the nodes in
their internal order, which is rarely what you want unless the graph was
built that way. A seriation or a first principal coordinate is the usual
source.

The graph must be undirected and valued, and no edge value may be
negative or missing.

## A correctness note

TDA 6.4's `ucl` did not work. It optimised over the identity ordering
with its first two nodes swapped, counted a state with no legal split as
costing nothing, searched only the first few split positions, and read
past the end of its own cost matrix. On the example shipped with TDA it
returned a partition of diameter 7 where one of 4 exists. The number of
clusters and the size bounds were fixed at 3, 1 and 3 whatever you asked
for.

All of that is repaired here, and the result is checked in this
package's tests against a brute force over every contiguous partition.
The fixes are ours, not Rohwer's, and are listed in the package's
`README`.

## References

Alpert, C. J. and Kahng, A. B. (1997). Splitting an ordering into a
partition to minimize diameter. *Journal of Classification* **14**,
51–74.

## Examples

``` r
# six points on a line, distance growing with the gap between them
d <- as.dist(outer(1:6, 1:6, function(i, j) 3 * abs(i - j) + 1))
g <- tda_graph(d)
p <- tda_ucl(g, clusters = 3)
p
#>   position node cluster label
#> 1        1    1       1     1
#> 2        2    2       1     2
#> 3        3    3       2     3
#> 4        4    4       2     4
#> 5        5    5       3     5
#> 6        6    6       3     6
attr(p, "diameter")   # the largest cluster diameter, minimised
#> [1] 4
```
