# Cut a dendrogram into groups

Builds a dendrogram with TDA's `hcls` and cuts it with `hclsp`.

## Usage

``` r
tda_cutree(
  d,
  nlev = 2,
  nodes = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- d:

  a `dist`, or a square dissimilarity matrix.

- nlev:

  number of levels at which to cut.

- nodes:

  an explicit sequence of node numbers to cut at, TDA's `cn=`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying `merge`, the dendrogram, and `table`, one row per
cluster in the cut (`node`, `leaves`, `parent`).

## Details

Unlike R's `cutree`, which returns one group label per observation,
`hclsp` reports the cut as a set of nodes: `table` has one row per
resulting cluster, its TDA node number, how many leaves (original
observations) it holds, and its parent node – not a per-observation
vector. Getting a label for each observation means walking `merge` (the
dendrogram `hcls` built, one row per step) down from a node in `table`
to the leaves under it; this wrapper does not do that walk for you.

## See also

Other clustering:
[`TDA_CLUSTER`](https://janmarvin.github.io/TDA/reference/tda_cluster.md),
[`TDA_MDS`](https://janmarvin.github.io/TDA/reference/tda_mds.md),
[`tda_conjoint()`](https://janmarvin.github.io/TDA/reference/tda_conjoint.md),
[`tda_dma()`](https://janmarvin.github.io/TDA/reference/tda_dma.md),
[`tda_pdatd()`](https://janmarvin.github.io/TDA/reference/tda_pdatd.md)

## Examples

``` r
set.seed(35)
pts <- rbind(matrix(rnorm(10, 0), 5, 2), matrix(rnorm(10, 8), 5, 2))
ct <- tda_cutree(dist(pts), nlev = 2)
ct$table    # two clusters of five leaves each
#>   node leaves parent
#> 1   19     10     NA
#> 2   17      5     19
#> 3   18      5     19
ct$merge    # the nine merge steps that built the dendrogram
#>   step    height a  b
#> 1    1  0.800665 2  4
#> 2    2  1.112126 1  3
#> 3    3  1.358110 7  9
#> 4    4  1.921337 6  8
#> 5    5  1.934383 1  2
#> 6    6  2.489013 7 10
#> 7    7  3.492846 1  5
#> 8    8  4.215649 6  7
#> 9    9 12.841900 1  6
```
