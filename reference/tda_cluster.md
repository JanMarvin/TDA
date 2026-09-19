# Clustering from a dissimilarity matrix

TDA's clustering commands, all of which take a dissimilarity and so
accept a `dist` – including the one
[`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md)
returns.

## Usage

``` r
TDA_CLUSTER

tda_cluster(
  d,
  method = "hierarchical",
  algorithm = c("centers", "diameter"),
  min_size = NULL,
  max_splits = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- d:

  a `dist`, or a square dissimilarity matrix.

- method:

  one of the names in `TDA_CLUSTER`, matched partially.

- algorithm:

  for `method = "hierarchical"` only: split by `"centers"` (default,
  maximally different cluster centers) or `"diameter"` (minimal-diameter
  partitions) – `hcld`'s `alg=`.

- min_size:

  for `algorithm = "centers"`, the minimum cluster size (`hcld`'s
  `min=`).

- max_splits:

  for `algorithm = "diameter"`, the maximum number of splits (`hcld`'s
  `max=`).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying a `table`, or for `"hierarchical"` `clusters` instead
(see Details). For `"hierarchical"` the object also carries `merges`:
one row per split, the two clusters, their sizes and their diameters.

## Details

The default, `"hierarchical"` (`hcld`), writes a ragged table – a node,
its member count, then that many member indices – which is not a
rectangular `data.frame`, so it is parsed separately into `clusters`, a
named list of member-index vectors, one per node of the divisive tree
(the root first, then each split). `"unidimensional"` and `"additive"`
report their result only as text (a final partition or a fitted tree),
never through `df=`; read `$run$output` for those.

## See also

Other clustering:
[`TDA_MDS`](https://janmarvin.github.io/TDA/reference/tda_mds.md),
[`tda_conjoint()`](https://janmarvin.github.io/TDA/reference/tda_conjoint.md),
[`tda_cutree()`](https://janmarvin.github.io/TDA/reference/tda_cutree.md),
[`tda_dma()`](https://janmarvin.github.io/TDA/reference/tda_dma.md),
[`tda_pdatd()`](https://janmarvin.github.io/TDA/reference/tda_pdatd.md)

## Examples

``` r
d <- dist(matrix(c(1, 2, 3, 8, 9, 10), ncol = 1))
tda_cluster(d)$clusters  # the divisive tree: root, then each split
#> $`1`
#> [1] 1 2 3 4 5 6
#> 
#> $`2`
#> [1] 1 2 3
#> 
#> $`3`
#> [1] 4 5 6
#> 
#> $`4`
#> [1] 1 2
#> 
#> $`5`
#> [1] 3
#> 
#> $`6`
#> [1] 4 5
#> 
#> $`7`
#> [1] 6
#> 
```
