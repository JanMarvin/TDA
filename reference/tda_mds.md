# Multidimensional scaling

Represent a dissimilarity matrix as points in a low-dimensional space,
so that the distances between the points approximate the
dissimilarities.

## Usage

``` r
TDA_MDS

tda_mds(
  d,
  method = "classical",
  ndim = 2,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- d:

  a `dist`, or a square dissimilarity matrix.

- method:

  one of the names in `TDA_MDS`, matched partially.

- ndim:

  number of dimensions to return.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying `points`, the coordinates (one row per object, one
column per dimension), and, for `"classical"`, `eigenvalues` (a data
frame of each eigenvalue and its percentage of the total); for the
iterative methods, `stress`, the criterion value at the returned
configuration.

## Details

`method = "classical"` (TDA's `mdsc`) is the eigenvalue decomposition
also known as principal coordinates analysis – the same method as
[`cmdscale`](https://rdrr.io/r/stats/cmdscale.html). TDA computes the
full solution, one dimension per positive eigenvalue, and never reads a
dimension count of its own; because the classical solution is nested,
the first `ndim` columns of the full solution *are* the
`ndim`-dimensional solution, so this wrapper simply keeps those. The
eigenvalues, with the share of the total each accounts for, come back in
`eigenvalues` – the usual guide to how many dimensions the data support.

The other three methods fit the configuration iteratively. `"metric"`
(`mdsm`) minimizes the raw stress, the sum of squared differences
between dissimilarities and fitted distances, with TDA's general
minimizer; it fits two dimensions only, and errors for any other `ndim`.
`"nonmetric"` (`mdsn`) minimizes Kruskal's stress-1, preserving only the
rank order of the dissimilarities, by Kruskal's (1964) gradient method
with his step length adaptation; it reaches the same or lower stress
than [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html) on the usual
test data, but its default `mxit = 100` iterations is often not enough –
pass `options = list(mxit = 300)`, and `ns =` random restarts to guard
against local minima. `"minmax"` (`mdsx`) fits a minmax criterion; TDA
labels that command experimental and its output is left as it comes.
Check `stress` before trusting any iterative configuration.

## See also

Other clustering:
[`TDA_CLUSTER`](https://janmarvin.github.io/TDA/reference/tda_cluster.md),
[`tda_conjoint()`](https://janmarvin.github.io/TDA/reference/tda_conjoint.md),
[`tda_cutree()`](https://janmarvin.github.io/TDA/reference/tda_cutree.md),
[`tda_dma()`](https://janmarvin.github.io/TDA/reference/tda_dma.md),
[`tda_pdatd()`](https://janmarvin.github.io/TDA/reference/tda_pdatd.md)

## Examples

``` r
set.seed(34)
pts <- rbind(matrix(rnorm(10, 0), 5, 2), matrix(rnorm(10, 6), 5, 2))
m <- tda_mds(dist(pts), method = "classical", ndim = 2)
m$points       # recovers the two clusters as two groups of coordinates
#>         dim1        dim2
#> 1   4.670634  0.33670914
#> 2   2.922285  0.39611247
#> 3   5.386370  0.22155993
#> 4   3.975889 -1.08635946
#> 5   4.463387 -0.07961882
#> 6  -4.819784 -1.62779302
#> 7  -4.758347  0.11665206
#> 8  -3.763238 -0.12215922
#> 9  -3.755112  1.71887876
#> 10 -4.322085  0.12601814
m$eigenvalues  # two dominant dimensions, as built
#>         value   percent
#> 1  187.912506 96.332262
#> 2    7.154548  3.667738
#> 3    0.000000  0.000000
#> 4    0.000000  0.000000
#> 5    0.000000  0.000000
#> 6    0.000000  0.000000
#> 7    0.000000  0.000000
#> 8    0.000000  0.000000
#> 9    0.000000  0.000000
#> 10   0.000000  0.000000
```
