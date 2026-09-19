# Distance matrix from raw variables

`pdatd`: builds a case-by-case distance matrix directly from a data
matrix, the way [`dist`](https://rdrr.io/r/stats/dist.html) does –
Euclidean or city-block distance, the count of variables on which two
cases differ, or a dissimilarity index. The result feeds directly into
[`tda_cluster`](https://janmarvin.github.io/TDA/reference/tda_cluster.md)
or [`tda_mds`](https://janmarvin.github.io/TDA/reference/tda_mds.md).
Verified against [`dist`](https://rdrr.io/r/stats/dist.html) for
Euclidean and city-block (Manhattan) distance: exact match.

## Usage

``` r
tda_pdatd(
  data,
  measure = 1,
  variables = NULL,
  weights = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- data:

  a data frame or matrix, one row per case.

- measure:

  the distance measure: `1` Euclidean (default), `2` city-block, `3`
  number of variables on which the two cases differ, `4` a dissimilarity
  index.

- variables:

  optional subset of `data`'s column names to use, instead of all of
  them.

- weights:

  optional positive per-variable weights, one per (selected) variable,
  all distinct – see Details for TDA's parser constraint and the exact
  weighting convention.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).
  `opt=` (the command's variant switch) passes through as well.

## Value

A [`dist`](https://rdrr.io/r/stats/dist.html) object.

## Details

As with
[`tda_dma`](https://janmarvin.github.io/TDA/reference/tda_dma.md), TDA's
`v=` needs variable names that start uppercase to be referenced
correctly; `variables` handles this automatically.

Per-variable weights work through `weights`, with one genuine TDA quirk
absorbed for you: `wt=` is parsed by TDA's *time-points* parser, which
accepts only a strictly increasing positive list – `wt=2,1,1` is a
syntax error while `wt=1,2,4` is fine, which is why this option long
looked dead. Since a distance does not care about variable order, the
wrapper sends the variables sorted by ascending weight. The weighting
convention, verified against R directly: euclidean uses \\d^2 = \sum_j
w_j \Delta_j^2\\ (so it matches `dist(sweep(x, 2, sqrt(w), "*"))`),
city-block \\d = \sum_j w_j \|\Delta_j\|\\. *Tied* weights cannot pass
TDA's parser at all; scale those columns yourself (\\\sqrt{w} x\\ for
euclidean, \\w x\\ for city-block) and call this unweighted – the result
is identical.

## See also

Other clustering:
[`TDA_CLUSTER`](https://janmarvin.github.io/TDA/reference/tda_cluster.md),
[`TDA_MDS`](https://janmarvin.github.io/TDA/reference/tda_mds.md),
[`tda_conjoint()`](https://janmarvin.github.io/TDA/reference/tda_conjoint.md),
[`tda_cutree()`](https://janmarvin.github.io/TDA/reference/tda_cutree.md),
[`tda_dma()`](https://janmarvin.github.io/TDA/reference/tda_dma.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x1 = rnorm(5), x2 = rnorm(5))
tda_pdatd(d)
#>          V1       V2       V3       V4
#> V2 1.538458                           
#> V3 1.572765 1.049697                  
#> V4 2.624046 1.414400 2.436338         
#> V5 1.085896 0.806124 1.564251 1.542284
tda_pdatd(d, measure = 2)   # city-block, not the euclidean default
#>           V1        V2        V3        V4
#> V2 2.1179946                              
#> V3 1.7679679 1.2701676                    
#> V4 3.6179843 1.4999898 2.5934528          
#> V5 1.4710416 0.9386819 2.2088495 2.1469428
tda_pdatd(d, variables = "x1")   # a subset of the variables
#>           V1        V2        V3        V4
#> V2 0.8100971                              
#> V3 0.2091748 1.0192719                    
#> V4 2.2217346 1.4116375 2.4309094          
#> V5 0.9559616 0.1458644 1.1651364 1.2657730
```
