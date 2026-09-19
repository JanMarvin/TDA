# Non-metric conjoint analysis

`nmca` is conjoint measurement: it takes a rank order over profiles and
the categorical factors describing them, and estimates the part-worth of
each factor level. It draws a Hasse diagram of the order as part of the
output.

## Usage

``` r
tda_conjoint(formula, data, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- formula:

  `rank ~ f1 + f2`, where the response is the rank order and the
  right-hand side names the factors.

- data:

  a data frame. Factor levels are numbered from 1 upwards, and a value
  outside that range is rejected by TDA with the record it occurred in.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object carrying the run; the diagram and the estimates are in
`$run$output`.

## See also

Other clustering: [`TDA_CLUSTER`](tda_cluster.md),
[`TDA_MDS`](tda_mds.md), [`tda_cutree()`](tda_cutree.md),
[`tda_dma()`](tda_dma.md), [`tda_pdatd()`](tda_pdatd.md)

## Examples

``` r
d <- expand.grid(price = 1:2, brand = 1:2)
d$rank <- c(4, 2, 3, 1)
cj <- tda_conjoint(rank ~ price + brand, d)
cat(tda_payload(cj$run), sep = "\n")   # the nmca section itself
#> nmca(...)=...
#> Non-metric conjoint analysis. Current memory: 390161 bytes.
#> Rank order variable: Vrank
#> Dimension  Variable   Categories
#>      1      Vprice       2 :  1  2 
#>      2      Vbrand       2 :  1  2 
#> Found: ( 2 , 1 ) ( 2 , 1 ) 
#> Found 1 partial orders.
```
