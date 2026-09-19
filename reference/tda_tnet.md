# Temporal-network summaries

TDA's `tnet` on a directed edge list.

## Usage

``` r
tda_tnet(
  edges,
  report = c("basic", "subgraphs", "indegrees", "outdegrees", "layers"),
  ...
)
```

## Arguments

- edges:

  data frame: from, to, then at least three value columns in the order
  VALIDITY (\>= 0 keeps the edge), then the edge's valid-from and
  valid-to dates with from \<= to – established from the source's check;
  the command requires a multigraph with three or more subgraphs.

- report:

  "basic", "subgraphs", "indegrees", "outdegrees", or "layers".

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the printed output.

## Examples

``` r
r <- tda_tnet(data.frame(i = c(1, 2, 3), j = c(2, 3, 4),
                         valid = c(1, 1, 1),
                         from = c(10, 12, 14), to = c(20, 22, 24)))
cat(head(tda_payload(r), 12), sep = "\n")
#> tnet(opt=1)
#> Temporal networks. Current memory: 410396 bytes.
#> Valid temporal network.
#> Number of nodes: 4
#> Number of edges: 3
#> Time axis: 10 to 24
#> Interpretation: temporal
```
