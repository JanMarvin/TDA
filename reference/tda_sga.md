# Scalogram (Guttman) analysis

TDA's `sga`: values \> 0 mean success on a task.

## Usage

``` r
tda_sga(x, ...)
```

## Arguments

- x:

  matrix, persons in rows, tasks in columns.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the printed output.

## Examples

``` r
r <- tda_sga(rbind(c(1, 1, 0), c(1, 0, 0), c(1, 1, 1)))
cat(head(tda_payload(r), 12), sep = "\n")
#>         1    2    3       3              0      0      0
#> Reproducibility
#> Err0: 1.0000
#> Err1: 1.0000
#> Err2: 1.0000
#> 0 1 0 
#> 0 0 1 
#> 0 0 0 
#> digraph g {
#> n1 [label=1];
#> n2 [label=2];
#> n3 [label=3];
```
