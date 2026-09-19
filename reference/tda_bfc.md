# Context dependencies in Boolean functions

TDA's `bfc`: Y as a Boolean function of binary X's.

## Usage

``` r
tda_bfc(y, x, undefined = c("dontcare", "as_one", "as_zero"), ...)
```

## Arguments

- y:

  binary outcome.

- x:

  binary matrix.

- undefined:

  how undefined arguments count: "dontcare", "as_one", or "as_zero".

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the printed output (implicants).

## Examples

``` r
x <- expand.grid(a = 0:1, b = 0:1)
r <- tda_bfc(y = as.integer(x$a & x$b), x = x)
cat(head(tda_payload(r), 12), sep = "\n")
#> bfc(...)=...
#> Context-dependencies in Boolean functions. Current memory: 402085 bytes.
#> Number of arguments: 2
#> Rows of reduced truth table: 2
#> Variable: X1
#> Not known: none
#> Positive effect: X2
#> Negative effect: none
#> No effect: X2'
#> Variable: X2
#> Not known: none
#> Positive effect: X1
```
