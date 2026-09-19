# Print all case pairs

TDA's `pdatr`: writes the kept variables for every ordered pair of
cases.

## Usage

``` r
tda_pdatr(data, keep = names(data), ...)
```

## Arguments

- data:

  data frame.

- keep:

  variables to keep (all by default).

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the output file read back as a data frame.

## Examples

``` r
tda_pdatr(data.frame(x = c(1, 2, 3)))
#>   i j x <NA>
#> 1 1 1 1    1
#> 2 1 2 1    2
#> 3 1 3 1    3
#> 4 2 1 2    1
#> 5 2 2 2    2
#> 6 2 3 2    3
#> 7 3 1 3    1
#> 8 3 2 3    2
#> 9 3 3 3    3
```
