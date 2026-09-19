# Derive variables with TDA's expression language

Runs `nvar` on a data frame with new variables defined by TDA
expressions, and returns the frame with them added. This is the way to
the operators of the manual's 5.2.6 that have no R counterpart: the
change and count operators (`change`, `cntch`), the aggregates over
blocks (`gcnt`, `grec`, `gmean`, `gstd`, ...), TDA's random numbers, and
so on.

## Usage

``` r
tda_derive(data, ..., block = NULL, options = list(), dir = tempfile("tda"))
```

## Arguments

- data:

  a data frame.

- ...:

  named TDA expressions, as strings: `CH = "change(S)"`.

- block:

  a column name: with it, `nvar` runs in block mode (`dblock=`), and the
  block operators aggregate within blocks of consecutive rows with the
  same value.

- options:

  a named list of further `nvar` options.

- dir:

  working directory.

## Value

`data` with the new columns, from the export channel.

## Examples

``` r
d <- data.frame(ID = c(1, 1, 2), S = c(1, 3, 2))
tda_derive(d, CH = "change(S)", CN = "cntch(S)")
#>   ID S CH CN
#> 1  1 1  1  3
#> 2  1 3  1  3
#> 3  2 2  1  3
tda_derive(d, CN = "cntch(S)", block = "ID")
#>   ID S CN
#> 1  1 1  2
#> 2  1 3  2
#> 3  2 2  1
```
