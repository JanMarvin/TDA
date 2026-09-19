# Write a data frame to SPSS, Stata or a TDA system file

These build TDA's internal data matrix from `data` and then hand it to
TDA's writers. `tda_write_spss` produces a portable file (`.por`) or,
with `format = "sav"`, an SPSS system file. `tda_write_stata` produces a
`.dta` for a chosen Stata release. `tda_write_sys` produces a TDA system
file, which [`tda_read_sys`](tda_read_sys.md) reads back.

## Usage

``` r
tda_write_spss(
  data,
  file,
  keep = NULL,
  drop = NULL,
  sort = NULL,
  format = c("portable", "sav"),
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_write_stata(
  data,
  file,
  keep = NULL,
  drop = NULL,
  sort = NULL,
  release = 10,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_write_sys(data, file, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- data:

  a data frame.

- file:

  path to write to.

- keep, drop:

  variables to keep or drop, as a character vector. Give at most one of
  them.

- sort:

  variables to sort the cases by, as a character vector.

- format:

  for `tda_write_spss`: `"portable"` (the default, TDA's `wspss`) or
  `"sav"` (`wspss1`).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

- release:

  for `tda_write_stata`: the Stata release to write for – 4, 6, 7 or 10
  (default).

## Value

`file`, invisibly. `attr(x, "run")` carries the run.

## Details

A TDA system file is TDA's format: it stores the data matrix together
with the variable definitions, so reading one back restores the
session's variables without repeating the `nvar` block.

## Variable names

TDA will not read a variable name that does not begin with a capital –
`nvar(id = c1)` is a syntax error to it, and the same parser handles
`keep`, `drop` and `sort`. Columns are therefore capitalised on the way
in, so `id` is written as `Id`. Give `keep`, `drop` and `sort` the names
as they are in `data`; the capitalisation is applied for you.

## See also

[`tda_read_spss`](tda_read_spss.md), [`tda_read_sys`](tda_read_sys.md)

## Examples

``` r
d <- data.frame(id = 1:5, y = c(2.5, 3, 1.5, 4, 2))
f <- tempfile(fileext = ".por")
tda_write_spss(d, f)
tda_read_spss(f, portable = TRUE)
#>   Id   Y
#> 1  1 2.5
#> 2  2 3.0
#> 3  3 1.5
#> 4  4 4.0
#> 5  5 2.0
```
