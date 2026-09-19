# Read a Stata file

Reads a Stata `.dta` file through TDA's `rstata`, which supports
releases 4, 6, 7 and 10 – the same set
[`tda_write_stata`](tda_write_foreign.md) writes, so a file written here
reads back here.

## Usage

``` r
tda_read_stata(
  file,
  n_records = NULL,
  missing = NULL,
  upper_names = FALSE,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- file:

  path to the `.dta` file.

- n_records:

  number of records to read. `NULL`, the default, reads all of them.

- missing:

  value to substitute for Stata's system missing. TDA's default is -5.

- upper_names:

  translate variable names to upper case (`rstata`'s `n=2`). Default
  `FALSE`.

- options:

  a named list of further TDA options, passed through. The reader also
  accepts, through `options`, TDA's `noc=`, `msys=`, `arcd=`, `dfa=`,
  `fn=`, `p=`, `pn=`, `zoo=` and `vdf=` switches.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame. `attr(x, "run")` carries the run.

## See also

[`tda_write_stata`](tda_write_foreign.md),
[`tda_read_spss`](tda_read_spss.md)

## Examples

``` r
d <- data.frame(id = 1:4, y = c(1.5, 2, 2.5, 3))
f <- tempfile(fileext = ".dta")
tda_write_stata(d, f)
tda_read_stata(f)
#>   Id   Y
#> 1  1 1.5
#> 2  2 2.0
#> 3  3 2.5
#> 4  4 3.0
```
