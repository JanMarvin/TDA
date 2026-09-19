# Read a UCINET file

Reads a UCINET dataset through TDA's `rucinet`. UCINET stores a dataset
as a pair of files, a header `.##h` and the data `.##d`; name the one
you want read.

## Usage

``` r
tda_read_ucinet(
  file,
  n_columns = NULL,
  form = c("matrix", "edges"),
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- file:

  path to the UCINET file, usually the `.##d`.

- n_columns:

  number of columns. `NULL`, the default, leaves TDA to work it out.

- form:

  `"matrix"` (default) reads it as a square matrix; `"edges"` reads it
  as an edge list, keeping only entries greater than zero.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame of the values TDA read. `attr(x, "run")` carries the run.

## Details

`#` begins a comment in TDA's command language, so a file called
`something.##d` truncates the command that names it and TDA reports
"Command file ends with an incomplete command". The file is therefore
copied to a `#`-free name for the run – which is why this takes a path
rather than leaving the caller to discover that.

## Examples

``` r
# a Ucinet ##d file is packed 4-byte floats, row-major
M <- rbind(c(0, 2, 0), c(2, 0, 1), c(0, 1, 0))
f <- tempfile(fileext = ".##d")
con <- file(f, "wb")
writeBin(as.numeric(t(M)), con, size = 4); close(con)
tda_read_ucinet(f, n_columns = 3)
#>   V1 V2 V3
#> 1  0  2  0
#> 2  2  0  1
#> 3  0  1  0
tda_read_ucinet(f, n_columns = 3, form = "edges")
#>   from to value
#> 1    1  2     2
#> 2    2  1     2
#> 3    2  3     1
#> 4    3  2     1
```
