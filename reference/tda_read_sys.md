# Read a TDA system file

Reads a system file written by [`tda_write_sys`](tda_write_foreign.md).
A system file carries the data matrix and the variable definitions
together, so the frame comes back with its column names without an
`nvar` block being repeated.

## Usage

``` r
tda_read_sys(file, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- file:

  path to the system file.

- options:

  a named list of further TDA options, passed through. The underlying
  `pdata` print accepts, through `options`, its selection and format
  switches: `keep=`, `drop=`, `sort=`, `noc=`, `nn=`, `nc=`, `nq=`,
  `ap=`, `l0=`, `sepc=`, `sd=` and `sdfmt=`.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame. `attr(x, "run")` carries the run.

## See also

[`tda_write_sys`](tda_write_foreign.md)

## Examples

``` r
d <- data.frame(id = 1:4, y = c(1.5, 2, 2.5, 3))
f <- tempfile(fileext = ".sys")
tda_write_sys(d, f)
tda_read_sys(f)
#>   Id   Y
#> 1  1 1.5
#> 2  2 2.0
#> 3  3 2.5
#> 4  4 3.0
```
