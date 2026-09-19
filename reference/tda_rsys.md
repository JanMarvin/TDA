# Read a TDA system file

TDA's `rsys`: restores a data matrix written by `wsys`. Verified in the
suite as a write/read round trip.

## Usage

``` r
tda_rsys(file, ...)
```

## Arguments

- file:

  the system file.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the restored data as a data frame, like
[`tda_read_sys`](https://janmarvin.github.io/TDA/reference/tda_read_sys.md),
which this calls.

## Examples

``` r
d0 <- data.frame(id = 1:3, x = c(1.5, 2, 3))
dr <- tempfile("tda"); dir.create(dr)
write.table(d0, file.path(dr, "d.dat"),
            row.names = FALSE, col.names = FALSE)
invisible(tda_run(c("nvar(dfile=d.dat, ID[4.0]=c1, X[8.2]=c2);",
                    "wsys = t.sys;"), dir = dr))
tda_rsys(file.path(dr, "t.sys"))   # the values restored
#>   ID   X
#> 1  1 1.5
#> 2  2 2.0
#> 3  3 3.0
```
