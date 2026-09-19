# Inspect a raw file

Small file-diagnostic utilities from TDA's toolkit, from an era before
assuming `wc`, `od` or `split` were available: `tda_ccnt` a character
frequency table, `tda_lcnt` a record-length frequency table, `tda_dump`
a hex dump, and `tda_dsplit` splits a file into fixed-size parts on
disk.

## Usage

``` r
tda_ccnt(file, dir = tempfile("tda"), ...)

tda_lcnt(file, noc = NULL, dir = tempfile("tda"), ...)

tda_dump(file, nc = NULL, offset = NULL, dir = tempfile("tda"), ...)

tda_dsplit(file, len = NULL, dir = tempfile("tda"), ...)
```

## Arguments

- file:

  path to the file to inspect.

- dir:

  working directory; for `tda_dsplit`, also where the parts are written,
  alongside the original file.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- noc:

  for `tda_lcnt`, the maximum number of records to read.

- nc:

  for `tda_dump`, how many bytes to show; all of them by default.

- offset:

  for `tda_dump`, where to start, in bytes.

- len:

  for `tda_dsplit`, the size of each part in bytes.

## Value

`tda_ccnt` and `tda_lcnt` return a data frame. `tda_dump` returns the
hex dump as a character vector, one line each. `tda_dsplit` returns the
paths of the parts it wrote.

## See also

Other TDA infrastructure:
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md),
[`tda_output()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md),
[`tda_read_table()`](https://janmarvin.github.io/TDA/reference/tda_read_table.md),
[`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md),
[`tda_write_data()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)

## Examples

``` r
f <- tempfile()
writeLines(c("hello world", "foo bar baz", "x"), f)
tda_ccnt(f)
#>    char hex count
#> 1  <NA>  0A     3
#> 2  <NA>  20     3
#> 3     a  61     2
#> 4     b  62     2
#> 5     d  64     1
#> 6     e  65     1
#> 7     f  66     1
#> 8     h  68     1
#> 9     l  6C     3
#> 10    o  6F     4
#> 11    r  72     2
#> 12    w  77     1
#> 13    x  78     1
#> 14    z  7A     1
tda_lcnt(f)
#>   index length frequency
#> 1     1      1         1
#> 2     2     11         2
cat(tda_dump(f), sep = "\n")
#> 00000000 : 68 65 6c 6c 6f 20 77 6f 72 6c 64 0a 66 6f 6f 20   hello world.foo 
#> 00000010 : 62 61 72 20 62 61 7a 0a 78 0a                     bar baz.x.

f2 <- tempfile()
writeLines(rep("0123456789", 10), f2)
tda_dsplit(f2, len = 20)
#> [1] "/tmp/RtmptmSFox/tda1ebb2771a566/file1ebb3c541183.a"
#> [2] "/tmp/RtmptmSFox/tda1ebb2771a566/file1ebb3c541183.b"
#> [3] "/tmp/RtmptmSFox/tda1ebb2771a566/file1ebb3c541183.c"
#> [4] "/tmp/RtmptmSFox/tda1ebb2771a566/file1ebb3c541183.d"
#> [5] "/tmp/RtmptmSFox/tda1ebb2771a566/file1ebb3c541183.e"
#> [6] "/tmp/RtmptmSFox/tda1ebb2771a566/file1ebb3c541183.f"
```
