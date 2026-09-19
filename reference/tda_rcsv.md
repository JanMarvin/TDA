# Read a CSV file through TDA

TDA's `rcsv`: parses a csv file and writes a whitespace table. Verified
in the suite as a round trip against `read.csv`.

## Usage

``` r
tda_rcsv(file, sep = ";", ...)
```

## Arguments

- file:

  csv input. TDA's reader splits on SEMICOLONS (the European dialect);
  with `sep = ","` the wrapper converts a comma file on the way in.

- sep:

  the input file's separator.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

a data frame of TDA's parse when the export layer is on; otherwise the
rendered lines.

With the export layer on (the default), the engine hands back each raw
semicolon-bounded line of the converted file and the result is TDA's
parse as a data frame, first line as header. Without it, only TDA's
padded rendering exists; it carries no separator guarantee (adjacent
one-character entries can touch), so the rendered lines are returned as
they are. The engine accepts a `df=` option for this command but does
not write to it.

## Examples

``` r
f <- tempfile(fileext = ".csv")
write.csv(data.frame(a = 1:3, b = c(2.5, 1, 4)), f, row.names = FALSE)
tda_rcsv(f, sep = ",")   # a data frame under the export layer
#>   a   b
#> 1 1 2.5
#> 2 2 1.0
#> 3 3 4.0
```
