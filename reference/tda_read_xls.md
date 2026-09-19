# Read an Excel .xls workbook

`tda_read_xls` reads a BIFF `.xls` workbook through TDA's `rxls` and
returns the numeric table it extracts. Every sheet is stacked into one
frame with a leading `sheet` column.

## Usage

``` r
tda_read_xls(
  file,
  skip_unknown = TRUE,
  numeric_only = FALSE,
  protocol = FALSE,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- file:

  path to a `.xls` file.

- skip_unknown:

  skip record types TDA does not know (default `TRUE`; see above).

- numeric_only:

  keep only columns TDA classed as numerical (TDA's `ns=1`).

- protocol:

  also return the workbook's string table, as `attr(x, "strings")`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame; `attr(x, "run")` carries the run.

## Details

**Unknown record types are skipped by default** (`skip_unknown = TRUE`,
TDA's `ni=1`). Excel 2007 and later write records into BIFF8 that TDA's
table does not list – `XFCRC` (`0x087c`), `STYLEEXT` (`0x0892`) and
others – so with TDA's default any `.xls` saved this century stops at
“unknown record type” before reaching a single cell. A BIFF record
carries its length, so skipping one is safe; naming them individually is
not, because the list keeps growing. Pass `skip_unknown = FALSE` for
TDA's original behaviour.

**Text columns come back as codes.** A column of text becomes an integer
index into the workbook's string table, not the text. TDA reports which
columns those are (the “Column / New / Type” table in the run's output),
and `protocol = TRUE` writes the strings themselves. This is TDA's
design, not a fault, but it means a column can come back as plausible
small integers that are codes: the `cyl` column of readxl's
`datasets.xls` is stored as text and arrives as 1, 2, 3 where the values
are 4, 6, 8. For anything where that matters, readxl reads the workbook
directly.

## See also

Other data:
[`tda_read_spss()`](https://janmarvin.github.io/TDA/reference/tda_read_spss.md),
[`tda_strings()`](https://janmarvin.github.io/TDA/reference/tda_strings.md)

## Examples

``` r
# readxl ships these workbooks, so this runs where readxl is installed
if (requireNamespace("readxl", quietly = TRUE)) {
  x <- tda_read_xls(readxl::readxl_example("datasets.xls"))
  # three sheets of different widths, stacked with a `sheet` column;
  # each is also available with its own column types
  head(attr(x, "sheets")[[1]])
}
#>     V1 V2  V3  V4   V5    V6    V7 V8 V9 V10 V11
#> 1 21.0  6 160 110 3.90 2.620 16.46  0  1   4   4
#> 2 21.0  6 160 110 3.90 2.875 17.02  0  1   4   4
#> 3 22.8  4 108  93 3.85 2.320 18.61  1  1   4   1
#> 4 21.4  6 258 110 3.08 3.215 19.44  1  0   3   1
#> 5 18.7  8 360 175 3.15 3.440 17.02  0  0   3   2
#> 6 18.1  6 225 105 2.76 3.460 20.22  1  0   3   1
```
