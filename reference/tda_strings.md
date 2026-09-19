# String variables from a run

`tda_strings` returns the string variables a run printed with `pdata`,
as a named list of character vectors.

## Usage

``` r
tda_strings(res, trim = TRUE)
```

## Arguments

- res:

  a run, as returned by
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- trim:

  trim the fixed-width padding (default `TRUE`).

## Value

A named list of character vectors, one per string variable, or an empty
list if the run printed none.

## Details

TDA stores strings fixed-width, so every value is padded to the
variable's declared length; the padding is trimmed here, since it is
storage rather than data. The values come through the export channel
rather than being re-read from the printed file: a string cannot travel
on the numeric channel, and parsing it back out of fixed-width text is
exactly what the exports exist to avoid.

A character or factor column in a data frame becomes a TDA string
variable (see
[`tda_write_data`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)),
and `rspss` and `rstata` create them from the file's string columns.
Those are the only ways TDA makes one.

## See also

Other data:
[`tda_read_spss()`](https://janmarvin.github.io/TDA/reference/tda_read_spss.md),
[`tda_read_xls()`](https://janmarvin.github.io/TDA/reference/tda_read_xls.md)

## Examples

``` r
d <- data.frame(ID = 1:3, S = c("pear", "fig", "date"),
                X = c(1.5, 2.5, 3.5))
r <- tda_run(c(tda_nvar(d), "pdata() = out.txt;"), data = d)
tda_strings(r)
#> $S
#> [1] "pear" "fig"  "date"
#> 
```
