# Read a PLZ-style CSV file

TDA's `rplz`: the postcode-data variant of the semicolon CSV reader;
like [`tda_rcsv`](https://janmarvin.github.io/TDA/reference/tda_rcsv.md)
it renders to the protocol.

## Usage

``` r
tda_rplz(file, ...)
```

## Arguments

- file:

  the input file.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

a data frame, one row per record: id, type, place count, the three
leading fields, and the place list.

## Examples

``` r
# the PLZ export shape: a header line, then three semicolon fields
# of metadata before the place list (commas separate places, colons
# start sub-lists)
f <- tempfile()
writeLines(c("plz;land;kreis;orte",
             "01067;SN;Dresden;Dresden,Altstadt",
             "01069;SN;Dresden;Dresden,Suedvorstadt"), f)
tda_rplz(f)   # one row per record, places joined
#>   id type n_places   plz region district                places
#> 1  1    1        1 01067     SN  Dresden     Dresden; Altstadt
#> 2  2    1        1 01069     SN  Dresden Dresden; Suedvorstadt
```
