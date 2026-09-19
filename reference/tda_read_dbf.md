# Read a dBase (.dbf) attribute table

`tda_read_dbf` reads a dBase file through TDA's `rdbf` and returns it as
a data frame. Character fields come back as character, numeric ones as
numeric, taken from the field types the file itself declares.

## Usage

``` r
tda_read_dbf(file, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- file:

  path to a `.dbf` file.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame; `attr(x, "run")` carries the run.

## Details

This is the attribute half of a shapefile.
[`tda_read_shapefile`](tda_read_shapefile.md) reads the geometry; the
`.dbf` beside it holds one row per shape with whatever columns the data
has, which for real data is usually where the names live.
[`tda_read_shapefile`](tda_read_shapefile.md) attaches it as
`attr(x, "attributes")`.

## See also

Other spatial analysis: [`plot.tda_spatial()`](plot.tda_spatial.md),
[`tda_map()`](tda_map.md), [`tda_polygons()`](tda_polygons.md),
[`tda_read_shapefile()`](tda_read_shapefile.md),
[`tda_read_spatial`](tda_read_spatial.md), [`tda_sd()`](tda_sd.md),
[`tda_sd_analyses`](tda_sd_analyses.md),
[`tda_sd_ps3()`](tda_sd_ps3.md), [`tda_spatial()`](tda_spatial.md)

## Examples

``` r
# a shapefile's .dbf, which sf ships alongside nc.shp
if (requireNamespace("sf", quietly = TRUE)) {
  f <- sub("[.]shp$", ".dbf", system.file("shape/nc.shp", package = "sf"))
  head(tda_read_dbf(f)[, c("NAME", "AREA")])
}
#>          NAME  AREA
#> 1        Ashe 0.114
#> 2   Alleghany 0.061
#> 3       Surry 0.143
#> 4   Currituck 0.070
#> 5 Northampton 0.153
#> 6    Hertford 0.097
```
