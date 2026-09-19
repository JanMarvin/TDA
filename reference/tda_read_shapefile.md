# Read a shapefile into a TDA spatial structure

`sdshp` reads an ESRI shapefile – the `.shp`, its `.shx` index and the
`.dbf` of attributes – and writes TDA's spatial data file, which is what
the rest of the spatial commands work on.

## Usage

``` r
tda_read_shapefile(file, dir = tempfile("tda"), ...)
```

## Arguments

- file:

  path to the shapefile, with or without the `.shp` extension. The
  `.shx` and `.dbf` beside it are read too.

- dir:

  working directory for the run.

- ...:

  further options for `sdshp`.

## Value

A `tda_spatial` object, as [`tda_spatial`](tda_spatial.md) returns.

## See also

Other spatial analysis: [`plot.tda_spatial()`](plot.tda_spatial.md),
[`tda_map()`](tda_map.md), [`tda_polygons()`](tda_polygons.md),
[`tda_read_dbf()`](tda_read_dbf.md),
[`tda_read_spatial`](tda_read_spatial.md), [`tda_sd()`](tda_sd.md),
[`tda_sd_analyses`](tda_sd_analyses.md),
[`tda_sd_ps3()`](tda_sd_ps3.md), [`tda_spatial()`](tda_spatial.md)

## Examples

``` r
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
  print(tda_sd_info(nc))
  plot(nc)
}
#> Call: tda_sd_info(nc)
#> 
#> Cases: 108 
#> 
#>  points lines polygons      xmin     ymin      xmax     ymax
#>       0     0      108 -84.32385 33.88199 -75.45698 36.58965
```
