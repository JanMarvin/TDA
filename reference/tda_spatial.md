# Spatial data for TDA's spatial commands

Builds TDA's spatial data structure from a data frame of coordinates.
Each object – a point, a line or a polygon – is a group of rows sharing
an id.

## Usage

``` r
tda_spatial(
  data,
  id = "id",
  x = "x",
  y = "y",
  attributes = NULL,
  type = 1,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- data:

  a data frame with an id column and two coordinate columns.

- id, x, y:

  the column names.

- attributes:

  further columns to carry along, TDA's `av=`.

- type:

  type of spatial object, TDA's `typ=`: 1 for points.

- options:

  a named list of further `sdgen` options.

- dir:

  working directory; the spatial file is written there and stays for the
  commands that read it.

## Value

A `tda_spatial`.

## See also

Other spatial analysis: [`plot.tda_spatial()`](plot.tda_spatial.md),
[`tda_map()`](tda_map.md), [`tda_polygons()`](tda_polygons.md),
[`tda_read_dbf()`](tda_read_dbf.md),
[`tda_read_shapefile()`](tda_read_shapefile.md),
[`tda_read_spatial`](tda_read_spatial.md), [`tda_sd()`](tda_sd.md),
[`tda_sd_analyses`](tda_sd_analyses.md), [`tda_sd_ps3()`](tda_sd_ps3.md)

## Examples

``` r
d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
s <- tda_spatial(d)
tda_sd_info(s)
#> Call: tda_sd_info(s)
#> 
#> Cases: 4 
#> 
#>  points lines polygons xmin ymin xmax ymax
#>       4     0        0    0    0    1    1
```
