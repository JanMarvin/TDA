# Maps

`tda_map` draws spatial objects in a geographical coordinate system – a
projection, rather than the plain x and y of
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). TDA builds
that system with `psetupg` and then draws with `sdpmap`; `graticule`
adds meridians and parallels with `sdpgrat`, and `points` adds
coordinates given directly with `sdpgeo`.

## Usage

``` r
tda_map(
  s,
  view,
  region,
  projection = 10,
  width = 150,
  graticule = list(),
  points = NULL,
  symbol = 1,
  newpage = TRUE,
  ...
)
```

## Arguments

- s:

  a [`tda_spatial`](tda_spatial.md) structure.

- view:

  centre of the projection, as longitude and latitude.

- region:

  half-width and half-height of the mapped area, in degrees.

- projection:

  10 cylindrical equidistant, 11 cylindrical equal-area, 12 Mercator, 20
  azimuthal orthographic. For projection 20 the longitude half-width is
  limited to 90 degrees, and the whole region rectangle must fit on the
  projection disc, so both halves of `region` need
  `sin(lon)^2 + sin(lat)^2 <= 1` at the view latitude.

- width:

  width of the plot in mm.

- graticule:

  a list of options for the meridians and parallels, or `NULL` for none.
  `lon` and `lat` take TDA's sequence notation, `"-84(2)-76"` meaning
  every 2 degrees from -84 to -76.

- points:

  a data frame of `lon` and `lat` to mark, or `NULL`.

- symbol:

  marker symbol for `points`; TDA's `s=`, which cannot be spelled that
  way here because `s` is the structure. 4 is an open circle, 5 a filled
  one, 8 a square, 15 a diamond.

- newpage:

  start a new grid page.

- ...:

  further options for `sdpmap`.

## Value

The map, drawn; the session is returned invisibly.

## Details

`region` is required and has no default: it is the half-width and
half-height of the mapped area in degrees, measured from `view`. So
`view = c(-80, 35), region = c(6, 3)` maps 6 degrees of longitude and 3
of latitude either side of that point.

## See also

Other spatial analysis: [`plot.tda_spatial()`](plot.tda_spatial.md),
[`tda_polygons()`](tda_polygons.md),
[`tda_read_dbf()`](tda_read_dbf.md),
[`tda_read_shapefile()`](tda_read_shapefile.md),
[`tda_read_spatial`](tda_read_spatial.md), [`tda_sd()`](tda_sd.md),
[`tda_sd_analyses`](tda_sd_analyses.md),
[`tda_sd_ps3()`](tda_sd_ps3.md), [`tda_spatial()`](tda_spatial.md)

## Examples

``` r
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
  tda_map(nc, view = c(-80, 35), region = c(6, 3),
          graticule = list(lon = "-84(2)-76", lat = "34(1)37",
                           fsx = 2, fsy = 2))
}
```
