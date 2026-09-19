# Read other spatial formats

The same shape as
[`tda_read_shapefile`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md):
TDA reads the file and writes its spatial data file, which the rest of
the family works on. `tda_read_e00` reads an uncompressed ARC/INFO
export, `tda_read_gshhs` a GSHHS coastline file, and `tda_read_dcw` a
Digital Chart of the World polygon point file.

## Usage

``` r
tda_read_e00(file, dir = tempfile("tda"), ...)

tda_read_gshhs(
  file,
  level = NULL,
  centre = c("greenwich", "pacific"),
  dir = tempfile("tda"),
  ...
)

tda_read_dcw(file, dir = tempfile("tda"), ...)
```

## Arguments

- file:

  path to the input file.

- dir:

  working directory for the run.

- ...:

  further options for the command. `sde00` takes `attr` for how many
  attribute variables to keep; `sdgshhs` takes `level` to select
  polygons; `sddcwp` takes `nc` to drop centroids.

- level:

  for `tda_read_gshhs`: which hierarchy levels to keep – 1 land, 2 lake,
  3 island in lake, 4 pond in island in lake, and since GSHHG 2.3.0 5
  the Antarctic ice front and 6 the Antarctic grounding line. Several
  may be given, `level = c(1, 2)` for land and lakes. `NULL`, the
  default, keeps all of them; level 1 alone is much the fastest way to
  draw a coastline.

  **Antarctica appears twice** in a 2.3.x file, once as level 5 and once
  as level 6, and they are about six degrees of latitude apart near the
  dateline. Keeping all levels draws both, which is why a world map
  shows two near-parallel southern coastlines; GSHHG intends you to pick
  one. Neither is an artefact – the ice front runs along roughly 78
  degrees south there (the Ross Ice Shelf edge, near-straight) and the
  grounding line along roughly 84.

- centre:

  for `tda_read_gshhs`: where the world map is cut. `"greenwich"`
  (default) gives the familiar view, longitudes -180..180, cut at the
  dateline. `"pacific"` cuts at the Greenwich meridian instead,
  longitudes 0..360. Polygons spanning whichever seam is in force are
  split either way.

## Value

A `tda_spatial` structure.

## Details

All three target older formats from GIS's pre-shapefile era, and only
shapefile really remains in wide, current use – e00 and DCW were both
effectively superseded by it, and DCW's distribution site has been gone
for years, which is worth knowing going in:

- `tda_read_e00` works, tested against a real ARC/INFO export (a point
  coverage).

- `tda_read_dcw` works, tested against a `.pnt` file rebuilt from a real
  DCW coverage read via the sf package (`st_read()` on the original
  coverage directory, coordinates regrouped into rings and written out
  in the plain text form `sddcwp` expects). Coordinates round-trip
  exactly: the intermediate this reader writes uses `fmt = "24.16"`, not
  TDA's default, which would cost six digits of each one.

- `tda_read_gshhs` reads both GSHHS header layouts: the pre-2.2 one and
  GSHHG 2.x, which is what you would download today. The layout is
  detected from the file itself and the version reported in the run's
  output.

  **Straight lines across a world map** would come from the dateline;
  that is handled for you, and there is nothing to configure. GSHHS
  stores the polygons that span the dateline with longitudes running
  past 180 (Eurasia reaches 190). Those are cut at the dateline and come
  back as separate objects, so each piece is drawn where it belongs and
  nothing is drawn between them; a polygon that does not span it is
  untouched. The Antarctic coastline, which runs from one edge of the
  map to the other, is returned as a line rather than a closed area for
  the same reason – closing it would draw a straight segment back across
  the world.

  This is what the GSHHG maintainers do in their own shapefile
  distribution, where the dateline-straddling polygons – the Antarctic
  cap chief among them – are split into east and west components. The
  native binary files are left unsplit for the reader to handle.

  `centre = "pacific"` cuts the map at Greenwich instead, giving
  longitudes on 0..360; the splitting then happens at that seam rather
  than at the dateline. Either way no polygon is drawn across the map.

## Further options

Anything else these commands accept can be passed through `...` under
TDA's option name. `tda_help("sdshp")`, `tda_help("sde00")`,
`tda_help("sddcwp")` and `tda_help("sdgshhs")` print the full lists. One
worth knowing is `nc`, which tells the DCW reader that the input has no
comment lines.

`df=`, `dtda=` and `fmt=` are set by the reader itself: the first two
name the intermediate files it writes and reads back, and `fmt` is fixed
at `"24.16"` so that no coordinate is rounded on the way through.

## See also

Other spatial analysis:
[`plot.tda_spatial()`](https://janmarvin.github.io/TDA/reference/plot.tda_spatial.md),
[`tda_map()`](https://janmarvin.github.io/TDA/reference/tda_map.md),
[`tda_polygons()`](https://janmarvin.github.io/TDA/reference/tda_polygons.md),
[`tda_read_dbf()`](https://janmarvin.github.io/TDA/reference/tda_read_dbf.md),
[`tda_read_shapefile()`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md),
[`tda_sd()`](https://janmarvin.github.io/TDA/reference/tda_sd.md),
[`tda_sd_analyses`](https://janmarvin.github.io/TDA/reference/tda_sd_analyses.md),
[`tda_sd_ps3()`](https://janmarvin.github.io/TDA/reference/tda_sd_ps3.md),
[`tda_spatial()`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)

## Examples

``` r
# nc.shp ships with sf, so this one runs where sf is installed.
# The .dbf beside the geometry holds the attribute table -- county
# names and areas -- and is attached as $attributes.
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
  head(nc$attributes[, c("NAME", "AREA")])
}
#>          NAME  AREA
#> 1        Ashe 0.114
#> 2   Alleghany 0.061
#> 3       Surry 0.143
#> 4   Currituck 0.070
#> 5 Northampton 0.153
#> 6    Hertford 0.097

if (FALSE) { # \dontrun{
# the other readers need files that do not ship with any package
tda_read_e00("world.e00")
tda_read_gshhs("gshhs_f.b", level = 1)  # either header layout
tda_read_dcw("polygon.pnt", nc = TRUE)
} # }
```
