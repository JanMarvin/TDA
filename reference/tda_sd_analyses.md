# Spatial analyses

One function per spatial command, each a thin wrapper: the command name,
the options it takes and the table TDA writes. They all work on a
[`tda_spatial`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)
structure, however it was built – from a data frame with `tda_spatial`,
from a shapefile with
[`tda_read_shapefile`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md),
or from line segments with
[`tda_polygons`](https://janmarvin.github.io/TDA/reference/tda_polygons.md).

## Usage

``` r
tda_sd_info(s, ...)

tda_sd_relations(s, ...)

tda_sd_neighbours(s, ...)

tda_sd_enclosing(s, shape = c("rectangle", "hull"), ...)

tda_sd_voronoi(s, ...)

tda_sd_select(s, ...)

tda_sd_intersect(s, ...)

tda_sd_lines(s, ...)

tda_sd_points(s, ...)

tda_sd_data(s, ...)

tda_sd_clip(s, rec, ...)
```

## Arguments

- s:

  a
  [`tda_spatial`](https://janmarvin.github.io/TDA/reference/tda_spatial.md).

- ...:

  options for the command.

- shape:

  for `tda_sd_enclosing`, `"rectangle"` (default) or `"hull"` (the
  convex hull); `sdencl`'s `opt=`.

- rec:

  the clipping rectangle, `xmin,ymin,xmax,ymax`.

## Value

An object carrying the run and, where there is one, a table.

## Details

`tda_sd_info` reports what a structure holds: how many points, lines and
polygons, and the bounding box. `tda_sd_voronoi` builds a Voronoi
diagram and the Delaunay triangulation of the points; `tda_sd_enclosing`
the enclosing rectangle or the convex hull, `tda_sd_neighbours`
nearest-neighbour lists, and `tda_sd_relations` the relations between
objects. `tda_sd_select` picks out objects in a region – like
`tda_sd_clip`, it needs `rec` (as an option:
`tda_sd_select(s, rec = "0,0,1,1")`), which TDA will not guess. The
polygon operations – `tda_sd_intersect`, `tda_sd_lines`, `tda_sd_points`
and `tda_sd_clip` – each write a record file that comes back as a table.

## See also

Other spatial analysis:
[`plot.tda_spatial()`](https://janmarvin.github.io/TDA/reference/plot.tda_spatial.md),
[`tda_map()`](https://janmarvin.github.io/TDA/reference/tda_map.md),
[`tda_polygons()`](https://janmarvin.github.io/TDA/reference/tda_polygons.md),
[`tda_read_dbf()`](https://janmarvin.github.io/TDA/reference/tda_read_dbf.md),
[`tda_read_shapefile()`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md),
[`tda_read_spatial`](https://janmarvin.github.io/TDA/reference/tda_read_spatial.md),
[`tda_sd()`](https://janmarvin.github.io/TDA/reference/tda_sd.md),
[`tda_sd_ps3()`](https://janmarvin.github.io/TDA/reference/tda_sd_ps3.md),
[`tda_spatial()`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)

## Examples

``` r
# four points, as a spatial structure
d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
s <- tda_spatial(d)

# what is in it
tda_sd_info(s)
#> Call: tda_sd_info(s)
#> 
#> Cases: 4 
#> 
#>  points lines polygons xmin ymin xmax ymax
#>       4     0        0    0    0    1    1

# the Voronoi diagram of those points: what sdvd writes is itself a
# spatial data file (points for the Voronoi vertices, lines for the
# edges), so plot() reads it back with the input's plot machinery
v <- tda_sd_voronoi(s, opt = 2)
plot(v, xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5))


# clip to a rectangle; rec is required, TDA will not guess one
tda_sd_clip(s, c(0, 0, 0.5, 0.5))
#> Call: tda_sd_clip(s, c(0, 0, 0.5, 0.5))
#> 
#> Cases: 4 
#> 
#> 2 records, lengths 2-5
#> [1]  1 1 1 1 1 
#> [2]  0 0 

# a select needs one too, the region to select within
tda_sd_select(s, rec = "0,0,1,1")$table
#> 8 records, lengths 2-5
#> [1]  1 1 1 1 1 
#> [2]  0 0 
#> [3]  2 1 1 2 1 
#> [4]  1 0 
#> [5]  3 1 1 3 1 
#> [6]  1 1 
#> [7]  4 1 1 4 1 
#> [8]  0 1 

# relations and nearest neighbours, on the same points
tda_sd_relations(s)
#> Call: tda_sd_relations(s)
#> 
#> Cases: 4 
#> 
#> Elementary topological relations. Current memory: 406449 bytes.
#> 0 records written to: out.txt
tda_sd_neighbours(s)
#> Call: tda_sd_neighbours(s)
#> 
#> Cases: 4 
#> 
#> Creating a network from lines. Current memory: 406449 bytes.
#> Number of lines: 0
#> Number of segments: 0
#> Tolerance: 0.0001
#> Maximal number of segments: 0
#> Number of nodes: 0
#> Number of edges: 0
#> 0 records written to: out.txt
tda_sd_enclosing(s, shape = "hull")$table
#> 5 records, lengths 1-2
#> [1]  4 
#> [2]  0 1 
#> [3]  0 0 
#> [4]  1 0 
#> [5]  1 1 

# polygon data: typ = 3 groups consecutive rows sharing an id into one
# polygon each -- here two squares
poly <- data.frame(id = c(1, 1, 1, 1, 2, 2, 2, 2),
                   x  = c(0, 2, 2, 0, 3, 5, 5, 3),
                   y  = c(0, 0, 2, 2, 0, 0, 2, 2))
sp <- tda_spatial(poly, type = 3)
tda_sd_data(sp)$table   # the polygons' own vertex coordinates back out
#>   V1 V2 V3
#> 1  1  0  0
#> 2  1  2  0
#> 3  1  2  2
#> 4  1  0  2
#> 5  2  3  0
#> 6  2  5  0
#> 7  2  5  2
#> 8  2  3  2

# tda_sd_intersect/tda_sd_lines/tda_sd_points compare polygons against a
# *different* object type (points or lines) already in the same spatial
# structure -- with only the two polygons above and nothing else, there
# is nothing for them to find, which is what they correctly report
tda_sd_lines(sp)$table
#> NULL
tda_sd_intersect(sp)$table
#>   V1 V2 V3
#> 1  1  1 -1
#> 2  2  2 -1
tda_sd_points(sp)$table
#> NULL

# a real shapefile, if one is available: sf's example data, North
# Carolina counties -- everything above works on synthetic points or two
# squares; this is what these look like on real, complicated polygons
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
  print(tda_sd_info(nc))                 # 108 counties, 2421 vertices

  # relations' real answer is printed, not written to a file -- $text
  # holds the boundary comparisons TDA's console output shows
  head(tda_sd_relations(nc)$text, 4)

  # a neighbour network needs line data; a set of polygons alone
  # legitimately has none to report -- 0 nodes, 0 edges, not a bug
  tda_sd_neighbours(nc)$text

  tda_sd_enclosing(nc)$table             # the whole state's bounding box

  # select/clip write one record per vertex, not all the same length --
  # a polygon's header row, then that many coordinate pairs, the
  # same ragged shape as sdvd's edge list above -- a proper rectangular
  # table is what tda_sd_data() gives instead
  tda_sd_select(nc, rec = "-79,35,-77,36")$table   # counties in a region
  tda_sd_clip(nc, c(-79, 35, -77, 36))$table       # the same, cut to it
  head(tda_sd_data(nc)$table)            # every county's vertices,
                                         # as an ordinary data frame
}
#> Call: tda_sd_info(nc)
#> 
#> Cases: 108 
#> 
#>  points lines polygons      xmin     ymin      xmax     ymax
#>       0     0      108 -84.32385 33.88199 -75.45698 36.58965
#>   V1        V2       V3
#> 1  1 -81.47276 36.23436
#> 2  1 -81.54084 36.27251
#> 3  1 -81.56198 36.27359
#> 4  1 -81.63306 36.34069
#> 5  1 -81.74107 36.39178
#> 6  1 -81.69828 36.47178
```
