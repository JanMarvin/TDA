# Read a GTOPO30 elevation profile

TDA's `gtopo`: reads a raw GTOPO30 digital-elevation tile (16-bit
big-endian integers, row-major from the upper-left corner) and extracts
a PROFILE – one latitude across a longitude range, or one longitude
across a latitude range. A rectangular window is not supported by the
command (established from the source: the one-value-plus-range shapes
are the only accepted forms).

## Usage

``` r
tda_gtopo(file, rows, cols, upper_left, pixel_size, lon, lat, ...)
```

## Arguments

- file:

  the raw DEM file.

- rows, cols:

  tile dimensions.

- upper_left:

  c(x, y) of the upper-left corner.

- pixel_size:

  c(dx, dy).

- lon, lat:

  the profile: one of them a single value, the other a c(from, to)
  range.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the profile records (a numeric matrix under direct exports, the file's
lines otherwise).

## Examples

``` r
f <- tempfile(fileext = ".dem")   # synthesize a tiny 3x4 tile
con <- file(f, "wb")
writeBin(as.integer(100 + 1:12), con, size = 2, endian = "big")
close(con)
tda_gtopo(f, rows = 3, cols = 4, upper_left = c(10, 50),
          pixel_size = c(0.5, 0.5), lat = 49.5, lon = c(10, 11.5))
#>      [,1] [,2] [,3]
#> [1,] 10.0 49.5  105
#> [2,] 10.5 49.5  106
#> [3,] 11.0 49.5  107
#> [4,] 11.5 49.5  108
```
