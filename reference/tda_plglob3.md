# Globe in a 3-d plot

TDA's `plglob3`: longitude/latitude grid on a sphere.

## Usage

``` r
tda_plglob3(radius, lon = seq(-150, 150, 30), lat = seq(-60, 60, 30), ...)
```

## Arguments

- radius:

  sphere radius.

- lon, lat:

  grid sequences (degrees).

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
invisible(tda_plglob3(1, lon = c(-90, 0, 90), lat = c(-45, 0, 45)))
```
