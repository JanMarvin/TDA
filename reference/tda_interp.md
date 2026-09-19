# Interpolate a surface over a grid

Interpolate a surface over a grid

## Usage

``` r
tda_interp(x, y, z, rx, ry, options = list(), dir = tempfile("tda"))
```

## Arguments

- x, y, z:

  the surface: a point cloud in three variables.

- rx, ry:

  the grids to interpolate over.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying a `table`: one row per grid point, `x`, `y` and the
interpolated `z`.

## See also

Other smoothing:
[`tda_integrate()`](https://janmarvin.github.io/TDA/reference/tda_integrate.md),
[`tda_isotonic()`](https://janmarvin.github.io/TDA/reference/tda_isotonic.md),
[`tda_mat()`](https://janmarvin.github.io/TDA/reference/tda_mat.md),
[`tda_smd()`](https://janmarvin.github.io/TDA/reference/tda_smd.md),
[`tda_spl()`](https://janmarvin.github.io/TDA/reference/tda_spl.md)

## Examples

``` r
# a plane z = x + 2y, evaluated back on its own grid
g <- expand.grid(x = 1:5, y = 1:5)
g$z <- g$x + 2 * g$y
head(tda_interp(g$x, g$y, g$z, rx = 1:5, ry = 1:5)$table)
#>   x y  z
#> 1 1 1  3
#> 2 1 2  5
#> 3 1 3  7
#> 4 1 4  9
#> 5 1 5 11
#> 6 2 1  4
```
