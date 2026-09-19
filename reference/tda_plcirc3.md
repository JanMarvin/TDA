# Circle in a 3-d plot

TDA's `plcirc3` on a `psetup3` system.

## Usage

``` r
tda_plcirc3(radius, center = c(0, 0, 0), normal = c(0, 0, 1), ...)
```

## Arguments

- radius:

  circle radius.

- center:

  xyz center.

- normal:

  direction vector of the circle plane.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
invisible(tda_plcirc3(1, center = c(0, 0, 0), normal = c(0, 0, 1)))
```
