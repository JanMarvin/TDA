# Run a spatial command

`tda_sd` reaches any of TDA's spatial commands; the named functions are
the ones with an obvious meaning. The spatial data is opened with
`sdnvar` first, which is what the commands operate on.

## Usage

``` r
tda_sd(s, cmd, options = list(), rhs = NULL, setup = character(), ...)
```

## Arguments

- s:

  a
  [`tda_spatial`](https://janmarvin.github.io/TDA/reference/tda_spatial.md).

- cmd:

  the command name, e.g. `"sdinf"`, `"sdvd"`.

- options:

  a named list of options for the command.

- rhs:

  the command's right-hand side, when it takes one.

- setup:

  extra commands to run before it, such as a PostScript setup.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying the run and, where the command writes one, a `table`.

## See also

Other spatial analysis:
[`plot.tda_spatial()`](https://janmarvin.github.io/TDA/reference/plot.tda_spatial.md),
[`tda_map()`](https://janmarvin.github.io/TDA/reference/tda_map.md),
[`tda_polygons()`](https://janmarvin.github.io/TDA/reference/tda_polygons.md),
[`tda_read_dbf()`](https://janmarvin.github.io/TDA/reference/tda_read_dbf.md),
[`tda_read_shapefile()`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md),
[`tda_read_spatial`](https://janmarvin.github.io/TDA/reference/tda_read_spatial.md),
[`tda_sd_analyses`](https://janmarvin.github.io/TDA/reference/tda_sd_analyses.md),
[`tda_sd_ps3()`](https://janmarvin.github.io/TDA/reference/tda_sd_ps3.md),
[`tda_spatial()`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)

## Examples

``` r
d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
s <- tda_spatial(d)
r <- tda_sd(s, "sdinf")
cat(r$run$output, sep = "\n")
#> TDA. Analysis of Transition Data (6.4q). Sat Sep 19 20:00:52 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: commands
#> ============================================================================
#> sdnvar(...)
#> Reading a spatial data file. Current memory: 390032 bytes.
#> 
#> Idx Variable  T   S  PFmt  Definition
#> -------------------------------------
#>   1 SDID      3   4   0.0  c1
#>   2 SDTyp     3   4   0.0  c2
#>   3 SDN       3   4   0.0  c3
#>   4 SDPtr     3   5  11.0  rd
#> 
#> Creating a new data matrix.
#> Maximum number of cases: 1000
#> 
#> Using data file(s): spatial.sd
#> Free format. Separation character(s): default.
#> Reading data file: spatial.sd
#> Read records: 4 
#> 
#> Created a new data matrix.
#> Number of cases: 4
#> Number of variables: 4
#> Missing values in data file(s): none.
#> 
#> Number of points: 4
#> Number of lines: 0
#> Number of polygons: 0
#> Number of unknown objects: 0
#> 
#> Maximal number of points in spatial objects: 1
#> Total number of points in spatial objects: 4
#> Note: spatial.sd remains opened for further access.
#> 
#> X values. Minimum:         0.0000  Maximum:         1.0000
#> Y values. Minimum:         0.0000  Maximum:         1.0000
#> 
#> End of creating new variables. Current memory: 406449 bytes.
#> ----------------------------------------------------------------------------
#> Information about spatial data. Current memory: 406449 bytes.
#> 
#> Number of points: 4
#> XMin:     0.000000000000   YMin:     0.000000000000
#> XMax:     1.000000000000   YMax:     1.000000000000
#> 
#> Number of lines: 0
#> 
#> Number of polygons: 0
#> ----------------------------------------------------------------------------
#> Current memory: 390032 bytes. Max memory used: 426151 bytes.
#> End of program. Sat Sep 19 20:00:52 2026
```
