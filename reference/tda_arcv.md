# Print archive variables

TDA's `arcv` on an open archive.

## Usage

``` r
tda_arcv(out, data_file = NULL, archive = NULL, ...)
```

## Arguments

- out:

  output file name.

- data_file:

  optional data-file filter.

- archive:

  the archive description (`.zad`) to load first. `arcv` reads an
  archive that `arcd` has loaded, and TDA keeps no state between runs,
  so this is needed unless the caller issues `arcd` in the same run
  themselves.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

  Requires an open archive (see
  [`tda_arcd`](https://janmarvin.github.io/TDA/reference/tda_arcd.md))
  containing a type-2 (variable description) member.

## Examples

``` r
d <- tempfile(); dir.create(d)
for (f in c("tda.zad", "tda.zoo"))
    file.copy(system.file("extdata", f, package = "tdaR"), d)
tda_arcv("vars.vd", archive = "tda.zad", dir = d)
readLines(file.path(d, "vars.vd"))
#> [1] "# arcd = tda.zad;"                                 
#> [2] "# nvar( "                                          
#> [3] "# V1<5>[8.0] = A:V1, # [adata.dat] first variable" 
#> [4] "# V2<5>[8.0] = A:V2, # [adata.dat] second variable"
#> [5] "# V3<5>[8.0] = A:V3, # [adata.dat] third variable" 
#> [6] "# );"                                              
```
