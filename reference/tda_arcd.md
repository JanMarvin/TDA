# Open or close a TDA data archive

TDA's `arcd`: opens access described by an archive description file, or
closes it with `file = NULL`.

## Usage

``` r
tda_arcd(file = NULL, ...)
```

## Arguments

- file:

  archive description file, or NULL to close.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the printed output, invisibly.

## Examples

``` r
# an archive description file names the zoo archive and then one
# line per member: number, name, type (1 data, 2 variable
# descriptions), record length, records, variables
dr <- tempfile("tda"); dir.create(dr)
zoo <- system.file("extdata", "tda.zoo", package = "tdaR")
file.copy(zoo, file.path(dr, "tda.zoo"))
#> [1] TRUE
writeLines(c("tda.zoo", "1 adata.dat 1 24 20 3", "2 avar.dat 2 40 3 0"),
           file.path(dr, "arc.ad"))
out <- tda_arcd("arc.ad", dir = dr)
cat(grep("archive|Checking|adata", out, value = TRUE), sep = "\n")
#> Reading archive description file: arc.ad
#> ZOO data archive: tda.zoo
#> Checking definition of files in archive.
#>   1   1     24       20     3       480  1  adata.dat
```
