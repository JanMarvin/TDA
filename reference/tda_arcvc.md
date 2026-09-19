# Check a variable description file

TDA's `arcvc`.

## Usage

``` r
tda_arcvc(file, rewrite = NULL, suffix = c("number", "lower", "upper"), ...)
```

## Arguments

- file:

  the description file.

- rewrite:

  optional path for a rewritten file with unique names.

- suffix:

  how names are made unique: "number", "lower", or "upper" (file number
  or file name in either case).

- ...:

  passed to [`tda_run`](tda_run.md).

  Requires an open archive (see [`tda_arcd`](tda_arcd.md)) containing a
  type-2 (variable description) member.

## Examples

``` r
# the file arcv writes, checked for duplicate variable names
d <- tempfile(); dir.create(d)
for (f in c("tda.zad", "tda.zoo"))
    file.copy(system.file("extdata", f, package = "tdaR"), d)
tda_arcv("vars.vd", archive = "tda.zad", dir = d)
# arcv writes a template with every definition commented out, for the
# caller to uncomment; arcvc reports its verdict on it either way
out <- tda_arcvc("vars.vd", dir = d)
cat(grep("variable|once", out, value = TRUE), sep = "\n")
#> Checking variable description file: vars.vd
#> Found 0 variables.
#> Number of variable names used more than once: 0
```
