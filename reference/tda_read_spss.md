# Read an SPSS data file

`tda_read_spss` reads an SPSS system file (`.sav`, TDA's `rspss1`) or a
portable file (`.por`, TDA's `rspss`) and returns its data as a data
frame. The format is taken from the file extension unless `portable`
says otherwise.

## Usage

``` r
tda_read_spss(
  file,
  portable = grepl("\\.por$", file, ignore.case = TRUE),
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- file:

  path to a `.sav` or `.por` file.

- portable:

  read as a portable (`.por`) file. Defaults to whether `file` ends in
  `.por`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame, with `attr(x, "run")` carrying the run.

## Details

These two commands were in TDA all along but had no wrapper, and no
test: the only case in the suite fed `rspss1` a text file and pinned its
refusal. Both are now checked against real files written by readspss,
and both round-trip their values exactly.

String columns *do* come through: `rspss`/`rspss1` create real TDA
string variables (type 1) and `pdata` prints them, so a character column
arrives as text rather than a code. This is one of the few ways to get a
string variable into TDA at all –
[`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md)
refuses a character column in a data frame, because TDA's data-frame
path has no type for it.

Two consequences worth knowing. TDA renames any variable whose name does
not start with a capital (see *Variable names must start with a
capital*), so a lowercase SPSS name comes back as `VAR1`, `VAR2`, ...;
write the file with capitalised names if you need them preserved. And
once a string variable exists, TDA's string operators work on it in a
later `nvar` block of the same run – `strlen(S)`, `strsp(S)` (the
alphabetical sort position, matching R's `rank`), `strv(S)` and
`strvp(S,n,m)`.

For value labels and missing-value codes, use readspss directly: TDA
keeps the data and the names, so a label or a user-missing code is lost
on the way in, and this wrapper does not pretend otherwise.

## See also

Other data:
[`tda_read_xls()`](https://janmarvin.github.io/TDA/reference/tda_read_xls.md),
[`tda_strings()`](https://janmarvin.github.io/TDA/reference/tda_strings.md)

## Examples

``` r
# readspss can write the fixture, so this runs where it is installed
if (requireNamespace("readspss", quietly = TRUE)) {
  f <- file.path(tempdir(), "example.sav")
  readspss::write.sav(data.frame(ID = 1:3, X = c(1.5, 2.5, 3.5),
                                 S = c("a", "bb", "ccc")), f)
  tda_read_spss(f)
}
#>   ID   X   S
#> 1  1 1.5   a
#> 2  2 2.5  bb
#> 3  3 3.5 ccc
```
