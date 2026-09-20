# Read and write Zoo archives

`unzoo` extracts a `.zoo` archive, and `zoo` creates one –
[`unzip`](https://rdrr.io/r/utils/unzip.html) and
[`zip`](https://rdrr.io/r/utils/zip.html) for Zoo rather than Zip. This
is tdaR's code, not part of TDA (TDA can only load an archive's contents
into its data matrix, via a command that needs a companion `.zad`
description file, and has no general extract or create capability at
all).

## Usage

``` r
unzoo(zoofile, files = NULL, list = FALSE, exdir = ".", overwrite = TRUE)

zoo(zoofile, files, method = 1L, names = basename(files))
```

## Arguments

- zoofile:

  path to a `.zoo` file, or an `http(s)://` URL, which is downloaded to
  a temporary file first.

- files:

  paths of the files to add.

- list:

  if `TRUE`, return a data frame describing the archive's contents
  instead of extracting anything.

- exdir:

  directory to extract into, created if it does not exist.

- overwrite:

  whether to replace files already present in `exdir`.

- method:

  packing method for the members: `1` (LZD, the default) or `0` (stored,
  uncompressed). TDA's `arcd` reads all of 0, 1 and 2, but other readers
  of TDA-written archives may not, and LZD is what `zoo` itself has
  always produced.

- names:

  the names the members get in the archive, one per file: by default the
  files' base names. A name may carry a directory, `"sub/x.dat"`, with
  `/` as the separator; it must be relative and may not contain `..`.

## Value

For `list = TRUE`, a data frame with one row per member: `name` (with
its directory, if the entry stores one), `method` (0 stored, 1 LZD, 2
LZH), `size` (original, uncompressed), and `ok` (whether this reads
correctly – always `TRUE` for methods 0/1/2, since those are the only
methods Zoo defines). Otherwise, the paths of the files written,
invisibly the same as [`unzip`](https://rdrr.io/r/utils/unzip.html).

## Details

`unzoo` decodes both of Zoo's packing methods, LZD and LZH, ported from
TDA's reference implementation – the same code path TDA itself has used
since the 1990s to read `.zoo` archives, so archives TDA itself could
read, this reads too. Both methods are tested against a historical TDA
archive (`deha1.zoo` from the TDA teaching pages, kept outside the
package as an external test fixture): every decoded byte matches the
CRC-16 the original archive itself stored, and separately, every file
`zoo` writes here was extracted correctly, with a passing CRC-16, by the
independent, unrelated `zoo`/`unzoo` reference tools
(<https://github.com/troglobit/zoo>) – not just read back correctly by
this same code, which would only prove self-consistency. `zoo` writes
LZD-compressed entries by default, or stored (uncompressed) ones with
`method = 0`, each with its CRC-16 checksum, readable by any Zoo
implementation including TDA's.

Member names are Zoo 2.1's: a long file name and a directory name are
stored in the directory entry beside the DOS-style 8.3 short name, and
both `unzoo` and `zoo` use them, so `sub/long_name.dat` round-trips as
itself. `unzoo` creates the directories it needs under `exdir`. An
archive `zoo` writes is, field for field, what the reference `zoo` tool
writes for the same members (its header, every directory entry's names
and version bytes, the file leaders, the terminal record; only the
members' dates differ, which `zoo` leaves unset), and TDA's own `arcd`
finds a member by its full stored name. TDA's archive description file
separates its fields by blanks, so a member TDA is to read may not have
a blank in its name.

Zoo can, in principle, keep more than one stored version of a file under
the same name (its version-history feature) – `unzoo` does not collapse
or deduplicate by name, so if an archive has several entries called the
same thing, all of them come back as separate rows/files rather than
only the latest.

## Examples

``` r
archive <- system.file("extdata", "tda.zoo", package = "tdaR")
unzoo(archive, list = TRUE)
#>        name method size   ok
#> 1 adata.dat      1  480 TRUE
#> 2  avar.dat      1  120 TRUE

ex <- tempfile()
paths <- unzoo(archive, files = "avar.dat", exdir = ex)
readLines(paths)
#> [1] "V1 1 0 8.0 first variable              "
#> [2] "V2 1 8 8.0 second variable             "
#> [3] "V3 1 16 8.0 third variable             "

# round-trip: write what was just extracted back out, and read it again
out <- tempfile(fileext = ".zoo")
zoo(out, paths)
unzoo(out, list = TRUE)
#>       name method size   ok
#> 1 avar.dat      1  120 TRUE

# names longer than 8.3 and a directory are kept as given
zoo(out, paths, names = "results/avar_first_run.dat")
unzoo(out, list = TRUE)
#>                         name method size   ok
#> 1 results/avar_first_run.dat      1  120 TRUE
ex2 <- tempfile()
unzoo(out, exdir = ex2)
list.files(ex2, recursive = TRUE)
#> [1] "results/avar_first_run.dat"
if (FALSE) { # \dontrun{
unzoo("https://example.com/archive.zoo", exdir = tempfile())
} # }
```
