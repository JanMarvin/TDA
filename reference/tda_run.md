# Run TDA commands

The escape hatch: anything the typed interfaces do not cover can be
written as TDA commands and run directly. `tda_run` writes a data frame
out and runs a script against it; `tda_run_cf` runs an existing command
file.

## Usage

``` r
tda_run(
  commands,
  data = NULL,
  dir = tempfile("tda"),
  data_file = "data.dat",
  ...
)

tda_run_cf(
  file,
  maxnv = 2000L,
  maxmat = 200L,
  args = character(),
  echo = FALSE,
  data = NULL
)
```

## Arguments

- commands:

  a character vector of TDA commands.

- data:

  optional data frame, written to `data_file` in `dir` so that the
  commands can read it with `dfile=`.

- dir:

  working directory; TDA's output files are left there.

- data_file:

  name for the written data file.

- ...:

  passed to `tda_run_cf`.

- file:

  an existing command file.

- maxnv, maxmat:

  limits passed to TDA at start-up.

- args:

  extra command line arguments.

- echo:

  print TDA's output as it is read.

## Value

An object of class `tda_result`: `$output` and `$stderr` (the captured
text), `$dir`, `$cf` and `$files` (the working directory, the command
file, and the files the run produced), `$commands` (the command file's
lines, kept on the object), `$diagnostics` (every error, warning or note
line TDA printed, see `tda_diagnostics`), `$exports` (TDA's numbers,
handed over as matrices in parallel with the printed output), `$errors`
(how many errors TDA reported) and `$status` (TDA's exit code, which is
0 for a run that reported errors as well – test `$errors`, not this).

## See also

Other TDA infrastructure:
[`tda_ccnt()`](https://janmarvin.github.io/TDA/reference/tda_ccnt.md),
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md),
[`tda_output()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md),
[`tda_read_table()`](https://janmarvin.github.io/TDA/reference/tda_read_table.md),
[`tda_write_data()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)

## Examples

``` r
tda_run("mem;")
#> TDA. Analysis of Transition Data (6.4q). Sun Sep 20 03:44:58 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: commands
#> ============================================================================
#> Currently requested memory: 390032 (390032) bytes.
#> Current memory: 390032 bytes. Max memory used: 390032 bytes.
#> End of program. Sun Sep 20 03:44:58 2026

# tda_run_cf: the same thing, but the command file already exists on disk
f <- tempfile()
writeLines("mem;", f)
r <- tda_run_cf(f)
cat(r$output, sep = "\n")
#> TDA. Analysis of Transition Data (6.4q). Sun Sep 20 03:44:58 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: file1fe472ad083a
#> ============================================================================
#> Currently requested memory: 390032 (390032) bytes.
#> Current memory: 390032 bytes. Max memory used: 390032 bytes.
#> End of program. Sun Sep 20 03:44:58 2026
```
