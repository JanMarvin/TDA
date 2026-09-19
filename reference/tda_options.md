# TDA options not exposed as named arguments

Most wrappers accept an `options = list(...)` argument that is passed
straight through to the TDA command they build. TDA has several hundred
options and promoting each to a named R argument would be unreadable, so
the named arguments cover what is needed in ordinary use and `options=`
reaches the rest.

## Which options a command accepts

TDA's manual is installed with the package and is the authoritative
list. [`tda_help`](tda_help.md) prints the section for any command,
including its full option table:


      tda_help("rate")      # every option the rate command takes
      tda_help("sdgshhs")   # the GSHHS reader's options

The name in `options=` is TDA's name, exactly as that table spells it,
so `options = list(grp = "SEX")` becomes `grp = SEX,` in the generated
command.

## options= and ... are not the same thing

Where a wrapper has both, `options=` goes to the *TDA command* and `...`
goes to [`tda_run`](tda_run.md) (working directory, echo, and so on).
Passing a TDA option in `...` is therefore an error rather than a silent
no-op:


      tda_rate(..., grp = "SEX")                   # unused argument
      tda_rate(..., options = list(grp = "SEX"))   # reaches the command

A few readers – [`tda_read_gshhs`](tda_read_spatial.md) and the other
spatial readers among them – have no separate `options=` and forward
`...` into the command directly, so for those
`tda_read_gshhs(f, opt = 2)` is right.

## Caveat

An option reaching the command is not the same as the option working.
TDA validates its parameters and some are only meaningful for particular
models: `grp=` on `rate`, for instance, reaches the command and is
accepted, but the fit then returns nothing unless the model is one of
the stratified ones. Check the run's output when an option seems to have
no effect.

## Examples

``` r
# what does the GSHHS reader accept?
tda_help("sdgshhs")
#> sdgshhs: This command can be used to convert a GSHHS binary file into a TDA
#> spatial data file. The syntax is:
#> 
#>     sdgshhs(
#>         df=...,     output file (required)
#>         opt=...,    option, def. 1
#>                     1 = standard coding of longitudes
#>                     2 = longitudes unchanged
#>         level=...,  select polygons with specified level, def. 0
#>         fmt=...,    floating point print format, def. 10.5
#>         dtda=...,   create TDA description file
#> 
#>     ) = gshhs_input_file;
#> 
#> The input file must be a binary file from the GSHHS data base. For more
#> information see Section 9.4.3.2 of the TDA Manual.

# opt = 2 keeps longitudes on 0..360, which stops dateline-crossing
# polygons being drawn straight across a world map
if (FALSE) { # \dontrun{
tda_read_gshhs("gshhs_i.b", opt = 2, level = 1)
} # }
```
