# Build TDA's input

`tda_write_data` writes a data frame as the whitespace-delimited file
TDA reads; `tda_nvar` generates the variable block that describes it;
`tda_block` builds any other command. The typed interfaces use these,
and they are exported so a hand-written script can too.

## Usage

``` r
tda_write_data(data, file, na = ".")

tda_nvar(
  data,
  file = NA_character_,
  vars = colnames(data),
  mpnt = -5,
  fmt = NULL,
  extra = NULL,
  ...
)

tda_block(name, ..., rhs = NULL)
```

## Arguments

- data:

  a data frame.

- file:

  output file name.

- na:

  the marker written for a missing value. It has to be a marker rather
  than a number: a number is data to TDA, counted and used by the fits.
  "." is the point marker `nvar` understands.

- vars:

  variable names, defaulting to the column names.

- mpnt:

  the value the point marker stands for. TDA substitutes it and reports
  how many there were.

- fmt:

  print formats, one per variable.

- extra:

  further variable definitions, as text.

- ...:

  further options, as `name = value`.

- name:

  a command name.

- rhs:

  the command's right-hand side.

## Value

A character vector of commands, or the file name, invisibly.

## Numeric precision, storage, and the exports

Every wrapper hands its data over through `rdataframe`, which stores R
double columns as TDA `<8>` double arrays and integer columns as `<5>`
int arrays – nothing is squeezed through a narrower type on the way in.
All of TDA's estimation arithmetic is double, and the direct-export
channel (`$run$exports`) copies doubles into R doubles, so an exported
value is bit-for-bit the number TDA computed; only the *printed* output
rounds, at `tfmt`/`mfmt`. Two deliberate exceptions to know about: TDA's
default storage for variables declared in a hand-written command file
(`tda_run(cf=)`, or `tda_nvar` with `file=` and `fmt=`) is a 4-byte
float, faithful to the original program – declare `<8>` in the
definition if you need double storage there; and a few of TDA's internal
summaries (the episode table's weighted counts, for one) pass through
float by TDA's design, which the package reproduces rather than
corrects.

## The `tdaR.use_exports` option

Every reader takes its numbers from the export channel and falls back to
parsing TDA's printed output or its written files only if the export is
absent or does not fit. Setting `options(tdaR.use_exports = FALSE)`
forces that fallback everywhere, which is useful for two things:
comparing the two paths when a value looks wrong, and reproducing what
an older version of the package returned.

Expect the parser path to give you *less precision*, not different
answers – it reads values TDA has already rounded for printing,
typically to four or five decimals. A few readers also differ in shape,
because a text table cannot always express what the export can:
`tda_atab` keeps the open first and last classes (the cases outside your
breaks) on the export path and drops them on the parser path, since TDA
prints those rows with a blank bound that the fixed-width reader cannot
parse.

The package's test suite runs under both settings, so the fallback stays
working rather than bit-rotting.

## Variable names must start with a capital

TDA's parser accepts a variable name only if it starts with an uppercase
letter, `_`, `@` or `$`; a lowercase start is read as some other token.
This is not an arbitrary restriction that could be patched out: the
first character's case is how TDA's expression lexer tells variables
apart from everything spelled in lowercase – functions (`log`, `exp`,
`if`), comparison operators (`eq`, `le`, `gt`), the episode accessors
(`ts`, `tf`, `org`, `des`, `time`), random draws (`rd`, `rdn`), and the
estimated parameters (`b1`, `b2`, ...). Allowing a lowercase variable
would make a column called `time` or `b1` ambiguous against those, with
no rule to break the tie, so the case requirement stays and the R side
translates instead.

The typed wrappers do that translation invisibly: a lowercase name is
sent to TDA as `V<name>` (`.tda_names()`, internal) and the caller's
spelling is put back on everything returned – coefficient tables,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) dimnames, table columns,
labels. The `V` spelling appears only in raw TDA console output
(`$run$output`) and in command files, where it is TDA's name for the
column. Only when writing `tda_nvar`/`tda_block` commands by hand does
the requirement reach you: use names that start with a capital there, or
pass them through the same convention.

The symptom varies by exactly where the bad name ends up, and none of
the variants mention case, which makes this easy to lose an hour to (it
has happened more than once writing this package): a lowercase name used
to *reference* an existing variable fails with “Syntax error or
undefined variables”, pointing at the command, not the name; a lowercase
name used to *declare* a new one – `x = c4` rather than `X = c4` in a
hand-built `sdnvar` or `nvar` block, say – fails differently and more
confusingly, with “Syntax error: x=c4,)”, which looks like a problem
with the column reference, the trailing comma, or a missing type/format
specifier, and is none of those.

## String variables and the string operators

A character or factor column becomes a TDA *string variable*: `tda_nvar`
declares it with `= str(n,m)` (TDA's only way of making one) and
`tda_write_data` pads every field to a fixed width so those character
positions are known. Strings are written last in the file whatever order
your columns are in, because TDA cannot reach a `= cK` field lying after
a string one; the declaration keeps your order, so nothing changes on
the R side.
[`tda_strings`](https://janmarvin.github.io/TDA/reference/tda_strings.md)
returns them as character.

Four operators work on a string variable, and `extra=` accepts a
[`{ }`](https://rdrr.io/r/base/Paren.html) block of plain R, so they
read as R rather than TDA syntax:


      tda_nvar(d, extra = {
          Len  = strlen(S)        # storage width, not visible characters
          Rank = strsp(S)         # alphabetical position, as R's rank()
          Num  = strv(Code)       # digits to a number, -1 if not all digits
          Sub  = strvp(Code, 2, 3)  # characters 2..3 as a number
      })

`strlen` returns the declared storage size, so a padded `"fig"` in a
width-4 variable is 4, not 3. `strv` and `strvp` return `-1` when the
text is not all digits, rather than failing.

## See also

Other TDA infrastructure:
[`tda_ccnt()`](https://janmarvin.github.io/TDA/reference/tda_ccnt.md),
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md),
[`tda_output()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md),
[`tda_read_table()`](https://janmarvin.github.io/TDA/reference/tda_read_table.md),
[`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md)

## Examples

``` r
d <- data.frame(Vx = c(1, 2, 3, 4), Vy = c(2, 4, 5, 9))

# the plain data file tda_run() writes and TDA reads -- fixed-width
# columns, no header, exactly what dfile= in an nvar block points at
f <- tempfile()
tda_write_data(d, f)
cat(readLines(f), sep = "\n")
#> 1 2
#> 2 4
#> 3 5
#> 4 9

cmds <- c(tda_nvar(d), tda_block("dstat", rhs = "Vx,Vy"))
cat(cmds, sep = "\n")
#> rdataframe;
#> dstat(
#> ) = Vx,Vy;
res <- tda_run(cmds, data = d)
cat(res$output, sep = "\n")
#> TDA. Analysis of Transition Data (6.4q). Sat Sep 19 18:53:25 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: commands
#> ============================================================================
#> Idx Variable  T   S  PFmt  Definition
#> -------------------------------------
#>   1 Vx        3   8  24.16 c1
#>   2 Vy        3   8  24.16 c2
#> 
#> Reading a data frame to create internal data matrix.
#> Maximum number of cases: 4
#> Created a data matrix with 2 variables and 4 cases.
#> ----------------------------------------------------------------------------
#> dstat()=Vx,Vy
#> Descriptive statistics. Current memory: 390132 bytes.
#> 
#> Variable    Minimum    Maximum       Mean   Std.Dev.    Sum of values
#> ---------------------------------------------------------------------
#> Vx           1.0000     4.0000     2.5000     1.2910          10.0000
#> Vy           2.0000     9.0000     5.0000     2.9439          20.0000
#> ----------------------------------------------------------------------------
#> Current memory: 390032 bytes. Max memory used: 390176 bytes.
#> End of program. Sat Sep 19 18:53:25 2026
```
