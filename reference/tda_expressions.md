# Expand an expression, or parse a matrix expression

`tda_expand` runs TDA's `expm`, which expands an expression into its
fully written-out form – a model formula with interaction shorthand
becomes the explicit list of terms. `tda_mparse` runs `mparse`, which
parses a matrix expression and reports how it read it, without
evaluating anything.

## Usage

``` r
tda_expand(expr, setup = NULL, options = list(), dir = tempfile("tda"), ...)

tda_mparse(
  expr,
  setup = NULL,
  mfmt = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- expr:

  the expression, as a single string.

- setup:

  TDA commands to run first, as a character vector. A matrix expression
  needs its matrices to exist, so this is where `mdef`/`mdefi` lines go.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- mfmt:

  for `tda_mparse`: print format, TDA's default is `12.4`.

## Value

A character vector: the lines TDA printed in response. `attr(x, "run")`
carries the run.

## Details

Both are diagnostics: they answer "what did TDA make of what I wrote",
which is otherwise only visible in an error message.

## Examples

``` r
tda_expand("A*B")
#> [1] "expm=A*B"
#> attr(,"run")
#> TDA. Analysis of Transition Data (6.4q). Sun Sep 20 09:48:14 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: commands
#> ============================================================================
#> expm=...
#> Expand model description. Current memory: 390032 bytes.
#> 
#> expm=A*B
#> ----------------------------------------------------------------------------
#> Current memory: 390032 bytes. Max memory used: 390036 bytes.
#> End of program. Sun Sep 20 09:48:14 2026
tda_mparse("A+B", setup = c("mdefi(2,2,A);", "mdefi(2,2,B);"))
#> [1] " Cnt        Typ        Val   Dimension"                   
#> [2] "--------------------------------------"                   
#> [3] "   0       8000     0.0000      2 x 2"                    
#> [4] "   1       8001     1.0000      2 x 2"                    
#> [5] "   2          1     2.0000  "                             
#> [6] "Dimension: 2 x 2. Max stack length: 2. Max matrix size: 4"
#> [7] "      2.0000       0.0000 "                               
#> [8] "      0.0000       2.0000 "                               
#> attr(,"run")
#> TDA. Analysis of Transition Data (6.4q). Sun Sep 20 09:48:14 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: commands
#> ============================================================================
#> mparse()=A+B
#> 
#>  Cnt        Typ        Val   Dimension
#> --------------------------------------
#>    0       8000     0.0000      2 x 2
#>    1       8001     1.0000      2 x 2
#>    2          1     2.0000  
#> 
#> Dimension: 2 x 2. Max stack length: 2. Max matrix size: 4
#> 
#>       2.0000       0.0000 
#>       0.0000       2.0000 
#> 
#> ----------------------------------------------------------------------------
#> Current memory: 390032 bytes. Max memory used: 402220 bytes.
#> End of program. Sun Sep 20 09:48:14 2026
```
