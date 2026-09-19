# Matrix operations

Runs any of TDA's matrix commands on R matrices and returns R matrices.
This is a generic interface rather than one function per operation: R
already has `%*%`, `solve` and `eigen`, so named wrappers would add
nothing, but the whole of TDA's matrix language stays reachable.

## Usage

``` r
tda_mat(op, ..., out = "R", options = list(), dir = tempfile("tda"))
```

## Arguments

- op:

  the command name. TDA's matrix language has about eighty: `mmul`
  multiplies, `minvs` and `mginv` invert, `mtransp` transposes, `mev`
  and `mevs` give eigenvalues, `msvd` the singular value decomposition,
  `mchol` a Cholesky factor, `mrank` the rank. Note the names: inversion
  is `minvs`, not `minv`, and transposition `mtransp`, not `mtra`.

- ...:

  input matrices, in the order the command expects.

- out:

  names of the result matrices to read back.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

A matrix if one result was asked for, otherwise a named list of them.

## How it works

The inputs are written as TDA matrices named `A`, `B`, ... in order, the
command is run, and every matrix named in `out` is read back. TDA's
convention is that results come last in the argument list, so
`mmul(A,B,R)` multiplies `A` by `B` into `R`.


    tda_mat("mmul", A, B, out = "R")            # A 
    tda_mat("mev", A, out = c("ER", "EI"))      # eigenvalues, real and imaginary

## See also

Other smoothing:
[`tda_integrate()`](https://janmarvin.github.io/TDA/reference/tda_integrate.md),
[`tda_interp()`](https://janmarvin.github.io/TDA/reference/tda_interp.md),
[`tda_isotonic()`](https://janmarvin.github.io/TDA/reference/tda_isotonic.md),
[`tda_sma()`](https://janmarvin.github.io/TDA/reference/tda_sma.md),
[`tda_smd()`](https://janmarvin.github.io/TDA/reference/tda_smd.md),
[`tda_spl()`](https://janmarvin.github.io/TDA/reference/tda_spl.md)

## Examples

``` r
A <- matrix(1:4, 2)
B <- diag(2)
tda_mat("mmul", A, B, out = "R")
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
```
