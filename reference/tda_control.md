# Convergence control

Convergence control

## Usage

``` r
tda_control(
  algorithm = NULL,
  tolsg = NULL,
  tolsp = NULL,
  tolg = NULL,
  tolf = NULL,
  maxit = NULL,
  tol_param = NULL,
  tol_var = NULL,
  tol_reduction = NULL,
  step_length = NULL,
  step_reduction = NULL,
  step_min = NULL,
  maxit_line = NULL,
  criterion = NULL,
  derivatives = NULL
)

tda_strict(tol = 1e-12)

tda_converged(x)
```

## Arguments

- algorithm:

  TDA minimiser number, its `mina=`: 1 direct search, 2 simplex, 3
  conjugate gradients, 4 BFGS, 5 Newton (I, TDA's own default), 6 Newton
  (II), 7 CES with a quadratic model, 8 CES with a tensor model. 7 and 8
  use the scale-invariant criterion.

- tolsg, tolsp:

  tolerances for the scaled gradient and scaled parameter change, used
  by algorithms 7 and 8.

- tolg, tolf:

  tolerances for the gradient norm and the function value, used by the
  default algorithm.

- maxit:

  iteration cap, TDA's `mxit=`.

- tol_param, tol_var, tol_reduction:

  three more tolerances TDA's minimiser reads for other
  algorithms/situations – `tolp=` (parameter change, direct search and
  others), `tolv=` (Simplex), `tols=` (reduction factor, direct search).

- step_length, step_reduction, step_min:

  step-size controls for TDA's line search – `slen=`, `sred=`, `smin=`.

- maxit_line:

  cap on line-search iterations within one step – TDA's `mxitl=`,
  default 50.

- criterion:

  which convergence test to apply – TDA's `crit=`, 1 to 3. Confirmed
  limited: it cannot select the scale-invariant Dennis-Schnabel test
  described under Criteria above, which is reachable only through
  `algorithm = 7` or `8` directly, not through this argument.

- derivatives:

  `0` for analytical derivatives (TDA's default when unset), `1` or `2`
  for numerical approximations of varying cost – TDA's `dopt=`. Applies
  only to `algorithm = 7` or `8`; harmless but inert for every other
  algorithm.

- tol:

  tolerance for `tda_strict`, applied to both `tolsg` and `tolsp`.

- x:

  a fitted model.

## Value

A list of options for `control=`; `tda_converged` returns `TRUE`,
`FALSE`, or `NA` where TDA reports nothing.

## Criteria

TDA's default minimiser stops on the norm of the gradient against an
absolute tolerance, plus a loose disjunction that also stops when the
change in function value is small. Algorithms 7 and 8 switch to
`optstp()`, the Dennis-Schnabel test on the largest *relative* gradient
and relative parameter step. That is scale-invariant, closer to what R's
optimisers do, and gives results that reproduce across platforms. It is
reachable only through `algorithm`: TDA's `crit=` accepts 1 to 3 and
cannot select it.

The default is TDA's, so a plain call reproduces a plain command file.
`tda_strict()` opts in.

## A caution

Algorithms 7 and 8 also lower the iteration cap to 20 when derivatives
are numerical, against 100 otherwise, so asking for a tighter tolerance
can turn a converged fit into a truncated one. Every fit reports how it
stopped, warns if it did not converge, and `tda_converged()` reads that
back.

## See also

Other rate models:
[`tda_constrain()`](https://janmarvin.github.io/TDA/reference/tda_constrain.md),
[`tda_dple()`](https://janmarvin.github.io/TDA/reference/tda_ltb.md),
[`tda_rates()`](https://janmarvin.github.io/TDA/reference/tda_rates.md),
[`tda_split()`](https://janmarvin.github.io/TDA/reference/tda_split.md),
[`tda_survivor()`](https://janmarvin.github.io/TDA/reference/tda_survivor.md),
[`tda_transitions()`](https://janmarvin.github.io/TDA/reference/tda_transitions.md),
[`vcov.tda_fit()`](https://janmarvin.github.io/TDA/reference/tda_rate.md)

## Examples

``` r
set.seed(36)
d <- data.frame(x = rnorm(200))
d$y <- rbinom(200, 1, plogis(0.3 + 0.8 * d$x))
f <- tda_qreg(y ~ x, d, control = tda_control(maxit = 50, tolg = 1e-8))
coef(f)
#> Intercept         x 
#> 0.3851079 0.9164042 
tda_converged(f)
#> [1] TRUE

# tda_strict: the same idea, one shared tolerance for both algorithm 7/8
# criteria at once
f2 <- tda_qreg(y ~ x, d, control = tda_strict())
coef(f2)
#> Intercept         x 
#> 0.3851079 0.9164042 
```
