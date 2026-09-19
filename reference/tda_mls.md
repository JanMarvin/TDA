# Constrained least squares

Least squares with optional linear equality constraints, linear
inequality constraints, and non-negativity – TDA's `mls`, `mlse`,
`mlsi`, `mlsei`, `mlsei1`, and `mnls`, each as its thin wrapper sharing
one implementation. The problem solved is \$\$\min\\ \\Ax - b\\\_2 \quad
\mathrm{s.t.}\\ A_e x = b_e,\\ A_i x \ge b_i,\\ (x \ge 0)\$\$ These are
TDA's documented semantics (`tda.hlp`, entries
`##mls`/`##mlse`/`##mlsi`/`##mlsei`/`##mnls`): an earlier working note
had assumed `mlse`/`mlsi` took their constraint counts, which the source
does not support – `tda_mlse` and `tda_mls` turn out to make the
identical call (confirmed empirically: an overdetermined, inconsistent
instance returns the plain least-squares fit rather than an error), and
`tda_mlsi` treats its entire input as inequality rows. `tda_mlsei` and
`tda_mlsei1` solve the same general equality+least-squares+inequality
problem via two different TDA algorithms (Lawson-Hanson vs. an
undocumented active-set QP, confirmed to agree to solver tolerance on
two independent instances, not merely assumed from the similar name).
`tda_mnls` is the only one of the six with non-negativity, applied to
every variable (TDA's `k` split point fixed at 0) and cannot be combined
with `inequalities` (`mnls` has no general inequality block of its own).

## Usage

``` r
tda_mls(A, b, dir = tempfile("tda"), ...)

tda_mlse(A, b, dir = tempfile("tda"), ...)

tda_mlsi(inequalities, inequalities_bounds, dir = tempfile("tda"), ...)

tda_mlsei(
  A = NULL,
  b = NULL,
  equalities = NULL,
  equalities_bounds = NULL,
  inequalities = NULL,
  inequalities_bounds = NULL,
  dir = tempfile("tda"),
  ...
)

tda_mlsei1(
  A = NULL,
  b = NULL,
  equalities = NULL,
  equalities_bounds = NULL,
  inequalities = NULL,
  inequalities_bounds = NULL,
  dir = tempfile("tda"),
  ...
)

tda_mnls(
  A,
  b,
  equalities = NULL,
  equalities_bounds = NULL,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- A:

  the least-squares data matrix (`ma` x n). Required for
  `tda_mls`/`tda_mlse`/`tda_mnls`; optional for `tda_mlsei`/`tda_mlsei1`
  (`NULL` for a pure feasibility problem, together with `equalities`
  and/or `inequalities`); not used by `tda_mlsi`.

- b:

  the target vector for `A x ~ b`, length `ma`.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

- inequalities:

  an inequality constraint matrix `Ai x >= inequalities_bounds`.
  `tda_mlsi` (required, its only input) and `tda_mlsei`/`tda_mlsei1`
  (optional).

- inequalities_bounds:

  the right-hand sides for `inequalities`.

- equalities:

  an equality constraint matrix `Ae x = equalities_bounds`. `tda_mlsei`/
  `tda_mlsei1`/`tda_mnls` only.

- equalities_bounds:

  the right-hand sides for `equalities`.

## Value

A list: `x`, the solution; `residual`, the Euclidean norm of
`A %*% x - b` (`NA` when `A` is absent), computed directly in R from
`A`, `x` and `b` rather than trusted from TDA's text, since `tda_mlsei1`
prints no residual at all; `rank`, the reported rank of `A` (`NA` for
`tda_mnls` and `tda_mlsei1`, neither of which report one).

## Details

All three also accept the manual's single-matrix form, `mls(S, B)` with
`S = [X, y]`: when the right-hand side is not given, the last column of
the first argument is split off as `y`.

## See also

Other optimization: [`tda_boolean_min()`](tda_boolean_min.md),
[`tda_mlp()`](tda_mlp.md), [`tda_mlpi()`](tda_mlpi.md),
[`tda_mqap()`](tda_mqap.md), [`tda_mqp()`](tda_mqp.md)

## Examples

``` r
# unconstrained: fit x1, x2 to two noisy targets
tda_mls(diag(2), c(2.1, 2.9))
#> $x
#> [1] 2.1 2.9
#> 
#> $residual
#> [1] 0
#> 
#> $rank
#> [1] 2
#> 
# the manual's augmented form S = [X, y] works too
tda_mls(cbind(diag(2), c(2.1, 2.9)))
#> $x
#> [1] 2.1 2.9
#> 
#> $residual
#> [1] 0
#> 
#> $rank
#> [1] 2
#> 
# x1 + x2 = 3, x1 close to 2, x2 close to 3
tda_mlsei(diag(2), c(2, 3), equalities = rbind(c(1, 1)),
         equalities_bounds = 3)
#> $x
#> [1] 1 2
#> 
#> $residual
#> [1] 1.414214
#> 
#> $rank
#> [1] 1
#> 
```
