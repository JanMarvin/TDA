# Quadratic programming

Minimize a convex quadratic subject to box bounds and/or linear
equality/inequality constraints – TDA's `mqp`, `mqpb`, and `mqpc`
(undocumented in `tda.hlp`; identified from `qld`'s header, the
Powell/Schittkowski active-set solver `ZQPCVX` shared with
[`tda_mlsei1`](https://janmarvin.github.io/TDA/reference/tda_mls.md)).
The problem solved is \$\$\min\\ d'x + \tfrac12 x'Cx \quad
\mathrm{s.t.}\\ A_e x = b_e,\\ A_i x \ge b_i,\\ l \le x \le u\$\$
`tda_mqp` is unconstrained (box bounds fixed at +-huge internally by
TDA); `tda_mqpb` adds box bounds; `tda_mqpc` adds general
equality/inequality constraints instead of bounds – no single TDA
command accepts both bounds and general constraints, confirmed by
reading `m_mqp`'s three branches. Verified against hand-solved instances
for all three commands: an unconstrained case, a box-clipped case, and a
mixed equality+inequality case – see `examples/coverage/mqpops2.cf`.

## Usage

``` r
tda_mqp(C, d, dir = tempfile("tda"), ...)

tda_mqpb(C, d, lower, upper, dir = tempfile("tda"), ...)

tda_mqpc(
  C,
  d,
  equalities = NULL,
  equalities_bounds = NULL,
  inequalities = NULL,
  inequalities_bounds = NULL,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- C:

  the symmetric objective matrix (n x n).

- d:

  the linear objective term, length n.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- lower, upper:

  box bounds, each length n. `tda_mqpb` only.

- equalities:

  an equality constraint matrix `Ae x = equalities_bounds`. `tda_mqpc`
  only.

- equalities_bounds:

  the right-hand sides for `equalities`.

- inequalities:

  an inequality constraint matrix `Ai x >= inequalities_bounds`.
  `tda_mqpc` only.

- inequalities_bounds:

  the right-hand sides for `inequalities`.

## Value

A list: `x`, the solution; `value`, the objective `d'x + 0.5 x'Cx` at
`x`, computed directly in R since this command family prints no value of
its own.

## See also

Other optimization:
[`tda_boolean_min()`](https://janmarvin.github.io/TDA/reference/tda_boolean_min.md),
[`tda_mlp()`](https://janmarvin.github.io/TDA/reference/tda_mlp.md),
[`tda_mlpi()`](https://janmarvin.github.io/TDA/reference/tda_mlpi.md),
[`tda_mls()`](https://janmarvin.github.io/TDA/reference/tda_mls.md),
[`tda_mqap()`](https://janmarvin.github.io/TDA/reference/tda_mqap.md)

## Examples

``` r
# minimize x1^2 + x2^2 - 2x1 - 4x2 (unconstrained optimum (1, 2))
tda_mqp(diag(c(2, 2)), c(-2, -4))
#> $x
#> [1] 1 2
#> 
#> $value
#> [1] -5
#> 
tda_mqpb(diag(c(2, 2)), c(-2, -4), lower = c(0, 0), upper = c(0.5, 0.5))
#> $x
#> [1] 0.5 0.5
#> 
#> $value
#> [1] -2.5
#> 
```
