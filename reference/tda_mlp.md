# Linear programming

Maximize (or minimize) a linear objective subject to linear inequality
and/or equality constraints – TDA's `mlp` and `mlp1` (Salazar & Sen's
MINIT, CACM algorithm 333). `tda_mlp` takes only inequality constraints,
matching `mlp(T,X,Y)`; `tda_mlp1` adds equality constraints, matching
`mlp1(T,p,X,Y)` – the two are the same underlying call with `p` fixed at
0 for `tda_mlp` (confirmed by reading `m_mlp`'s `opt` branch). The
problem solved, after negating the objective internally for
`direction = "min"`, is \$\$\max\\ c'x \quad \mathrm{s.t.}\\ x \ge 0,\\
Ax \le b,\\ A_e x = b_e\$\$ with the dual \\\min\\ b'y \\\mathrm{s.t.}\\
y \ge 0,\\ A'y \ge c\\ (this is all stated directly in `lpf1`'s header
in `t_lp.c`; nothing here was guessed). Equality-constraint rows report
a dual value of 0 – TDA's convention, not something this wrapper
computes. Verified by hand for a plain and an equality instance, and
against brute-force vertex enumeration in `test-r-comparisons.R`.

## Usage

``` r
tda_mlp(
  objective,
  constraints,
  bounds,
  direction = c("max", "min"),
  dir = tempfile("tda"),
  ...
)

tda_mlp1(
  objective,
  constraints = NULL,
  bounds = NULL,
  equalities = NULL,
  equalities_bounds = NULL,
  direction = c("max", "min"),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- objective:

  coefficients `c`, one per variable (at least one).

- constraints:

  the inequality constraint matrix `A` (`Ax <= bounds`), one row per
  constraint, one column per variable. For `tda_mlp1`, may be omitted if
  `equalities` is given instead.

- bounds:

  the right-hand sides `b` for `constraints`.

- direction:

  `"max"` (default) or `"min"`.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- equalities:

  an equality constraint matrix `Ae` (`Ae x = equalities_bounds`), same
  number of columns as `constraints`. `tda_mlp1` only.

- equalities_bounds:

  the right-hand sides `be` for `equalities`. `tda_mlp1` only.

## Value

A list: `value`, the objective at the optimum; `x`, the primal solution;
`y`, the dual solution (one entry per constraint row, inequalities then
equalities, equality rows always 0).

## Details

`tda_mlp` also accepts TDA's single-tableau form, `mlp(T, X, Y)`: the
first row is the objective `(c, 0)` and the remaining rows are the
constraints `[A, b]`, exactly the layout the solver documents.
`tda_mlp(T)` decomposes it accordingly.

The manual's Box example writes `mlp(T, X, Y, 1)` with a fourth operand.
That form does not exist in this source tree: the shipped 6.4q parser
(checked against a pristine copy of the upstream source) accepts exactly
`mlp(T, X, Y)` and `mlp1(T, p, X, Y)`, so the box was produced by a
different release. On the box's problem the three-operand form gives the
same value (1) and the same dual vector (1, 0); only the primal solution
differs – (1, 0) here against the manual's (1, 1) – and both are
vertices of the same optimal face, since maximising x1 under x1 \<= 1,
x2 \<= 1 leaves x2 free at the optimum. Nothing disagrees beyond the
tie-break.

## See also

Other optimization:
[`tda_boolean_min()`](https://janmarvin.github.io/TDA/reference/tda_boolean_min.md),
[`tda_mlpi()`](https://janmarvin.github.io/TDA/reference/tda_mlpi.md),
[`tda_mls()`](https://janmarvin.github.io/TDA/reference/tda_mls.md),
[`tda_mqap()`](https://janmarvin.github.io/TDA/reference/tda_mqap.md),
[`tda_mqp()`](https://janmarvin.github.io/TDA/reference/tda_mqp.md)

## Examples

``` r
# maximize 3x1 + 2x2 s.t. x1 + x2 <= 4, x1 + 3x2 <= 6
tda_mlp(c(3, 2), rbind(c(1, 1), c(1, 3)), c(4, 6))
#> $value
#> [1] 12
#> 
#> $x
#> [1] 4 0
#> 
#> $y
#> [1] 3 0
#> 
# same, plus x1 + x2 = 3 exactly
tda_mlp1(c(3, 2), rbind(c(1, 3)), c(6), equalities = rbind(c(1, 1)),
         equalities_bounds = 3)
#> $value
#> [1] 9
#> 
#> $x
#> [1] 3 0
#> 
#> $y
#> [1] 0 0
#> 
```
