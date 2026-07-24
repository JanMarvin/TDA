# Changes to the C

TDA 6.4p is Rohwer and Pötter's last release; the version banner reads
6.4q so a patched build is distinguishable from it. Every entry below is
a change to the C itself, in both `src/` and `tdaR/src/`, made only
where TDA's own output was demonstrably wrong -- verified against
independent computation, TDA's manual, or its own references. Each entry
names its control file. The changes are our responsibility, never
attributed to TDA.

- `prn1_coeff()` (`t_pgen.c`): a standard error was printed as `---`
  whenever the coefficient's variance fell below `EPSI1`, which
  suppressed real, legitimately tiny standard errors for large-scaled
  predictors. The threshold is now `DBL_EPSILON`; verified against
  `lm()`/`glm()`, with the constrained-parameter case (where `---` is
  genuinely correct) confirmed unchanged and pinned by dedicated
  regression files.
- `cont_draw()` (`t_cplot.c`): contour fill fragments were written in
  grid-tracing discovery order, so a band could be painted underneath a
  geometrically larger fragment and disappear from the finished plot.
  Fragments are now buffered and flushed sorted by their own measured
  (shoelace) area, largest first.
- `qr_resida()` (`t_qrmod.c`): the block that finalises a covariate's
  standard deviation sat *inside* the loop that accumulates it, so it
  ran after every case -- subtracting the running mean, taking a square
  root, and letting the next iteration add to an already-rooted value.
  Every `Std.Dev.` collapsed to 0 and every `Exp(C*SD)` to 1, for all
  data, making `qreg`'s standardized-coefficient table non-functional.
  Verified against independent computation (the weighted SD of
  `Log10Dose` in `examples/exam/qr1.dat` is `0.27370532608834974`, which
  the fixed build reproduces to 17 digits) and against three further
  datasets that all showed the same collapse. The Intercept row still
  reports 0, correctly: a constant has no standard deviation.
  `examples/exam/qr1a.ref` is the one reference that had this baked in
  and has been regenerated; it is therefore the one block the pinning
  suite no longer checks against Rohwer's 6.4p binary.
- `rmod()` (`t_qrmod.c`): an unconditional `goto RMFin;` sat directly
  after the conditional one that handles a failed minimisation, so
  everything below it was unreachable -- the result printing, the
  additional options, and the `err = 0` that marks success. `rmod`
  estimated a model and then discarded it, reporting failure. The stray
  goto is removed. Found by sweeping the ported sources against Rohwer's
  for lost control flow; it had been fixed in the code but recorded only
  in a comment there.
- `edat_cmp()` (`t_edat.c`): the comparator `qsort` uses to order
  episodes returned 0 for two fully tied entries, which leaves their
  order to the library and so to the platform. It now breaks the tie on
  the original case index, making the sort deterministic. Same entry
  point as above: a real change that lived only in a source comment.
- `sdgshhs()` (`t_sdx.c`): GSHHS has two incompatible header layouts and
  TDA only understood the older one, so any file from GSHHG 2.x (2009
  onwards) stopped with `Error: while reading the file` -- the first
  polygon's points were read from eight bytes inside its own header. The
  header is now read word by word: both layouts open with the same eight
  ints, and the third decides -- in v1 it is `level`, which is 1..4; in
  v2 it is a packed flag whose version byte puts it at 256 or above.
  v2's three trailing ints (`area_full`, `container`, `ancestor`) are
  consumed, and `level`, `greenwich` and `source` are unpacked from the
  flag. The detected version is now reported in the output. Verified on
  a current `gshhs_l.b` (format version 15, 10717 polygons): polygon 1
  comes back with exactly the 6851 points its header declares and a
  latitude range matching the header's own `south`/`north` words,
  decoded independently.

  The point buffers also grow on demand now. `sdgshhs` allocated a fixed
  `GSHHSNP` (1435084) points -- whatever the largest polygon was in the
  release this was built against -- then wrote `AcXF[i]` for `i < h.n`
  with no check. GSHHG has grown since, so a full-resolution file can
  carry a polygon past that limit; AddressSanitizer confirms the
  unguarded code writes off the end of the allocation (`t_sdx.c:1729`,
  buffer from `t_sdx.c:1618`), which in this package means the same dead
  R session as the `rxls` bug above. A low-resolution file cannot reach
  it -- its largest polygon is 6851 points -- so this was found by
  constructing an oversized polygon, not by reading a real file.

  The v1 path is pinned with a synthesised file, since a real one is not
  ours to ship.
- `rxls()` (`t_xls.c`): every cell write went into `AcN`/`AcX` by
  ABSOLUTE row and column, but those arrays are sized from the row and
  column maxima the *first* pass measured. A record naming a cell beyond
  those maxima wrote outside the allocation. readxl's `type-me.xls` does
  exactly that: heap corruption and `munmap_chunk(): invalid pointer`,
  abort. That matters more here than in the standalone program, because
  the R package runs TDA **in-process** -- the abort took the whole R
  session down, with no condition to catch and nothing for `try()` to
  see. Every write now goes through `xls_cellok()`, which range-checks
  the cell and counts what it skips; the count is reported (`Cells
  outside the measured table (skipped): N`). `type-me.xls` reads instead
  of aborting, and the workbooks that already worked are unchanged --
  `datasets.xls` sheet 1 still matches R's own `mtcars` exactly on every
  numeric column.
- `sdgshhs()` (`t_sdx.c`): `opt=2`, which TDA's own help calls
  "longitudes unchanged", left one of two wraps running. The first fires
  whenever a polygon carries the greenwich flag and a point lies east of
  `max_east` -- and `max_east` drops from 270 to 180 degrees after the
  first polygon, so from the second polygon on, every point past 180 was
  moved to the far side regardless of `opt`. A polygon spanning the
  dateline then has consecutive points 350 degrees apart and is drawn as
  a straight line across the whole map. Both wraps are gated by `opt`
  now, so `opt=2` means what it says.

  With `opt=1` a spanning polygon is also **split** at the dateline
  rather than drawn across the map. Its dateline-free runs are emitted
  as separate objects, so each is drawn where it belongs and nothing is
  drawn between them. A split piece is no longer a closed area, so it
  goes out as a line (type 2) rather than a polygon; a run of a single
  point goes out as a point object rather than being dropped, so no
  vertex is lost. A polygon that does not cross is untouched and still
  produces one type 3 object.
- `sdgshhs()` (`t_sdx.c`): `level=` rejected anything above 4, silently
  resetting it to 0 -- "all levels". GSHHG 2.3.0 added level 5 (the
  Antarctic ice front) and level 6 (the grounding line), so asking for
  either drew the entire file instead. The bound is 6 now, and levels
  are additionally selectable as a SET, encoded as a bitmask offset by
  64; a plain `level=1` still means exactly that, so nothing that worked
  before changes. Reachable as `tda_read_gshhs(f, level = c(1, 2))`.
- `ucl()`/`udcl()`/`udcl1()` (`t_cl.c`): the command did not work, and
  the repairs are **ours, not Rohwer's**. It splits an ordering of a
  graph's nodes into contiguous clusters so the largest cluster diameter
  is smallest, by the dynamic program of Alpert and Kahng (*Splitting an
  Ordering into a Partition to Minimize Diameter*, J. of Classification
  14, 1997), which `udcl`'s own comment cites. Five separate faults:

  - the ordering was built as the identity and then had its first two
    entries swapped, hard-coded, so it optimised over an ordering nobody
    asked for;
  - `nc`, and the minimum and maximum cluster size, were fixed at 3, 1
    and 3 whatever the command said;
  - the search over split positions ran over the range of cluster
    *sizes* instead, so it could not place a part beyond the first few
    positions;
  - the cost matrix was zero-initialised, so a state with no legal split
    counted as free and was preferred, giving partitions that violate
    the size bounds;
  - `udcl1` indexed its cost matrix with `(i + l) % n`, wrapping blocks
    around the end of the ordering and reading a row past the matrix.

  On TDA's own example it returned a partition of diameter 7 where one
  of 4 exists. Fixed, and checked in this package's tests against a
  brute force over every contiguous partition: 28 of 28 random cases (n
  = 4..8, 2 to 4 clusters, assorted size bounds) now agree exactly.
  `nc=`, `min=`, `max=` and `cn=` (the ordering) are command parameters,
  the result is a record table rather than a debugging trace, and it is
  exported as `ucl.part`. Reachable as `tda_ucl()`.
- `dsum()` (`t_gcmd.c`): its own comment says it exists "to force a and
  b to be stored prior to doing the addition ... where optimizers might
  hold one of these in a register", but the body was a plain `return a +
  b`, which every modern compiler inlines and folds. The
  machine-characteristics probe that depends on it was therefore
  measuring the compiler, not the hardware. Now uses `volatile`, the
  same device LAPACK uses in `dlamc3`, which defeats folding,
  fused-multiply-add contraction and reassociation alike.
- `S_XWIN` removed (`tda.h`, `tda_const.h`, `t_cmd.c`, `t_xstub.c`):
  `tda.h` defaulted it to 1 and `tda_const.h` to 0, which is where the
  `"S_XWIN" redefined` warning on every compiled file came from. It
  guarded an X11 front end, `t_xwin.c`, that is not in this source tree
  and never was — the only definition of `x_show` has always been a stub
  — and the guarded call passed no argument while the stub takes a
  context, so it could not have compiled had anyone enabled it. The
  `xshow` command remains, since TDA's help documents it, and now
  reports that the X11 front end is not part of the build instead of
  silently doing nothing.

  The original `t_xwin.c` still exists outside this tree, and the
  decision was checked against it rather than against its absence: it is
  X11-only (its Windows branch needs `t_xnt.h`, which is also missing),
  it predates the `TDAContext` refactor throughout, and it reads a
  PostScript file back and replays the primitives -- which is what
  `tda_read_ps()` and `tda_plot_ps()` already do, into R's own device,
  on every platform R runs on.
- `gdf()`'s joint estimation (`t_gdf.c`): **now usable from R, plus one
  small guard added.** The joint methods estimate a distribution over
  several related measurements per unit — say each person's time to two
  different events. The data has to be in "long" form (one row per
  person and measurement) with two grouping columns: who the row belongs
  to, and which of the measurements it is. The R wrapper only ever
  filled in the first, so any data with more than one row per person
  failed with a confusing TDA error. `tda_gdf()` now takes a `dimension`
  argument for the second column, and with everything in place the
  estimate checks out exactly: on uncensored data it reproduces the
  plain empirical distribution to machine precision. One thing to know
  when writing command files by hand: the name on the right-hand side of
  `gdf(...) = name` is where results are *written* — putting your input
  file's name there overwrites it. That mistake is also how we found
  that running `gdf` on an empty data matrix crashed TDA outright
  instead of giving an error; it now refuses cleanly ("no cases in the
  data matrix", `examples/coverage/err18.cf`).
- `mreg()` and `conj()` printed pages of leftover debugging output on
  every run — internal working values, dumped once per case per
  iteration, so a small monotone regression produced several hundred
  lines of noise around a dozen lines of actual results. Rohwer had
  commented most of these debug printouts out before release; a handful
  were missed. They are disabled now (kept in the source, switched off),
  the same treatment `zreg1`'s leftover debugging block got. The results
  themselves are unchanged; the two reference files were regenerated and
  differ only in the removed noise.
- `indep()` (`t_com.c`): **fixed — it crashed on ordinary input.** The
  command measures how far two variables are from independent: for each
  group of Y-values, the largest gap between the distribution of X
  within that group and overall. Another command left mid-study (no help
  entry, the general case never written). Its frequency table was
  addressed with the wrong row length and without subtracting the
  smallest values, so almost any real data wrote past the end of the
  table and brought the program down — inside R, the whole session.
  Indexing corrected throughout and the results now match direct
  computation in R exactly, including data whose values start well away
  from 1 (`examples/coverage/indep.cf`, comparison test in
  `test-r-comparisons.R`). Leftover debug tables switched off.
- `sdrel()` (`t_top.c`): **fixed — polygon containment never worked.**
  The command reports which spatial objects lie inside which others.
  Points inside polygons were found, but a polygon inside another
  polygon never was: the case simply wasn't written (it fell through to
  "not contained"), and two further slips in the boundary handling meant
  that even after wiring it up, an object touching the container's
  boundary was rejected under the option that is supposed to allow
  exactly that. All three corrected; the control case has one square
  strictly inside another (found by both options) and one square with an
  edge on the boundary (found only when boundaries count) —
  `examples/coverage/sdrelops.cf`.
- Coverage census (`tools/coverage_census.py`): a small tool that lists
  every command the program dispatches, checks which have no test case
  anywhere, and prints each missing command's own usage notes from the
  source — so test cases can be written in batches from one listing. The
  first batch added cases for `dump`, `dsplit`, `indep`, `nmca`, `gtopo`
  (with a small synthetic elevation grid), `rplz` (synthetic postal-code
  file), `sdrel`, `plsurf3d` and `plcirc3` — and found the `indep` and
  `sdrel` defects above on first contact. Of what remains uncovered,
  four commands are X11-screen-only (stubbed out in this build) and
  `rdxf` is commented out of the program itself.
- Option sweep (`tools/option_sweep.py`): runs every documented value of
  every enumerated option of every command that has a test case (158
  variants), one run per variant, and reports anything that aborts,
  prints new errors, or produces suspiciously little output. First full
  sweep: thirteen flags, twelve of them legitimate refusals (options
  that need data the fixture doesn't have — a symmetric matrix, grouped
  observations, more categories), and one real defect: `bfa(alg=0)`, the
  command's own documented default, was rejected by the shared option
  parser with the raw command text echoed as the error message. The
  parser accepts 0 now (0 is also what the option resets to, so nothing
  else changes behaviour); pinned in `boolfa.cf`, and `err14.cf` updated
  to the corrected behaviour.
- The X11 screen commands, replicated in R: `xreg` added a fitted curve
  to a screen scatterplot (lowess, least squares, or least absolute
  deviations, one per group) and `xconh` added convex hulls; their
  screen front end was never part of this source tree. Both are now
  available as `tda_pl_regression()` and `tda_pl_hull()` in the ordinary
  plot session. The replication is exact, not approximate: TDA's lowess
  and R's own `lowess(f, iter = 2, delta = 0)` are the same Cleveland
  algorithm, and the test draws the same curve through both engines and
  compares the resulting PostScript line point-for-point — the
  difference is zero. The least-squares line is ordinary regression; the
  L1 line uses TDA's own `l1reg` through its wrapper. The third screen
  command, `xplotf`, plotted values straight from a file, which
  `tda_pl_scatter()`/`tda_pl_lines()` with raw vectors already do.
- `ivreg()` and larger samples: **fixed as far as the mathematics
  allows.** The exact optimization (`opt=3`) searches over every
  observation the first two steps couldn't pin down, and that search
  grows exponentially with the number of such observations — inherent to
  exact interval regression, no algorithm changes that. What *was*
  wrong: the default effort limits were a flat `mxit=100, nbox=100`,
  which made the command give up even on problems it certifies in under
  a second (50 observations, say), printing an uncertified value with a
  bare `***` where the bound should be. The defaults now scale with the
  number of undecided observations — a handful certifies instantly,
  around a dozen within a minute, and hopeless problems still terminate
  (capped) instead of hanging. Values you set yourself are respected
  exactly as before. And whenever the search does stop at its limits —
  in `ivreg`, `ivreg1/2`'s exact step, and the other interval optimizers
  sharing the reporter, like `ivar` — the output now says plainly that
  the value is the best point found, not proven optimal, and that
  raising `mxit=`/`nbox=` may certify it, instead of the bare `***`. The
  same scaled defaults now apply to `ivreg1` and `ivreg2`'s exact step,
  whose shipped example had also been giving up uncertified at the old
  flat limits. Controls: `examples/coverage/ivhops.cf` (certified runs
  under the scaled defaults, plus a forced tiny-limit run pinning the
  message and that user limits are honored), `ivreg1d.cf` (ivreg1 under
  the scaled defaults; separate file because `mxit=` is a session
  setting and the main `ivreg1.cf` pins that user-set limits persist),
  regenerated `ivreg1.ref` and `ivar.ref`. A follow-up correction from a
  user example: `ivreg` searches ONE direction per run (its printed
  "bounds" enclose the minimum of the slope; `ns=1` gives the maximum),
  and `tda_ivreg()` had been labelling that single-endpoint enclosure
  "bounds on the slope" — an interval that need not even contain the
  true set of slopes. The wrapper now runs both directions for every
  method and reports a genuine outer interval; with `method = "exact"`
  the certified sharp bounds, flagged as such, and verified against
  exhaustive search. One more paper cut from the same user example:
  numeric options written into the generated command file used R's
  scientific notation (`nbox = 200000` arrived as `2e+05`, which TDA
  rejects with an unknown-parameter error), so large values needed an
  `as.integer()` workaround. Numbers are now formatted plainly at both
  entry points (`tda_block` and the shared `options=` normalizer), for
  every wrapper at once.
- `doc/matrix-wrapper-the maintainer's notes`: the matrix-language
  wrapper queue (mqap, mlp, the constrained least-squares family, mqp)
  written up for delegation -- semantics established so far, the
  exporter/doc/test requirements in full, the traps, and the gate list.
  The base-R-equivalent commands stay unwrapped by decision.
- Cross-platform certification (found by a user's macOS run): whether
  any box passes the final width/tolerance acceptance sits within a few
  ulps of `tolbw` and flips between compiler/FPU combinations, so the
  same search certified on one machine and printed `***` on another --
  with the identical, correct optimum. The reporter now recognizes the
  exhausted state (`no unprocessed boxes` and `no survivors`) for what
  it provably is: every region was excluded against the running best,
  which only improved afterwards and is attained, so the best value IS
  the certified optimum -- "(certified by exhaustion)". The suite's own
  `ivar.cf` turned out to sit in exactly this state and now certifies at
  20.1369 / 65.2075 instead of `***`; reference regenerated, diff shows
  only the new certificate lines. The numeric-options regression test is
  rescoped to what it is about -- both entry paths producing identical
  runs -- rather than the platform-marginal certificate.
- `tda_interval_wages()`: a reproducible example dataset for the whole
  interval family -- schooling asked in whole years (a point), wages in
  widening brackets with a top-coded highest class, the classic survey
  mix. On it, `tda_ivreg`'s exact method certifies sharp bounds
  instantly BECAUSE the regressor is exact: the fixing steps decide
  every observation, which the wrapper now recognizes as a certificate
  ("Remaining coordinates: 0" with a collapsed bracket). The same shape
  scales: 10,000 observations with 15 top-coded incomes certify in 0.1 s
  -- the realistic large-n case for the exact method, since the search
  is exponential only in the number of interval-valued observations, not
  in n. The dataset also cross-validates the engine: with a point
  regressor, `ivreg`'s exact path and `ilsreg`'s closed form compute the
  same estimand through unrelated code, and the suite pins their
  agreement. Two more honesty fixes fell out of building it:
  `icov`/`icorr` and `ivar`'s shared `f_range` used to abort on box
  overflow without reporting the best point found; both now report it,
  uncertified, and continue with the other direction.
- Three more wrappers for the newly covered commands, exporter-first:
  `tda_locate_line()` (`gloc`, one-dimensional multifacility location:
  weight matrices in, optimal fixed-point assignment out; new
  `gloc.assignment` exporter beside the print), `tda_boolean_min()`
  (`bfa`, Quine-McCluskey boolean minimization of an outcome over up to
  15 binary conditions -- the machinery behind crisp-set QCA; new
  `bfa.selection` exporter encodes each cover term as 1/0/NA per
  condition), and `tda_mlrc_design()` (`mldes` via the existing
  `mpr.matrix` exporter). Each is verified against an independent
  computation in the suite: brute-force assignment enumeration,
  truth-table reproduction of the minimal cover, and the traced
  pairwise-product formula. A fuzzing harness (`tools/fuzz_wrappers.R`)
  sweeps the newer wrappers across option grids and edge-shaped inputs
  looking for dead ends; its one finding -- infeasible `mlpi`
  constraints surfacing as a raw command echo -- now comes back as a
  clean refusal ("the constraints cannot be satisfied by any 0-1
  vector").
- User interrupts: a long branch-and-bound run (`igmin`, `ivreg_min` --
  the loops behind `ivreg(opt=3)`, `ivar`, `icov`, `icorr` and friends)
  can now be stopped with Esc/Ctrl-C from R. The package-only helper
  `tda_interrupt.c` detects a pending interrupt with the
  `R_ToplevelExec`/`R_CheckUserInterrupt` idiom -- detection without the
  longjmp that would skip cleanup -- and the loops end as if the
  iteration limit were reached: ordinary bookkeeping, everything freed,
  an "Interrupted by user." line in the output, which `tda_run()`
  re-raises as a real R interrupt condition (class `tdaInterrupt`), so
  `tryCatch(interrupt = ...)` works as for any long R computation. The
  standalone binary is untouched (Ctrl-C from the OS, as always). Tested
  live: SIGINT lands within the polling throttle, and the session works
  normally afterwards.
- New commands `icov` and `icorr` (requested; explicitly marked as
  additions, not Rohwer's): covariance and correlation of two
  interval-valued variables, same estimand and machinery as `ivar`
  (igmin branch-and-bound, typ 4 and 5), with interval quotient-rule
  derivatives for the monotonicity tests. Measured envelope: covariance
  certifies at real sizes -- the patients example is exact against all 4
  million corner combinations; correlation certifies small problems (n=3
  instantly, exact against a dense grid) and beyond that stops at its
  limits with the honest message -- its direct-form inclusion is the
  limit, a centered form would be the next step. Wrappers
  `tda_icov()`/`tda_icorr()` (`~ iv(xlo,xhi) + iv(ylo,yhi)`;
  command-style names matching `tda_imean`/`tda_ivar`/`tda_igini` (the
  variance names now follow the commands: `tda_ivar` is `ivar`, the
  tunable branch-and-bound -- identical to `tda_ivariance`, which stays
  -- and `tda_ivar1` is `ivar1`, the newer algorithm, whose wrapper
  previously returned NULL fields because nothing parsed its output
  format; it now fills `$bounds`, and the suite pins that the two
  algorithms agree on the patients data); naming the same pair twice
  triggers a warning that the two copies vary independently -- that is
  not the variance), coverage case `examples/coverage/icovops.cf`,
  brute-force pins in the suite. The text parser question that came with
  this: the parser is the deliberate exporters-off fallback, exercised
  by the flag-off suite to prove the wrappers agree with a stock TDA
  binary's printed output; exporters are always primary.
- Interval-wrapper corrections from user testing: the heuristic method
  returned an iteration table's coordinate COUNT as the lower slope
  bound (its 3-column table was parsed with the two-step layout); it now
  reports its achieved values through a new `iv.heuristic` exporter and
  is labelled for what it is -- an inner approximation, nested the other
  way around from two-step's outer interval, which the suite now
  asserts. TDA's "no certified bound" and warning lines are raised as R
  warnings by `tda_ivreg` and the interval statistics instead of sitting
  unseen in `$run$output`. `tda_ivariance` gains `$sd` (square root
  transfers to the bounds; TDA itself has no interval sd, covariance or
  correlation command). `ilsreg` and `inpreg` were audited for the ivreg
  pathology and are clean: ilsreg computes both directions in closed
  form, inpreg is a fixed-point envelope with no branch-and-bound.
- Three new wrappers with clean interfaces, exporter-first as always:
  `tda_independence()` (the `indep` delta-independence measure, per
  y-subset), `tda_isotonic()` (`mmp`/`mmp1`/`mmp2`; plain PAVA equals
  `stats::isoreg` exactly, and the two Kruskal tie approaches are named
  `ties = "primary"/"secondary"` instead of digits), and `tda_mlpi()`
  (`mlpi` as maximize-subject-to with named arguments, no `opt=`
  anywhere). New C exporters `indep.deltas`, `mlpi.value` and
  `mlpi.solutions` sit directly beside the prints so the two can never
  drift; text parsing remains only as the exporters-off fallback, which
  the flag-off suite exercises. Each wrapper is verified against an
  independent computation in `test-r-comparisons.R`: `isoreg`,
  exhaustive 0-1 search, and the direct conditional-distribution
  formula.
- `ivar` is pinned against an exhaustive corner search on the two
  interval-variance datasets printed in Gioia & Lauro (2005), in
  `examples/coverage/gioia05.cf`; TDA's bounds are exact to every digit.
- Verification and documentation pass for the interval work: the `ivreg`
  bounds are now pinned against exhaustive search (the slope is linear
  in each y for fixed x, so y-extrema sit at interval corners; x is
  searched on a grid) — TDA's certified bounds enclose every brute-force
  slope and sit within grid resolution of the sharp extrema, on
  randomized instances, in `test-r-comparisons.R`. Three things this
  work established that the source never said are now written into it:
  the `[lo,hi]` interval literals in matrix expressions (a second value
  plane, printed as pairs — was completely undocumented), `mlpi`'s
  objective (maximize, with the first rows of A and B as objective and
  constant — the source stated none), and `m_mmp`'s header (the shipped
  comment misnumbered its own options and swapped mmp1/mmp2).
  `ivreg_h`'s header now states the slope formula and marks `ivreg_h1`
  as the dead preliminary draft it is.
- `rcsv()` (`t_sdx.c`): `err` was set to `-1` at the top of the function
  and never cleared, so the command reported failure on every run,
  including the ones that converted the file correctly. Set to 0 once
  the whole input has been read. Note also that `df=` is documented as
  an output file but the conversion is written through `tda_out()`, so
  no file appears — left as Rohwer wrote it, recorded here.
- `gloc()` (`t_ass.c`): **finished — it never worked as shipped.** The
  command is meant to solve a facility-location puzzle: given some fixed
  places along a line and weights saying how strongly each new facility
  should sit near each place (and near the other new facilities), find
  the best spot for every new facility. The solver underneath was
  complete and correct, but the command never read any input — it parsed
  its options, threw them away, and solved one small example hard-coded
  into the source, so every call returned the same answer whatever you
  gave it. It now reads the problem from a graph: the first `n=` nodes
  are the new facilities, the remaining nodes the fixed places in their
  order along the line, with the weights on the edges. Verified against
  exhaustive search on randomized problems, at two different spacings of
  the fixed places (the algorithm's claim that only their *order*
  matters holds). Control files: `examples/coverage/glocops.cf` (working
  cases, including the single-facility shortcut), `err19.cf` (refusal
  without `n=`), the rewritten smoke case `examples/tests/gloc.cf`, and
  the brute-force comparison in `tdaR`'s `test-r-comparisons.R`. Its
  leftover debug printouts were switched off along the way, like
  `mreg`'s.
- ZOO archives, both directions:
  - `lzd_encode()` added to `tda_zoo.c` — a port of `lzc()` from zoo's
    own `lzc.c` (Rahul Dhesi, 1986, released by him to the public
    domain). `zoo()` gains `method =`: 1 (LZD) by default, 0 to store.
    TDA's `arcd` refuses anything outside its accepted set, so stored
    archives were the one kind TDA could not open — writing LZD fixes
    that from our side.
  - `lzs()` added to `t_zoo.c` and method 0 accepted in `t_rzoo.c`. TDA
    only ever consumed archives made by Dhesi's `zoo`, which compresses,
    so it had never met a stored member and rejected it as an "unknown
    packing method" — while every other zoo reader takes one. All three
    methods are now accepted: 0 stored, 1 LZD, 2 LZH.
  - LZH decoding itself (`decode`, `lzh_decode`, `make_table`,
    `getbits`, and the rest of that family in `t_zoo.c`) sat at 0%
    coverage even after the above, because `tools/make_archive.py` (this
    project's fixture generator) can only write stored or LZD archives.
    It is now exercised with `deha1.zoo`, the example archive from the
    TDA teaching pages at stat.rub.de, which is not part of this
    repository: it lives with the other external fixtures
    (`tests/fixtures/README.md`), together with `arcd_lzh.cf`, which
    reads a real LZH member through `arcd`/`arcc`/`arcv`/`arcvc`, and
    `build_deha1_aug.py`, which appends the type-2 variable-description
    member `arcd` needs to a byte-for-byte untouched copy of the archive
    (`frd_dir()`/`b_to_dir()` in `t_zoo.c` parse only a handful of
    fixed-offset fields, and `next` is an absolute file position, so the
    appended entry needs no directory rewrite). With that fixture
    `t_zoo.c` goes 40.8% → 73.5% and `t_rzoo.c` (the archive-reading
    side, via `nvar(... = A:VarName, ...)`) 40.5% → 61.7%; `get_astr` is
    the one function left at 0%, and it has no caller anywhere in the
    sources, the same class of dead code as `bfa_lm2` elsewhere in this
    document.
- `psetupg()` projection 20 (the globe-view projection), `t_map.c`:
  **fixed.** This projection shows the earth as a globe seen from space,
  so at most half of it is ever visible. TDA accepted map regions wider
  than that. Ask for a region 150 degrees wide and it quietly drew a map
  only 30 degrees wide instead; ask for 180 and the drawing came out as
  astronomically large garbage numbers. Now any region wider than 90
  degrees to each side is refused with an error, which is the most a
  globe view can show anyway. An earlier note here claimed TDA also
  wrongly *refused* some valid regions — that was a misreading: those
  regions really don't fit on the visible globe, and refusing them is
  correct. The other three map projections were fine all along.
  Regression tests in `test-api.R`.
- `tda_ivreg2(method = "minimizer")`: **our bug, now fixed — not
  TDA's.** It returned `NA` coefficients, which I first reported as a
  TDA defect. TDA had in fact converged: it prints a point estimate

      Convergence reached in 4 iterations.
      FMIN=-1.25
      par=      0.500000000000       0.000000000000

  rather than the bounds table the other methods print, and the wrapper
  parsed only the table. `.iv_beta()` reads `par=` now, and the slope
  comes back correct (0.5 on the documented example, matching `search`).
  `alpha` is deliberately left unavailable for this method: `par[2]` is
  0 where the intercept is 2, so it is not the intercept, and a wrong
  number is worse than an honest `NA`.

  `contour` needing a current PostScript file (it draws) is expected
  rather than a defect, but is not documented as such.
- Stale documentation found while raising coverage, recorded not fixed:
  `intp`'s header comment documents an `ncp=` parameter ("must be >= 2
  and <= 25") but the function assigns `ncp = 3` as a literal and
  rejects the option as a syntax error. Likewise `gcnt(X,G)` and a bare
  `mr` are syntax errors in the form the `t_eval2.c` opcode comments
  imply. In each case either the comment is stale or the parser is.
  `m_cat`'s own syntax box (`t_matc.c`) names its opt=2 command
  `mdcathv`; the dispatcher (`t_mat.c`) registers it as `mcathv`. Pinned
  as `mcathv` — the working name — in `examples/coverage/matc4.cf`.
  `m_mmp`'s own header comment (`t_matc.c`) numbers its opt values 1/2/3
  and names opt=3 `mmp1` a second time (should be `mmp2`); the
  dispatcher (`t_mat.c`) actually uses opt=0/1/2 for `mmp`/`mmp1`/
  `mmp2` respectively. Pinned with the correct names in
  `examples/coverage/mmpops.cf`. `gcni`'s manual entry (tda.hlp)
  documents its secondary output as `df=...`; the real keyword, found in
  `t_parm.c` (`gcni` itself is undocumented in tda.hlp at all -- found
  only via its own source header comment), is `prot=`. Pinned as `prot=`
  in `examples/coverage/gcniops.cf`.
- `mdsn1` was unreachable, and it is a real command — **two bugs, both
  fixed**. It is the censored counterpart of `mdsn`, as `lsreg1` is of
  `lsreg`.
  - `t_cmd.c` had no dispatcher entry for it, and the test for `mdsn` is
    a four-character prefix match, so `mdsn1(...)` matched `mdsn`, ran
    the wrong command, and then reported a syntax error on the leftover
    `1(...)`. `lsreg1` is ordered before `lsreg` a few lines above;
    `mdsn1` was simply missing. Entry added, ordered first.
  - `mdsn1` then parsed its options from `CmdBuf + 4` rather than `+ 5`,
    leaving the trailing `1` in front of the option block, so every call
    was still a syntax error. Offset corrected.

  It now runs: 0% -> 36.1% of its 238 lines. Found by combining gcov
  with a call-site count while raising coverage, and by not accepting
  "dead code" as the explanation.
- Dead code, recorded not removed: `bfa_lm2` (248 lines) has no call
  site anywhere in the sources, so no input can execute it. Its
  neighbours `bfa_lm` and `bfa_lm1` are live and now covered (89% and
  84%): they are only reached when the essential prime implicants fail
  to cover the chart, which needs a *cyclic* function such as `ab' + bc'
  + ca'` — see `examples/coverage/boolcyc.cf`. `bfa_lm2` stays at 0%
  against that same fixture, which is what makes it dead rather than
  merely untested.
- Out-of-range option values: **not a defect after all — the earlier
  note here was wrong.** When an option gets a value TDA doesn't know
  (say `alg=99` where only 1–7 exist), TDA doesn't error and doesn't
  silently do nothing either: it falls back to the default and says so
  in the output ("Algorithm 1: principal components"). That's how Rohwer
  built it, consistently across commands, and it stays that way. The R
  side is stricter where it can be: `tda_dma()` was the one wrapper
  letting a raw `alg` value through and now rejects anything outside
  1–7. The cases where TDA does refuse bad input are pinned in
  `examples/coverage/err*.cf`.
- `dnum()` (`t_spss.c`), portable-format reader: **known, not changed.**
  The base-30 mantissa is accumulated into a `double` and rescaled at
  the end. A double holds integers exactly only to 2^53 = 9.0e15, while
  the thirteen digits the function itself allows reach 30^13 = 1.6e19,
  so the low bits are lost before the rescale. Costs a few ulp on some
  values; `1e6` reads back about 2e-10 low.

  `uint64_t` covers 1.8e19 and fixes it exactly — but that type did not
  exist when this was written (C99, 1999; the file is 1989-97), and the
  affected digits are below anything a statistical result depends on.
  Left as Rohwer wrote it. The same code was ported into the readspss R
  package, where the fix has been measured: exact round-trips go from
  1136/3003 to 2301/3003 and the worst relative error from 7.5e-16 to
  3.2e-16.
- Portability, behaviour-neutral: `open_memstream()` gets a `tmpfile()`
  fallback on Windows/Cygwin, and `qsort_r` variants are wrapped by
  `tda_compat.h`.
- `mag(A,C,R,B)` (matrix aggregation, `t_matc.c`): **fixed — gave the
  wrong answer depending on what else had run before it.** `C` and `R`
  are the row- and column-group vectors the command has just parsed into
  its own scratch buffers, but the code read them back from the
  persistent named-matrix table instead — specifically, from whichever
  matrix happened to be the 2nd and 3rd one ever defined in the whole
  session, not from `C` and `R` themselves. Call `mag()` first thing and
  it can work by pure coincidence (if `A`, `C`, `R` are the first three
  matrices `mdef`'d); run any other matrix command first and it silently
  returns zeros instead of the aggregate. Found because `mag` had 0%
  coverage before this session — nothing had ever exercised it, so the
  bug had never been seen. Fixed to read `C` and `R` from its own
  scratch buffers, matching how every other command in the file reads
  its arguments; checked against a hand aggregation, pinned in
  `examples/coverage/mag.cf`.
- `mdsn` and `mdsn1` (nonmetric MDS, `t_mds.c`): **fixed — did nothing
  at all unless you passed `ns=` yourself, and then printed garbage.**
  Both fit their configuration inside a `for (ii = 1; ii <= PMNS; ...)`
  repeat loop, and `PMNS` (the number of random restarts) defaults to
  `-1` until an `ns=` option sets it. The code that should clamp it to
  at least 1 — present and working in the sibling `mdsc`/`mdsm`
  functions — was sitting inside a large commented-out block in both
  `mdsn` and `mdsn1`, so by default the loop bound was `-1`: the loop
  never ran, no optimization ever happened, and the function printed
  "Number of random repeats: -1" followed by whatever garbage happened
  to be on the stack for the stress value and the final coordinates —
  the same uninitialized value every time in one build (`4.66308e-310`),
  a different one in another. Found because `mdsn` and `mdsn1` had 0%
  coverage before this session — nothing had ever called them without
  `ns=`, so nothing had ever seen the -1. Fixed by restoring the clamp
  (`if (PMNS < 1) PMNS = 1`), matching the pattern already active in
  `mdsc`/`mdsm`; the redundant `get_dsv()` call the original disabled
  block also had is not restored — neither function otherwise wires up
  custom starting values, and reintroducing an untested call wired to
  locals that no longer exist in `mdsn`'s current signature would trade
  one guess for another. Confirmed by running: real, distinct stress
  values from an actual gradient search, not the same repeated garbage
  float across every call. Pinned in `examples/coverage/mdsops.cf`.
- `ivreg2` (interval-data regression, `t_ireg.c`): **two bugs, both
  fixed.**
  - The final summary — "Best function value", "Alpha (center-radius)",
    "Beta (center-radius)" — used three local variables (`f`, `bc`,
    `br`) that are only ever assigned inside the `opt=2` branch. For
    `opt=1` (direct calculation, the default and by far the most common
    case), execution reaches that summary with all three still
    uninitialized: confirmed directly, a literal `-nan` printed as the
    "Best function value", and a `Beta` of `0.0` that flatly
    contradicted the correctly-computed `Parameter: 0.51` the command
    had printed two lines earlier in the same run. `opt=1`'s real answer
    was never missing — it's the per-domain table just above, with the
    valid domain(s) marked `*` — this summary was a garbage duplicate
    nobody had reason to trust. Nothing had caught it because the R
    wrapper's own test for `method = "direct"` only checks that the run
    doesn't error, never these specific numbers. Fixed by scoping the
    summary to `opt=2` only, the one path that actually computes it.
  - `opt=3` ("mina") called a bare `exit(0)` on its own *success* path —
    not an error path, the ordinary case where the simplex minimizer
    converges cleanly. That killed the whole TDA process on the spot:
    not just this command, but silently discarding every command after
    it in the same or any enclosing command file, with no message, and
    skipping normal cleanup entirely. Two sibling debug `exit(0)` calls
    elsewhere in the same file are already properly disabled inside `/**
    ... **/` comment blocks; this one was simply missed. Confirmed by
    running a plain `time;` command right after `ivreg2(opt=3)`: it
    never executed before the fix, and does after. Fixed by returning to
    the caller (`err = 0; goto IVREGFin;`) instead of exiting.
- `unf`/`mdsn_grad`/`mds_m` (`t_mds.c`): **fixed — three separate
  leftover debug `exit(0)` calls, all killing the whole process
  mid-command.** Same failure class as `ivreg2(opt=3)` above, found
  while chasing a coverage gap rather than by symptom.
  - `unf1()` (called from `unf`, the rank-unfolding command, for every
    case it processes, not just when `df=` requests a second output
    file) had `tda_out("vor unf2\n"); ...; exit(0);` sitting directly
    before an already-correct `if (d > c) {...}` block, so it died on
    the very first case, every time. The project's own existing pinned
    test (`examples/tests/unf.cf`) had captured this truncated, crashed
    output as its "expected" reference — the same silent-truncation trap
    as the `dmet1` bug elsewhere in this document. Confirmed by running
    the fix: before, no completion trailer at all; after, `unf`
    correctly reports "Best permutation" and finishes normally. Fixed by
    deleting the debug prints and the `exit(0)`, leaving the `if (d >
    c)` check that was already there.
  - `mdsn_grad()` had `tda_out("HIER\n"); exit(0);` directly before an
    already-written `continue;` guarding a degenerate-distance case —
    the `exit(0)` made that `continue` permanently unreachable. Fixed by
    deleting the debug print and the `exit(0)`.
  - `mds_m()` had two `exit(0)` calls on `ffmin()` failure, in a
    function that returns `int` and whose own caller (`mdsxr_m`) was
    *already written* to handle a nonzero return correctly (`if (rn) {
    printf1(...); goto MDSXFin; }`) — confirmed by reading that caller
    directly before touching anything. Fixed both to `return(-1);`,
    matching both the caller's existing expectations and the `goto
    ...Fin` convention this same file already uses for the identical
    `ffmin()` check elsewhere (`t_mds.c` line ~2621). One of the two
    `exit(0)` calls sat inside a block guarded by `fflag = 0;` set
    unconditionally two lines above it — dead code regardless, fixed
    anyway for consistency since the change was already in hand.
- `ivreg` (interval-data regression, `t_ireg.c`): **fixed — `opt=2`
  ("heuristic optimization") hung forever, confirmed reproducible on
  every `n=` block size tried and on two unrelated datasets.**
  `ivreg_z()` tracks which of a parameter's two slots (center, radius)
  is still open with a single encoded integer: `-i` for the center slot,
  `i` for the radius slot. At `i = 0` both encodings produce the same
  value, plain `0`, and every place that decodes it (`if (j <= 0) ...
  else ...`) always reads a `0` as the center slot. So whenever the
  radius slot at `i = 0` was the one that still needed narrowing, it
  silently got treated as the (already-narrow) center slot instead and
  never actually got fixed — confirmed directly with a temporary trace:
  the open-parameter count dropped 6, 5, 4, 3, 2, then stuck at 1
  forever, with that "stuck" parameter's own bounds already
  bit-identical (i.e. some other, uncounted parameter was the real
  holdout). The outer loop's own exit condition never saw that count
  reach zero, so it never stopped. Fixed by shifting the radius encoding
  to `i + 1` (never `0`, so it can no longer collide with the center
  slot's `-i`) and updating the three places that decode it to match.
  Confirmed fixed across `n = 1, 2, 3, 4` and on both datasets that hung
  identically before.

  Both checked against data taken directly from `tdaR`'s own roxygen
  example (`R/interval.R`), pinned in `examples/coverage/ivreg2.cf`
  (`opt=1,2,3`) and `examples/coverage/ivreg1.cf` (the box-search
  variant, which runs cleanly but — as tdaR's own docs already warn —
  doesn't complete under TDA's default iteration/box limits).
- First of the matrix-language wrappers: `tda_mqap()` (`mqap`, quadratic
  assignment, CACM algorithm 608). The objective was already stated in
  `qap_w()`'s own header, so no semantics work was needed there; a new
  `mqap.value` exporter sits beside the existing print, matching
  `mlpi.value`'s pattern, while the permutation itself goes through the
  already-generic `mpr` exporter, same as `tda_isotonic`. Verified two
  ways in `test-r-comparisons.R`: the returned permutation's objective,
  recomputed independently in R, matches TDA's own value on every random
  instance; and at `n` up to 6, exhaustive enumeration over all
  permutations confirms the heuristic never beats the true optimum, and
  reaches it on the pinned seed. Pinned separately in
  `examples/coverage/mqapops.cf`, including the nonzero-diagonal warning
  path. An all-degenerate-input case (flows, distances and costs all
  zero, or costs zero and one of flows/distances zero) is refused
  directly in R instead of leaking TDA's own "matrices are zero" text,
  caught by the fuzz harness (`tools/fuzz_wrappers.R`), which now also
  sweeps `mqap` across several sizes, with and without costs, and the
  refusal cases.
- Second: `tda_mlp()`/`tda_mlp1()` (`mlp`/`mlp1`, Salazar & Sen's MINIT,
  CACM 333). Unlike `mqap`, the tableau layout and both the primal and
  dual problems are already fully stated in `lpf1()`'s own German header
  comment in `t_lp.c` -- confirmed by reading the actual indexing code
  rather than trusting the comment blindly, since the header's own
  stride formula has a typo (`(m+1)` where the code actually uses `n1 =
  n+1`). `mlp1(T,p,X,Y)` treats the *last* `p` constraint rows as
  equalities, whose dual entries TDA fixes at 0 -- also stated in the
  header and confirmed by hand on a small instance. New scalar exporter
  `mlp.value` sits beside the existing "Value is:" print; that print
  itself uses the wrong printf level (`printf2`, visible only under
  `silent=-1`, unlike every sibling command's `printf1`) -- left as
  stock TDA's own behaviour rather than changed, with the wrapper
  sending `silent=-1` itself so its exports-off fallback still sees the
  line. Calling `mpr` twice in one run (for `X` then `Y`) exports under
  `mpr.matrix` and `mpr.matrix.2` respectively -- TDA's own generic
  per-run numbering, not a new mechanism. `mlp`'s own failure messages
  ("No solution.", "...is unbounded.") never start with "Error" and
  bypass the usual `m_cmdmsg` convention entirely, so the wrapper
  matches on the specific text instead, most specific first (the
  unbounded cases both also contain "No solution."). Verified against
  brute-force vertex enumeration over random small LPs in
  `test-r-comparisons.R` (every vertex of `{x>=0, Ax<=b}` sets `n` of
  the `m+n` constraints tight; solve and filter for feasibility), plus
  two hand-checked instances (plain and one equality) pinned in
  `examples/coverage/mlpops.cf`. The fuzz harness needed its own
  deliberate-refusal detector widened to recognize
  "infeasible"/"unbounded" phrasing -- a fix to the harness, not the
  wrapper, once it was clear the three flagged cases were genuine
  refusals worded differently from earlier wrappers'.
- Third, the constrained least-squares family:
  `tda_mls()`/`tda_mlse()`/`tda_mlsi()`/`tda_mlsei()`/`tda_mlsei1()`/`tda_mnls()`,
  unifying `mls`, `mlse`, `mlsi`, `mlsei`/`mlsei1`, and `mnls`. An
  earlier working note (`doc/matrix-wrapper-the maintainer's notes`) had
  assumed `mlse(S,me,R)`/`mlsi(S,mi,R)` took their own constraint
  counts; `tda.hlp`'s own `##mlse`/`##mlsi` entries say otherwise, and
  reading `m_mls`'s dispatch confirms it: `mlse` and `mls` are the
  literal same `lsei()` call (confirmed empirically too -- an
  overdetermined, inconsistent instance returns the plain least-squares
  fit rather than erroring), and `mlsi` treats its whole input as
  inequality rows. Every combination the wrapper exposes collapses to
  one `mlsei(S,me,mi,R)` call (`me`/`mi` zero when a block is absent),
  except `nonneg = TRUE`, which needs `mnls` (TDA's own non-negativity
  command, with no general inequality block of its own -- refused
  together with `inequalities`). `mlsei1`, undocumented in `tda.hlp`,
  turned out to be an alternate active-set QP solver (Powell/
  Schittkowski `qld`) for the identical E/A/H problem as `mlsei` --
  confirmed by matching both algorithms to solver tolerance on two
  independent instances, not assumed from the similar name, and exposed
  as `method = "active_set"`. New exporters `mls.rank` and `mls.rnorm`
  sit beside the existing prints in the shared `m_mls` path; `mlsei1`
  gets neither (it prints no residual or rank at all), so the wrapper
  always recomputes the residual itself in R from `A`, `x` and `b`
  rather than trusting any one command's own text. Verified three ways
  in `test-r-comparisons.R`, per the handover's own prescription:
  unconstrained against `qr.solve`; equality- constrained against the
  direct Lagrange/KKT block solve; inequality and `nonneg` against brute
  force over every active set (each subset of inequality rows tried as
  tight, solved, filtered for feasibility, best kept) -- plus a direct
  check that `method = "active_set"` agrees with the default on the same
  mixed instance. Pinned in `examples/coverage/mlsops.cf` (the four
  `m_mls`-dispatched commands) and `examples/coverage/mlsei1ops.cf` (the
  fifth, separately, since it is a different function entirely).
- Fourth and last of this batch: `tda_mqp()`/`tda_mqpb()`/`tda_mqpc()`
  (`mqp`/`mqpb`/`mqpc`), also undocumented in `tda.hlp` -- identified
  entirely from `qld()`'s own header in `t_qp.c`, the same
  Powell/Schittkowski `ZQPCVX` active-set solver as `mlsei1`. `qld`
  solves `min grad'x + 0.5 x'Gx s.t. A(k)x=b(k)` for the first `meq`
  rows, `A(k)x>=b(k)` for the rest, `xl<=x<=xu`; reading all three
  branches of `m_mqp` (rather than assuming from the similar names)
  shows `mqp(C,D,X)` is the unconstrained case (box bounds fixed at
  +-huge), `mqpb` adds user box bounds, and `mqpc` adds general
  equality/inequality constraints instead -- no command combines box
  bounds with general constraints, so the wrapper refuses that
  combination outright rather than picking one arbitrarily. No C changes
  were needed for this family: it prints nothing beyond an occasional
  Cholesky-enlargement diagnostic, so the wrapper computes the objective
  value itself in R from `C`, `d`, and the returned `x` rather than
  depending on any TDA output. Verified against hand- solved instances
  for all three commands (unconstrained, box-clipped, and mixed
  equality+inequality) pinned in `examples/coverage/mqpops2.cf`, plus
  against closed-form and active-set-brute-force solves in
  `test-r-comparisons.R` (same KKT-block recipe as the
  constrained-least-squares checks, adapted for a quadratic rather than
  a least-squares objective). This closes out every command family in
  `doc/matrix-wrapper-the maintainer's notes`.
- Renamed and split all of the above (plus the pre-existing
  `tda_binary_program()`) to match TDA's own command names one-to-one:
  `tda_mqap`, `tda_mlp`/`tda_mlp1`, `tda_mls`/`tda_mlse`/`tda_mlsi`/
  `tda_mlsei`/`tda_mlsei1`/`tda_mnls`, `tda_mqp`/`tda_mqpb`/`tda_mqpc`,
  `tda_mlpi`. The unified multi-mode functions above were a defensible
  design on their own, but left the package with two unrelated naming
  conventions side by side -- friendly names for some matrix-language
  commands, raw TDA tokens (via `tda_mat`) for the rest -- which made
  neither half easy to find from the other. Each new function keeps its
  predecessor's argument names and validation (never TDA's own tableau
  packing), just split one-function-per-TDA-command and named
  accordingly; the shared implementation moved into non-exported
  `.tda_ls_impl`/`.tda_qp_impl` so the six/three-way splits don't
  duplicate the file-writing and export-reading logic. Documentation is
  grouped by family via roxygen `@rdname` (one Rd page for `tda_mls`'s
  six functions, one for `tda_mqp`'s three, one for `tda_mlp`'s two)
  rather than one page per function. All previously hand-verified
  numeric results were re-checked unchanged under the new names before
  any test file was touched.
- Started wrapping the rest of TDA's matrix language one command at a
  time (per explicit request, beyond `doc/matrix-wrapper-the
  maintainer's notes`'s original scope): every command gets its own
  `tda_m...()` function named after it, grouped into a handful of Rd
  pages by family rather than one page per function, built as thin calls
  into the pre-existing generic `tda_mat(op, ...)` rather than
  reimplementing the file-writing/run/export-reading machinery per
  command. First batch, core linear algebra: `tda_mmul`, `tda_mtransp`,
  `tda_mchol`, `tda_minvs`, `tda_minvd`, `tda_mginv`, `tda_mnrow`,
  `tda_mncol`, `tda_mnorm`/`tda_mnorm1`/`tda_mnorm2`, `tda_mtrace`,
  `tda_mdiag`, `tda_mdiagd`. Every one of these has a base-R equivalent,
  so verification is direct comparison (`solve`, `t(chol(.))`, `t()`,
  `diag`, ...) rather than hand computation.
- **Found and fixed a real bug in `mnorm` along the way**, caught by
  this batch's own randomized verification rather than looked for:
  `tda.hlp` documents `mnorm` as "the absolute value of the matrix
  element that has largest absolute value", but `m_nrow`'s case-2 branch
  computed `dmax(0, raw_value)` -- no `fabs()`. For any matrix whose
  most-negative element outweighs its most-positive one, this silently
  returned the wrong number; for an all-negative matrix it silently
  returned 0, since the running max never starts below 0 and no negative
  value ever exceeds it. `mnorm1` (sum of absolute values) and `mnorm2`
  (Frobenius norm) were checked too and were already correct --
  confirmed by reading, not assumed from the shared dispatch function.
  Fixed by adding `fabs()`, in both `src/` and `tdaR/src/`. Pinned in
  `examples/coverage/mnrowops.cf`, including the all-negative case
  specifically (the more dangerous, previously silent failure mode, not
  just the milder mixed-sign one that first caught it).
- `tda_mchol`/`tda_minvs` (positive-definite required) and `tda_mginv`
  (at least as many rows as columns required) pre-validate in R with a
  clear message before ever calling TDA, using base R's own `chol` as
  the positive-definiteness check -- rather than letting a malformed
  call reach TDA and leak its terse `tda_mat`-generic "Error in command:
  ..." text, which the fuzz harness flagged as a real, if minor, UX gap.
- Batch 2, vectorizing/reshaping/sorting: `tda_mcvec`/`tda_mrvec`
  (`mcvec`/`mrvec`, the classic `vec(A)` operator and its row-major
  counterpart -- `as.vector(A)`/`as.vector(t(A))`), `tda_mivec`
  (`mivec`, `vec`'s inverse -- `matrix(x, nrow = n)`), `tda_mrsum`/
  `tda_mcsum`/`tda_mdrow`/`tda_mdcol` (`mrsum`/`mcsum`/`mdrow`/ `mdcol`
  -- row/column sums as vectors or as diagonal matrices),
  `tda_msort`/`tda_msort1`/`tda_mrank` (`msort`/`msort1`/`mrank`),
  `tda_mcath`/`tda_mcatv`/`tda_mcathv` (`mcath`/`mcatv`/`mcathv` --
  `cbind`/`rbind`/block-diagonal direct sum), `tda_mtrim` (`mtrim`, crop
  or zero-pad whole rows/columns at either edge). `mrank` is the one
  real naming trap in this batch: despite the name, it is **not**
  linear-algebra matrix rank -- `t_mat.c` dispatches it through the
  exact same function as `msort`/`msort1`, and reading (then hand-
  checking) that function shows it returns the sort *permutation* itself
  (`order()`), not each row's rank position. No TDA command returns a
  linear-algebra rank as a matrix value at all (`mginv` only prints its
  pseudorank, never returns it). `mivec`'s fill order was also worth
  double-checking rather than trusting a first read of the loop nesting:
  it turned out to be column-major (`vec`'s own inverse), not row-major
  as the loop structure suggested on a quick look -- caught by
  hand-checking a small instance before writing any documentation,
  exactly the case this project's house rules exist for.
- Batch 3, eigen decomposition and SVD: `tda_mevs` (`mevs`, symmetric
  eigen -- descending eigenvalues, matching `eigen(A, symmetric=TRUE)`
  exactly, including the column-vector convention for `$vectors`),
  `tda_msvd`/`tda_msvd1` (`msvd`/`msvd1` -- singular values alone, or
  the full `A = U diag(d) V'` decomposition, both matching base R's own
  `svd()` convention exactly). `mev`, the general (non-symmetric,
  complex-eigenvalue) sibling of `mevs`, is deliberately deferred -- its
  own dispatch comment in `t_mat.c` understates its arity
  (`mev(A,ER,EI,EV)`, four arguments) against what `m_mev` actually
  parses (`mev(A,ER,EI,EVR,EVI)`, five, split real/imaginary eigenvector
  parts) -- another comment-vs-code mismatch worth resolving properly
  with complex-number handling before wrapping, not folded hastily into
  this batch. `tda_mevs`'s own eigenvector orientation needed a real
  correction mid-batch: an initial 2x2 test instance satisfied the
  reconstruction check (`A %*% v == lambda * v`) equally well whether
  `v` was read as a row or a column of TDA's result, which was wrongly
  taken as confirming TDA uses rows (the opposite of `eigen()`'s own
  convention) -- a 3x3 instance immediately showed rows fail and columns
  succeed; the 2x2 case had only "confirmed" rows because that
  particular matrix's eigenvector matrix happened to be symmetric, an
  accident of the specific numbers chosen, not a property of TDA's
  convention. Caught before gates ran, not after, but recorded here
  since it's a clean example of why a single small test case --
  especially the smallest, most special-cased one someone reaches for
  first -- isn't enough to pin down an orientation convention.
- Batch 4, selecting, permuting, and aggregating: `tda_msrow`/
  `tda_mscol` (`msrow`/`mscol` -- row/column selection, repeats and
  omissions allowed, matching `A[rows,]`/`A[,cols]`), `tda_mprow`/
  `tda_mpcol`/`tda_mpsym` (`mprow`/`mpcol`/`mpsym` -- true permutations
  only, matching `A[p,]`/`A[,p]`/`A[p,p]`; their exact convention is
  stated directly in `m_mperm`'s own header comment, `b(i,j)=a(p(i),j)`
  etc., and confirmed by hand anyway), `tda_mag` (`mag`, sums `A`'s
  entries within each row-group/column-group block -- no base-R
  one-liner equivalent, verified against a hand-rolled double
  aggregate), `tda_mnc` (`mnc`, zero out elements passing a threshold
  test -- already fully documented in `tda.hlp`, confirmed by hand
  anyway since a documented convention can still be transcribed wrong).
- Batch 5, centering/standardizing/miscellaneous: `tda_msqrtd`/
  `tda_msqrti` (`msqrtd`/`msqrti` -- a diagonal matrix built from
  `sqrt(diag(A))`/`1/sqrt(diag(A))`, off-diagonal ignored, same spirit
  as `tda_minvd`), `tda_mcent`/`tda_mstand` (`mcent`/`mstand` --
  column-centering, and centering plus scaling by each column's
  *population* standard deviation, not base R's sample-SD `sd()` --
  confirmed by hand rather than assumed, since this is exactly the kind
  of denominator choice that looks identical to `scale()` until you
  check), `tda_mdcent` (`mdcent`, the classical Torgerson
  double-centering transformation multidimensional scaling applies to a
  matrix of squared distances), `tda_mcross` (`mcross`, `crossprod(A)`),
  `tda_mnum` (`mnum`, an arithmetic sequence -- `seq(x, by=d,
  length.out=n)`). Two commands, already fully wrapped under clearer
  names elsewhere in the package, were confirmed and deliberately
  skipped rather than duplicated: `mbrr` is `tda_brr()`
  (`descriptive.R`) and `mldes` is `tda_mlrc_design()`
  (`further-regression.R`) -- both predate this batch and already have
  real documentation and tests of their own.
- Batch 6, Kronecker product and a backward weighted mean: `tda_mkp`
  (`mkp`, matching base R's own `kronecker(A, B)` exactly),
  `tda_mwvec`/`tda_mwvec1` (`mwvec`/`mwvec1` -- for each position, the
  weighted mean of every later value, `mwvec1` letting an arbitrary
  order vector stand in for array position; no base-R builtin matches
  this, verified against a hand-written R equivalent), `tda_mscal1`
  (`mscal1`, `A / sum(A)`). `mnvar` and `mch` were looked at and set
  aside: `mnvar` turns out to be the primitive `tda_nvar()` (`tda.R`)
  already builds on internally to attach an R data frame as TDA's
  working dataset -- not a matrix operation to expose again; `mch`'s own
  header comment gives a stale argument list (`A(n,m),U(1,m),V(n,1)`)
  against what the code actually requires (`A` square), and the
  row/column scanning logic it implements is not yet fully understood
  well enough to name and document honestly -- deferred alongside `mev`
  rather than guessed at. `tda_mwvec`/`tda_mwvec1` originally used `w`
  for the weight vector; the package's own argument-grammar check
  (`test-api.R`, enforced across every exported function) reserves `w`
  for `tda_cwt_norm` alone, where the whole object literally is the
  weights -- caught by that check on the very first full-suite run,
  fixed by renaming to `weights`, matching the spelled-out convention
  every other TDA wrapper with a weight argument already uses
  (`tda_lsreg`, `tda_sma`, `tda_spl`, ...).
- Batch 7, permutations and edge lists: `tda_mpinv` (`mpinv`, matching
  base R's `order(p)` for a permutation `p`), `tda_mcel` (`mcel`, builds
  a 3-column edge list from an adjacency-style matrix). `mpz` (find a
  row permutation giving a zero-free diagonal, ACM 575) was read,
  wrapped, and then **un-shipped after its own verification failed**:
  even on a trivial permutation matrix with a unique, exact matching,
  the returned permutation did not produce the zero-free diagonal the
  command's own header promises -- traced the exact scatter relationship
  the C code implements (`B[p[i], ] == A[i, ]`, confirmed
  pixel-for-pixel against the raw `tda` output twice) and it is simply
  not enough on its own to reconstruct a zero-free-diagonal `B`, on a
  case where one plainly exists. Rather than ship a wrapper whose
  docstring would repeat a claim just shown false, or invent a corrected
  interpretation without being sure of it, `mpz` was pulled and set
  aside next to `mev` and `mch` for further investigation before it gets
  wrapped.
- **`mpz` resolved and shipped as a genuine C bug fix**, not a wrapper
  workaround: `m_rperm`'s own header (`t_mata.c`) states the permutation
  array's contract precisely -- `A[iperm[i],i] != 0` for every `i`, a
  *gather* index -- and confirming `P` itself against that formula
  (`A[P[1],1]`, `A[P[2],2]`, `A[P[3],3]` all nonzero on both earlier
  test matrices) showed the permutation was computed correctly. The bug
  was one line downstream, in how `m_mpz` used that correct `P` to build
  `B`: a scatter (`B[P[i],:] = A[i,:]`) instead of the gather its own
  algorithm's contract calls for (`B[i,:] = A[P[i],:]`). Fixed in both
  `src/` and `tdaR/src/`; the previously-failing permutation-matrix case
  now returns the identity exactly, and the general case returns a fully
  nonzero diagonal. `tda_mpz()` now ships, verified against the
  `A[p[i],i] != 0` contract directly on random matrices constructed to
  guarantee a solution exists (a random permutation plus random extra
  noise). Pinned in `examples/coverage/mpzops.cf`.
- The two engine fixes above were independently verified in this session
  against the pre-fix code: source reading (dmax in t_gf.c; m_rperm's
  gather contract in t_mata.c), the help entry, and both reproductions
  confirmed byte-for-byte before adoption. A second adversarial control,
  `examples/coverage/mnormpz.cf`, pins both with inputs chosen to defeat
  the accidental blind spots (an all-negative matrix; a 3-cycle
  permutation, since matc2's involution makes scatter and gather
  coincide).
- Batch 8, block-triangularization: `tda_mpbl`/`tda_mpbu` (`mpbl`/
  `mpbu`, Duff & Reid CACM 529 -- the same algorithm family that `mpz`
  had a bug in). Checked with particular care given that precedent:
  `m_perm`'s own header documents its permutation array's contract in
  the opposite direction from `m_rperm`'s (position `i` holds the
  *original* index, a gather already), and `m_mpb`'s output loop
  (`B[i,j] = A[AcM[i],AcM[j]]`) matches that contract exactly --
  confirmed pixel-for-pixel against raw `tda` output on two instances,
  one confirming the `B = A[p,p]` gather relationship on an
  already-triangular scrambled matrix, a second confirming actual
  block-triangular structure (not just index correctness) on a matrix
  with a genuine 2-element mutually-coupled block, checked both for
  `mpbl` (result strictly lower block-triangular) and `mpbu` (upper). No
  bug this time, but the extra scrutiny was the right call to make
  regardless of how it turned out.
- Batch 9: `tda_mpfit` (`mpfit`, iterative proportional fitting -- the
  RAS algorithm for adjusting a table's row and column sums to
  prescribed targets while preserving cross-product ratios). Fully
  documented in its own header (`t_mata.c`); verified anyway against a
  from-scratch RAS implementation rather than trusted from the name
  alone, on tables with randomly generated (but consistent) target
  margins. Built like the optimization-family wrappers (writes its five
  differently-shaped inputs to files and calls `tda_run` directly)
  rather than through the generic `tda_mat()`, since it needs to parse
  extra text output (iteration count, final accuracy) that `tda_mat()`'s
  matrix-only return doesn't expose.
- Batch 10: `tda_mpit`/`tda_mpit1` (`mpit`/`mpit1`, Leslie matrix
  population projection -- fully documented in its own header, verified
  against a from-scratch recursion), `tda_mkmet` (`mkmet`,
  Kemeny/Kendall-tau-with-ties distance between rankings, verified by
  hand), `tda_mple` (`mple`, the product-limit/Kaplan-Meier CDF via
  Efron's redistribute-to-the-right algorithm -- TDA's own lower-level
  building block behind `tda_ple`, confirmed to be a genuinely separate
  command, not a duplicate, before wrapping it). `mple` surfaced two
  real, non-bug algorithmic properties worth documenting precisely
  rather than glossed over: tied *events* (not censoring ties, which
  work as expected) are processed one at a time rather than as a single
  simultaneous risk-set reduction, so only the last-processed
  observation in a tied-event group reaches the value standard
  Kaplan-Meier software reports for that time; and `F` always reaches
  exactly 1 at the highest-time observation, even when it is censored, a
  direct consequence of the algorithm only ever moving probability mass
  forward and never destroying it -- standard software instead leaves
  the tail undefined beyond a censored maximum. Both confirmed by hand
  on small examples and pinned as explicit test cases, not folded into
  the random comparison against `survival::survfit` (which was itself
  adjusted to avoid the two edge cases it can't meaningfully check).
  `midf`/`midf1`/`midf2`/`midf3` (`t_imat.c`, a whole separate
  interval-arithmetic subsystem for interval-censored distribution
  estimation) were read far enough to see they need real, dedicated
  investigation -- deferred alongside `mev` and `mch` rather than
  rushed.
- **`mev` resolved and shipped**, the first of the three previously
  deferred commands: general (possibly non-symmetric, possibly
  complex-eigenvalue) eigen decomposition. Its true arity (five
  arguments, `mev(A,ER,EI,EVR,EVI)`, not the four the dispatch comment
  lists) was already established when it was deferred back in batch 3;
  what remained was pinning the eigenvector convention. A first
  verification attempt seemed to show *left* eigenvectors (`A' %*% v ==
  lambda * v`) rather than the right eigenvectors base R's `eigen()`
  returns -- but that conclusion came from a matrix constructed with a
  `byrow` mistake in the R-side verification code itself, not a real
  property of TDA's output. Rebuilding both test matrices correctly
  showed standard right eigenvectors (`A %*% v == lambda * v`)
  throughout, matching `eigen()`'s own convention exactly, confirmed on
  a real-eigenvalue case and a genuinely asymmetric complex-eigenvalue
  case (not simply plus/minus of each other, unlike the first,
  too-simple matrix that had initially masked the construction bug).
  `tda_mev()` returns eigenvalues and eigenvectors as R's native
  `complex` type.
- **`mch` resolved and shipped**, the last of the three deferred
  commands. Its header names only the vaguest description ("Check row
  and column sum conditions") with no formula and no textbook name to
  match against, so the exact rule was pinned down by tracing the loop
  body directly: for each index, an outward running-sum comparison
  against its mirror-image entries, flagging the index the first time
  the running total falls behind by more than a small tolerance. No
  independent formula exists to check this against, so the verification
  *is* a from-scratch line-for-line transcription of the same loops into
  R, cross-checked against TDA's own output on three hand-picked
  instances (an already-consistent upper-triangular matrix producing no
  flags at all, and two matrices with a deliberately out-of-order entry
  producing the expected pair of flags) before being trusted for the
  randomized test suite. `tda_mch()` documents the mechanical rule
  precisely rather than attach a name or interpretation this project
  can't back up. This closes out every command flagged for deferral this
  session.
- **The `midf` interval-arithmetic family, resolved and shipped**,
  closing out the one deferral (from batch 10) still open after
  `mev`/`mch`. All four commands (`midf`, `midf1`, `midf2`, `midf3`)
  decode cleanly once read carefully, no guessing needed: `midf`
  estimates a distribution function for interval-censored data (each
  observation known only to lie in `[lower,upper]`) at every distinct
  interval endpoint, returning a sure lower bound, a sure upper bound,
  and a uniform-within-interval point estimate -- its own header lists
  only three outputs, `midf(XL,XU,DL,DU,DM)`, but `m_midf` actually
  returns four (the breakpoints themselves are a result too), the same
  understated-arity pattern `mev`'s header had. `midf1` is the same
  point estimate evaluated at each observation's own lower endpoint
  instead of at every breakpoint. `midf2` computes each observation's
  restricted mean (the expected value given the true value exceeds that
  observation's own lower endpoint), built by integrating a
  piecewise-uniform density derived from `midf`'s own CDF estimate.
  `midf3` is unrelated to distribution estimation entirely -- a local
  smoothing that, for each interval, averages the lower (and separately,
  upper) endpoints of every interval whose corresponding endpoint falls
  inside it. All four confirmed cell- for-cell by hand on the same
  3-interval example before writing any wrapper, then verified in bulk
  against from-scratch R translations of the exact formulas (no
  independent textbook name applies to `midf2`/`midf3`, so the
  translation is the specification, the same situation `mag`/`mch` were
  in). The randomized test data was switched from rounded decimals to
  integers after a first attempt hit floating-point tie-breaking noise
  at exact-equality breakpoints introduced by the text round-trip
  through TDA's own formatted output -- a test-robustness fix, not a
  semantic one.
- `tda_qsort_r()` (`tda_compat.h`, `t_sort.c`): the index sorts went to
  the C library -- glibc `qsort_r`, BSD `qsort_r`, Windows `qsort_s`.
  glibc's is a stable merge sort and every reference output was made
  with it; the other two are not stable, so tied values came out in a
  different order and e.g. `mple` reported a different intermediate F
  for tied events on Windows. Replaced by a stable merge sort in
  `t_sort.c`; every platform now agrees with the references.
- `sdshp()` (`t_sdx.c`): `.shx`, `.shp` and `.dbf` were opened in text
  mode, which corrupts them on Windows ("can't read header of dbf
  file"). Opened binary, as the other binary readers already were.
- `rcsv()` (`t_sdx.c`): a CR before LF was kept as part of the last
  field, so CSV files written on Windows never round-tripped, on any
  platform. CR is now skipped; LF and CRLF files parse identically.
- `t_tri.c` free lists: recycled `Site`/`Edge`/`Halfedge` nodes were
  written through `struct Freenode *`, a strict-aliasing violation the
  build had been carrying `-fno-strict-aliasing` for. The link is now
  stored with `memcpy`; the flag is gone. No output change.
- `t_cl.c`, `t_freq.c`, `t_qp.c`: three locals (`qic`/`qjc`, `w`, `j`)
  initialised at declaration. Every read was already guarded by the
  condition that assigns them; the initialisers only quiet Rtools'
  `-Wmaybe-uninitialized`. No output change.

# Audit residue: what TDA prints that no channel carries

`tests/testit/test-audit-cf.R` runs all 634 command files of the C
suites through `tda_run()` and classifies every printed line. After the
echo/setting/notice/diagnostic/debug patterns and the export channel are
accounted for, these printouts remain -- results that exist only in
`$run$output` unless a wrapper parses them:

- `bfc`, `bfa` (Boolean context dependencies, covers, prime implicants):
  text results.
- `clu`, `spmod`, `dem`, `rfit`/`rfit1` (properties list, rank table),
  minmax `mds` configurations and trace: numeric tables.
- `ctab`: the small-cell counts under "Contingency measures".
- `rcsv`, `rdbf`: reader statistics (records, record length).
- `dump`: the hex dump; `pdata` with no `df=`: the data to stdout;
  `sga`: the dot text.

The wrappers that exist (`tda_bfc`, `tda_spmod`, `tda_rfit`, `tda_dump`,
`tda_sga`, `tda_rcsv`) parse the text, so their objects are complete;
the gap is the export channel, and the commands without a wrapper. Line
counts are pinned in the test; a rise means TDA started printing
something new.

- Leftover developer prints removed from the C (`t_mds.c` mdsn/mdsn1/
  mdsxr/unf, `t_sdr.c` sdcpol/sdrel/sdp lists, `t_sd.c`, `t_bfa.c`,
  `t_cl.c` clpyr, `t_clip.c`, `t_map.c` geodesic/projection): dumps of
  intermediate configurations, gradients, distance vectors, edge lists
  and German working notes (`HIER nach step`, `current XY nach
  startwert`, `Pointer und NE`, `NEUE DEGREES`, `Kantenliste`, ...).
  Output only; no computation touched (the gradient-norm loop in mdsn1
  kept its sum, lost its print). None of it appears in TDA's own
  `exam`/`ehhnew` references; 13 generated references in `tests/` and
  `coverage/` were regenerated, and each differs from before by removed
  lines only.  Same treatment as `mreg`/`conj` earlier.
- Package-side producers so that no output file is asked for or read
  (`tdaR/src` only, standalone untouched): `export_prn()` in `t_pgen.c`
  hands the matrix `prn_data()` would print to R as `lsreg.vcov`,
  `glm.vcov`, `l1reg.vcov` (`t_lsreg.c`, `t_glm.c`, `t_l1reg.c`) and
  `loglin.vcov`, `loglin.vcov.2`, ... (`t_loglin.c`), whether or not
  `pcov=` was given; `t_min.c`'s covariance block is entered
  unconditionally under the package build so `ml.vcov` is always
  produced. The file-write tap (`tda_export.c`) now also hands over each
  stream's file name (`file.<stream>.name`), its `#` comment lines
  (`file.<stream>.comments`), and integer literals written without an
  argument (the `1 3 4` object header of a spatial file), and skips
  numbers inside comment lines. tdaR asks for no `pcov=` anywhere and
  `tda_file()` builds every table from the tap by file name; the
  wrappers for TDA's file utilities (`esort`, `eskip`, `eselect`,
  `emerge`) are file-in, file-out and read nothing back.
- `rdataframe` (`tdaR/src/tda_rdf.c`, package only) did not set what
  `nvar` sets when it creates a data matrix: the case-weight state
  `WIVar = -1`, `WSum = NOC`, `WSumS = 0`, `WNorm = 1`, `WNormFlag = 0`.
  `WSum` stayed at zero, so every statistic normalised by the weight sum
  divided by zero -- `ple(csf)`'s Wilcoxon (Breslow) and (Tarone-Ware)
  came out `-nan` for any data frame handed over from R, while the same
  data read through `nvar(dfile=)` gave the manual's 2.6510 and 2.9767
  (6.5.4, `ple6.cf`). Found by comparing the vignette's box with the
  manual's; the standalone was never affected.
- `tda_rate(tp=)` (and `prate(tab=)`) wrote every vector of time points
  as TDA's `a (step) b` shorthand, which is only right for an equally
  spaced grid: `pl7.cf`'s `tp = 0,170.1,354.1,535.1` became `0 (170.1)
  535.1`, i.e. the periods 0, 170.1, 340.2, 510.3, and a different fit
  (goodness of fit 12.69 on 4 df instead of the manual's 10.2351 on 3).
  Unequal points are now written out in full. Found by comparing the
  vignette's box with the manual's.
- `gflow` with `opt=1` (`t_gio.c`): the per-source-node prefix (index,
  node) that belongs to the matrix form `opt=2` was written for both
  forms, so the first pair record of each source node carried two extra
  leading numbers and a node without records left them dangling. The
  manual's Box 1 of 7.2.8.1 shows the six-column record without them;
  the prefix is now written for `opt=2` only. Console output is
  unchanged, so no `.ref` moved.

## dplot laid out its grid from each input's own bounding box (t_dplot.c)

**The symptom.** A grid of panels combined with `dplot` came out with
its columns out of line: the second panel of one row started at a
different x from the second panel of the row above, and the panels were
not the same size.

**What caused it.** `dplot` sized each cell from that input file's own
`%%BoundingBox` -- `xorg1 += x2 - x1` per panel, `bx` the widest row.
Two inputs of the same physical size can still declare different boxes,
because the box has to reserve room for anything drawn outside the plot
frame, and a `pltext` label is drawn from its anchor without regard to
the frame. A panel carrying a longer label therefore declares a wider
box than its neighbour, takes a wider cell, and shifts everything after
it in its row.

**It is visible in TDA's own documented example.**
`examples/exam/plot7.cf` is the manual's illustration of `dplot`
(section 4.3, Figure 1). Its inputs `plot3.ps` and `plot4.ps` declare
315 and 326 points of width, so the two rows disagree: rendered at 150
dpi, row 1's second panel starts at x = 298 px and row 2's at x = 308.
The one example the feature ships with does not produce a square grid.

**The fix.** Read every input's box first, take the largest width and
height as the cell, and place panel *j* of each row at `j * cw`. Each
panel keeps its own offset within its cell (`x - x1`), so the frames
line up. Inputs of equal size give exactly the layout as before, which
is why `plot7.ref`'s own `New bounding box` and scale factors are
unchanged -- the only new line in it is the `Cell size` diagnostic.
Applied to both copies, `src/t_dplot.c` and `tdaR/src/t_dplot.c`,
confirmed identical after editing.

**Control.** `examples/tests/dplotbb1.cf`: two panels differing only in
the length of their `pltext` label, combined into two rows with the pair
in opposite order, so a per-file layout puts the second column in two
different places. Verified against a binary built from the pre-patch
source: before, cells of 283 and 347 points and the columns out of line;
after, a 283-point cell throughout. C suite 494/494.

Not patched, and worth recording so it is not re-raised: the bounding
box itself is correct. Measured against real ink with ghostscript, TDA's
per-character allowance is the width of the widest glyph -- a
48-character string of `W` declares 511.8 points and reaches 511.4 -- so
the box is a tight upper bound, exact in the worst case and generous for
ordinary text. Making it exact would need Adobe width metrics for
Times-Roman and Symbol and would change the box of every shipped `.ps`
carrying text.

## get_zoo read only the 13-byte short name of an archive member (t_zoo.c)

`b_to_dir()` copied the fixed part of a directory entry and stopped at
the DOS-style short name, so a member stored by zoo 2.1 or later under
a long file name, or inside a directory, could be named in an archive
description file only by its 8.3 form (`a_much_l.txt` for
`a_much_longer_file_name_than_dos_allows.txt`) and never by its
directory. `b_to_dir()` now also reads a type-2 entry's variable part --
the long name and the directory name, bounds-checked against the
declared length -- and `get_zoo()` matches a description-file name
against the full stored path (`data/sub/name.dat`) and, as before,
against the bare name. Type-1 archives and every TDA-era archive are
unaffected. Control: `examples/tests/arcd_long.cf` on `long.zoo`, an
archive written by the reference zoo tool (github.com/troglobit/zoo)
from TDA's own example data under long names in a subdirectory; `arcc`,
`arcv` and `nvar(... = A:V1)` all resolve the members. The package's
`zoo()` writes the same type-2 entries (the 8.3 short name derived as
zoo's `dosname()` derives it, the entry's own CRC-16) and `unzoo()`
reads them; `tools/check_zoo_ref.R` compares both directions and every
directory-entry name field with the reference tool, and
`test-zoo.R` reads `long.zoo` from `inst/extdata`.

## logLik of an unconverged fit was read from six-digit text (R/tda.R)

`logLik()` on a run read the printed "Log likelihood (final
estimates)" line, which TDA prints at the run's `tfmt` (full precision
in practice), with "Maximum of log likelihood" as the fallback. The
first line is printed only after convergence; the fallback is `%lg`,
six significant digits. A fit that stops at the iteration cap or
misses a tight tolerance therefore carried a rounded log-likelihood.
Found by r-universe's aarch64 build, where `tda_strict()`'s 1e-12
tolerance is not reached on the `tda_control` example (the fused
multiply-add arithmetic lands a few ulps away from x86) and the
example-audit test flagged the six-digit value. `logLik()` now takes
the `ml.logLik` export, which the C emits beside the six-digit line
whether or not the iteration converged; the text is the fallback for
the exporters-off path only. Pinned in `test-vignette-findings.R`
with a fit truncated at one iteration.

## tda_spl(rx=) named its columns one place out (R/smoothing.R)

`spl` reports the observed y beside the fitted value when it fits at the
data points, and omits it when `rx=` puts it on an interpolation grid,
where there is no observation to report -- six columns instead of seven.
The reader named them from a fixed seven-name vector either way, so
under `rx=` every column after `x` shifted: the fitted values were read
as `y`, the first derivative as the fitted values, and each derivative
as the one below it.

Found through the book: Figure 9.4.2 is the first derivative of a cubic
spline through the cumulative baseline rate (`ehi12.cf` smooths column 7
of `ehi10.bl` over `rx = 0(0.5)300`, `ehi13.cf` plots column 4 of the
result). Plotting `d1` gave a curve falling to zero by t = 20 instead of
the book's hump; the column actually being plotted was the second
derivative. Fixed by taking the names from the table's own width. With
`rx = seq(0, 300, 0.5)` the first row now reads `x 0, fitted
-0.02648526, d1 0.01191835, d2 0.001294091`, against `ehi12.bl`'s own
`0.0000 -0.0265 0.0119 0.0013`.

## prate= took one covariate constellation, TDA takes several (R/rate.R)

`ehi10.cf` asks for two baseline rate calculations in one `rate`
command, cohort 1 and cohort 3, both written to `ehi10.bl` and told
apart by the
sub-table ID in column 1; `ehi11.cf` then draws both curves. The wrapper
accepted a single named list, so only one constellation could be
requested and Figure 9.4.1 came out with one curve where the book has
two.  `prate` now also takes a list of named lists, one `prate()` block
each. Found by reading the command file against the figure, not by any
check that existed at the time -- `tools/check_vignette_plots.py` was
written afterwards and now compares every vignette plot's axes, tick
spacing, series count and labels against its own `.cf`.

## the episode table dropped the per-spell totals TDA prints (R/rate.R)

Multi-episode data prints one block per spell, each closed by a `Sum`
line carrying that spell's episode and weighted totals. The parser
skipped those lines, so the totals TDA printed were stored nowhere --
the components were there, the totals were not. They are now kept as the
`sums` attribute of `$episodes`, not as rows, because several tests
and wrappers index that frame by position. Found by comparing Box 2.2.8
against the vignette: the book prints 201, 162, 107, 62, 32, 20 for the
first six spells, which are the sums of the pairs we already had.

## print.tda_ple was defined twice, and the wrong one won (R/rate.R)

`rate-methods.R` carries the real method: it paginates with `n=`, names
each block's transition and group, and prints the median duration and
case counts under it. `rate.R` then had `print.tda_ple <-
print.tda_ltb`. R collates `rate-methods.R` before `rate.R`, so the
alias overwrote the method: every product-limit result printed through
the life-table
printer, which caps at ten rows and ignores `n=`. The real method had
never run.

Found while trying to print a product-limit table in full for the
vignette: `print(pl, n = nrow(pl$blocks[[1]]))` kept reporting "... 121
more rows". Alias removed; `n = 131` now prints 131 rows, and the book
comparison moved 7329 -> 7374 as the tables the book prints in full
became visible.

## NOT A FIX: the Cox goodness-of-fit table was already stored

Recorded so it is not "found" a third time.  A Cox model given `tp=`
prints a global goodness-of-fit table under the estimates. Comparing Box
9.3.5 against the vignette suggested nothing stored it, and a second
parser was written for it. That was wrong: `$gof` already carried the
table, and `summary()` already printed it -- the vignette simply never
called either. The duplicate parser has been removed and the vignette
uses `summary(ehi8)`.

What exposed the duplication was the structural pass in
`tools/check_vignette_audit.R`: it reported the block as covered after
the new parser was deleted, and the frame covering it turned out to be
`$gof`. The lesson is the cheaper check first -- grep the accessors for
an existing field before writing a parser for output that looks
unhandled.

## tda_spl()'s two column layouts cannot be told apart by width

A correction to the entry above, which got this wrong. `spl` writes six
columns either way, but they are not the same six:

    without rx:  index  x  y  fitted  d1  d2
    with rx:     index  x     fitted  d1  d2  d3

The grid has no observation to report, so `y` goes and `d3` takes the
free column. Keying the names on the table's width therefore reads the
observation as the fitted value and shifts every derivative -- which is
what the first attempt did, after a noiseless test that could not tell
`y` from `fitted` because the two coincide there. The names now key on
whether `rx` was given.

Caught on ds5.cf: it plots column 4 of its spl output and TDA draws a
smooth curve, while the vignette drew a jagged one through the noise.
Our column 3 matched the input Y exactly, which settled it. Verified
both ways: without `rx`, `y[1]` is the observation -0.0885 and
`fitted[1]` the smooth 0.0156, matching TDA's own df; with `rx`,
`fitted[1]` is -0.026485 and `d1[1]` 0.011918, matching ehi12.bl and
Figure 9.4.2.

## a plot with no clip region was fitted to its content, not its frame

`tda_plot_ps()` takes the view from TDA's clipping region where there is
one, and otherwise from the extent of what was drawn. A graph drawing
sets no clip, so a plot whose content sits in part of its declared
coordinate range was blown up to fill the frame. `gd31.cf` declares 0-8
by -1-5 and puts its four nodes in the middle third: TDA draws them
small and central, we drew them filling the picture with the node labels
correspondingly huge.

The file's own `%%BoundingBox` is TDA's statement of the area it drew
into, whitespace included, so it is now used as the view when there is
no
clip. Verified against `gd31.ps` rendered by ghostscript at the same
pixel size, and against `gd23.ps`, which was already correct because its
content nearly fills its range -- which is why this went unnoticed.

## OURS: a subscript increment lost in the port (t_gdf.c, gdf_edf3m)

Rohwer averages the exact observations in each box of method 1's grid,
one mean per dimension:

    for (j = 0; j < ndim; ++j)
        AcYF[k++] /= (double)AcS[i];

Rewriting the compound assignment for the explicit casts of the
TDAContext refactor moved the `k++` out of the subscript and past the
loop body, so `AcYF[k]` was divided ndim times and `AcYF[k+1]` never:
with two dimensions the first box mean came out over the count squared
and the second stayed a sum. Every expected value leaning on such a box
was then far too large.

Caught against the manual, which was right all along. Section 6.2.2 Box
12 tabulates the two joint methods over a hundred generated units;
method 2 matched and method 1 did not. With the increment restored,
`gdf4.cf` gives 37.6603, 37.6603 and 46.2104 for the rows whose first
dimension is 34, 36 and 43 -- the manual's own figures, digit for digit.
The protocol file shows the mechanism directly: for a box holding two
exact observations at 80 and 82 the means read 40.5000 and 162.0000
before, 81.0000 and 81.0000 after.

I first reported this as a fault in TDA itself, on the evidence of the
working tree's own copy -- which already contained the regression.  It
is ours. Checking the whole file for the same class afterwards: of the
six `array[idx++]` uses in Rohwer's t_gdf.c, five survived the port and
only this one lost its increment; the single `array[++idx]` survived.
That check wants running against the pristine sources for every other
file the refactor touched.

## NOT A DISCREPANCY: section 6.2.2 Box 12 (superseded)

Recorded and then withdrawn.  The manual's method 1 column was reported
here as disagreeing with the program; it does not.  The disagreement was
a regression of ours in `gdf_edf3m`, described in the entry above. The
manual is right and the program now agrees with it.

## the last observation of a product-limit table was dropped

TDA writes a transition's final line with a leading "#" when only
censoring happened at that time:

    #  0  131  428.00        0        8

and writes nothing after the two counts -- no risk set, no survivor, no
standard error, no cumulated rate.  There is no estimate to write: the
survivor function does not change where nothing fails. The file reader
skips "#" lines and the ple.table export never carried them, so the row
was lost, and with it the number of observations still running at the
end of observation.

It is now a row of `$blocks` with NA in the four columns TDA does not
compute, which is what the file says. `tda_survivor()` drops it -- it
returns the estimate, and there is no estimate there -- so plotting and
the confidence bands are unaffected. 428.00 with 8 censored for the
single-transition fit; 13, 8 and 10 at 428.00 for the three transitions
of the four-state data.

It was first kept on an attribute, on the strength of a comment in
t_ple.c saying that making it a row "put NA rows into $blocks, $table
and tda_survivor(), which broke plotting outright". That is a reason to
make the readers NA-aware, not to hide a row of the data: an attribute
is not somewhere anyone looks. Two tests compared $blocks against a
file-derived reference and needed the estimate rows selected; both now
say so.

## gcyc's tables came back unnamed (R/graph.R)

`tda_g_cycles()` returned V1, V2, ... for options 2 and 3, so finding
the edge in a row meant counting columns. Every option writes the cycle
number and an index first; option 2 then gives the edge and one
indicator per fundamental cycle, option 3 the same over every cycle
found, option 4 the edge alone. Named `cycle`, `index`, `from`, `to`,
`in1`... to match the manual's own headers (7.2.6.1, Boxes 1 to 4).

Option 1 (the nodes of each cycle) and option 4 stay ragged, correctly:
option 1's rows are genuinely of different lengths, and option 4 writes
two different row shapes into one file -- eleven edges, then one
indicator row per cycle. A rectangular ragged result is now converted to
a data frame, which leaves both of those alone.

## the estimated correlation matrix was never asked for (R/regression.R)

For the multivariate and simultaneous probits (TDA's QRPROB3 and
QRPROB4) `qr_corr()` writes the estimated correlation matrix among the
latent equations. It does so only when `ppar=` names a file, and it
writes the rows commented so that reading that file back as starting
values skips them. The wrapper never set `ppar=`, so TDA never computed
the matrix at all -- not dropped on the way out, never requested.

`tda_qreg()` now sets it for models 6 and 8 and returns the matrix as
`$correlation`.  Checked against TDA's own ppar file for qr7.cf, which
has 0.4134 off the diagonal: ours matches to the printed digits. The
vignette's qr6 box is the identity, correctly -- its own constraints
zero those terms.

Found by asking which other commands write commented rows into their
output files, after the product-limit case above. `ple` and this were
the only two; the other commented numeric lines in the shipped files are
input data where "#" is the author's own comment.

## gbcf's summary file was never asked for, and its table was unnamed

`gbcf` writes two files: the right-hand side takes the pair table -- one
row per controlling/controlled pair -- and `df=` a summary with one row
per node.  Without `df=` the summary is not written at all, so
`tda_g_backward()` never had it. It now asks by default and returns it
as `$summary`; both tables are named after the manual's own headers
(7.6.1.1, Box 4), where the pair table was coming back as V1 to V8.
Checked against the manual: the summary's first rows are 1/1/8, 10/10/1
and 11/11/0.

Noted while there: `becl`'s options 3 and 4 fill a non-edge with -1
unless `sc=` says otherwise, and `cl5.cf` passes `sc = 0`. The vignette
was calling them without it and printing a matrix of -1 where the manual
prints 0. Not a wrapper fault -- a parameter the vignette had not passed
-- but worth knowing that the default is not the manual's.

## secondary output files: a sweep, and three that were unreachable

TDA writes several of its results only when the caller names a file for
them -- `ppar=`, `pres=`, `pcov=`, `ptab=`, `df=`, `prot=`. If the
wrapper never sets one, the result is not written; it is not dropped in
transit, it is never computed. Swept every `PM<X>FDef` guard in the C
against the R side:

  pcov   covered: $vcov comes through the export channel instead, and
         matches TDA's own pcov file to the last digit (checked on
         lsreg: 275.5435943890666977 and the rest)
  ppar   was missing for the multivariate and simultaneous probits --
         fixed above, $correlation
  df     was missing for gbcf -- fixed above, $summary
  pres   covered for lsreg, loglin, l1reg, rate (residuals=) and glm
         (predictions=, the same file under another name); MISSING for
         nlreg and lsreg1, both now given residuals=
  ptab,
  prot   set where the wrappers need them

gmin also writes a pres= table, one row per case over the data matrix,
but it is a function minimiser that need not have a data matrix at all;
left alone deliberately.

## dgrp's per-group table was printed and not kept (R/regression.R)

With `dgroup=` TDA prints one row per group -- the indicator, the number
of cases in it and its weight -- and stores it nowhere. It is now
`$groups` on the fit (manual 6.9.1.3, Box 2).

The weight prints at four decimals, and the package's own audit caught
the first version of this reading 0.2500 from the text when the print
tap has the value in full. Taken from the tap, with the text kept only
as a fallback.

## gmin/gmax printed six counters and kept none (R/smoothing.R)

Above its accepted boxes `gmin` reports the iterations performed, the
function and inclusion-function evaluations, the boxes used, and the
temporarily and finally accepted counts (manual 8.4.1, Box 2).  None was
stored.  Now `$counts` on the result.

Found by classifying what the manual-comparison still missed rather than
by reading boxes one at a time -- the same residual walk the book got.

## hcld's second output file was never asked for (R/clustering.R)

`hcld` writes a second file with `df=`: one row per split, the two
clusters, their sizes and their diameters. Without `df=` it is not
written at all, and `tda_cluster()` never set it. Now returned as
`$merges`, named after the manual's own headers (7.5.2.1 Box 3 and
7.5.2.2 Box 2): 2/1/6/11/33/55 for the first split of the second
algorithm, matching the manual.

The third instance of this shape today, after qr_corr's ppar= and gbcf's
df=. The sweep in the entry above found them; this one was missed by it
because tda_cluster reaches hcld through .dist_cmd(), so the option
never appears in clustering.R's own text.

## seqpd had no wrapper (R/sequence.R)

Every sequence command had an R wrapper except `seqpd`, TDA's own "print
sequence data" -- the largest single gap the manual comparison found. It
writes the sequence data back out in one of four layouts (`m=`): every
structure side by side on one row per case; one row per case and
structure; one row per case and position; one row per spell with its
origin and destination state and the times it spans.

`tda_seqpd()` follows the same shape as its siblings, through the
`.seq_desc()` helper they share. Checked against the manual's own
section 3.4.3: Boxes 2 and 3 now read 240/240 and 316/316. The manual's
example defines two sequence structures and prints two sets of columns;
seq.d1 as the package ships it has one.

### tda_seqpd: the argument name, the example, the tests

Written first as `structure = 1`, which shadows `base::structure` and
tells the reader nothing. It is `layout = 1`, documented as a \describe
list of what each of the four produces, and a value outside 1:4 is
refused rather than passed to TDA to fail on.

The example was `\dontrun{}` around `ex("seq.d1")` -- neither runnable
nor checked. It is now four cases inline, like the rest of the family's
examples, and it runs: extracted with Rd2ex and executed against the
installed package.

Tests in tests/testit/test-sequence.R, against the manual's own section
3.4.3: case 1's sequence is 7 3 3 1 3 3 7 7, and as spells 7->3 over
0-1, 3->1 over 1-3, 1->3 over 3-4, 3->7 over 4-6, 7->-1 over 6-8. Six
assertions, including that a bad layout errors.

## everything added today, documented and tested

Prompted by the seqpd review: the new return fields had been added and
none of them documented as a return value or covered by a test. The word
appearing somewhere in an unrelated .Rd is not documentation.

Documented on their own pages now: tda_qreg's $correlation,
tda_cluster's $merges, tda_lsreg's $groups, tda_g_backward's $summary
(on tda_g_analyses, where that function is actually documented -- the
first attempt put it on the wrong block), tda_nlreg/tda_lsreg1's
residuals=.

Twelve assertions in tests/testit/test-vignette-findings.R, every value
from the manual's own boxes rather than from our output: gbcf's node
summary (1/10/11, 8/1/0), hcld's first split (2 1 6 11 33 55), qreg's
latent correlation (0.4134 off the diagonal), lsreg's four groups of
twelve -- including that the weight is full precision and not the
printed 0.2500 -- nlreg's residuals present with the option and absent
without it, and ple's censoring-only last row being in the table while
tda_survivor() leaves it out.

## only half of TDA's sequence data model was reachable (R/sequence.R)

TDA reads sequence data two ways. Type 1 is one column per time point.
Type 2 -- `seqdef`'s own `m=2` -- reads the columns as state/time pairs,
one pair per spell, and builds the time axis from the times rather than
from the column count. Every wrapper in the family emitted `seqdef =
cols;` and nothing ever emitted `m=2`, so type 2 could not be declared
from R at all.

Demonstrated against the binary before fixing: the same six columns of
seq.d2 read as type 1 give six variables over 0-5 with nine states, and
as type 2 three variables over 0-8 with three states, which is what the
manual's section 3.4.2 Box 6 prints.

`type = 1` now, on tda_seq_info and the descriptive family. Named for
the manual's own "Sequence type" -- the first attempt called it `pairs`,
which collides with seqm's own pair table in the same file.

Four assertions in tests/testit/test-sequence.R against Box 6's numbers.

## OURS: roxygenise dropped print.tda_fit, and objects spilled their run log

Regenerating the docs for tda_seqpd re-ran roxygen over the whole
package, which rewrites NAMESPACE. `print.tda_fit` carried no roxygen
tag -- its registration was a hand-written NAMESPACE line -- so it was
dropped: 46 S3method entries became 45.

Everything of class tda_fit then fell through to print.default, which
prints every element of the list including `$run`, the console log of
the TDA process. The vignette's section 6.12.2.1 showed the TDA banner,
the variable table and the memory figures where the coefficients should
have been. tda_fml returns c("tda_fml", "tda_fit") and has no print
method of its own, so it was the visible casualty.

Tagged with @exportS3Method base::print and regenerated; 46 again.
Checked against the archive's own NAMESPACE, which has the same 46.

Tests: every returned class has a registered print method, and printing
a tda_fml does not spill the run log. The vignette is swept for the
banner, "Reading command file:", "Current memory:" and "$run" -- zero
boxes.

## ejoin dropped every covariate (R/episodes.R)

`ejoin` carries every column past the six structural ones straight
through to the joined result, filling -3 where that side has no spell --
the manual's own example joins X1 and X2 from one file with Y1 from the
other (section 3.3.5, Box 2). The wrapper built its input files from the
six it computes and wrote nothing else, so the covariates never reached
TDA at all.

The extra numeric columns now go into the input files and their names
are carried onto the result. Box 2 reads 181/181, from 170.

Found by looking only at the distinctive misses -- values of three or
more significant digits, which cannot match by coincidence. All eleven
of that box's misses were the X2 column: 22.2, 33.3, 44.4, 55.5, 66.6.

## 6.12.5 Box 2: the standard errors differ in the fourth decimal

The multivariate probit's coefficients match the manual to every printed
digit -- 2.8003, -0.1915, 3.3703, -0.2320, 4.8489, -0.3740, 4.5686,
-0.3867. The standard errors do not: 0.5528 against the manual's 0.5530,
0.5454 against 0.5456, and so on, with the C/Error column following.

Running qr6.cf through the binary gives 0.5528, so the wrapper agrees
with the program. Whether the program agrees with what produced the
manual is a separate question and is NOT settled here: the gdf case
earlier in this file was exactly such a difference and turned out to be
ours. What can be said is that the two sweeps against Rohwer's own
sources (subscript side effects, control flow) report nothing in the
covariance path -- t_qrmod.c's only control-flow difference is the
removed unconditional goto documented above, which is in rmod, not here.

Recorded as open rather than dismissed.

## the manual's 6.12.5 standard errors: settled, not ours

Built Rohwer's own unmodified sources (tools/build_pristine.sh) and ran
qr6.cf through both binaries:

  Rohwer's own:  2.8003  0.5528  5.0655
  ours:          2.8003  0.5528  5.0655
  the manual:    2.8003  0.5530  5.0641

Identical to ours, digit for digit. The manual's figures come from
neither binary; it is out of date against the program Rohwer shipped.

Two more settled the same way in the same session:

  ple.2's "Duration times limited to" -- Rohwer's binary says 350, as
  ours does, not the manual's 428. the gdf joint method-1 fix --
  Rohwer's binary gives 37.6603, 37.6603, 46.2104, exactly what the
  patched build produces, so the fix restored his behaviour rather than
  inventing one.

This is the tool that should have existed from the start. "Our binary
says so" is not evidence about whether a difference is ours; the gdf
case proves it, where that reasoning would have closed a real bug as
someone else's problem.

### hcld stops after two splits by default

The manual's own run continues until the size bound; `max_splits=` is
how far it goes. The vignette's box was showing two splits against the
manual's nine rows (two rows per split, plus one where its run stopped
mid-split). 7.5.2.2's boxes read 81/82 and 71/72 with max_splits = 5.

## checked against R's own estimators, not only against TDA

Agreeing with the manual or with a .ref file shows the port reproduces
TDA. It does not show TDA is right. Where base R fits the identical
model there is an independent check, now in tests:

  logit   coefficients AND standard errors match stats::glm to five
          decimals (-3.22464, 5.97018; 0.88606, 1.44917)
  lsreg   coefficients match stats::lm, and $vcov matches vcov()
  probit  coefficients match to four decimals; the STANDARD ERRORS do
          not match glm's, and should not

The probit difference is real and explained, not waved away. TDA inverts
the observed information -- the Newton-Raphson Hessian -- where glm()
uses the expected information from Fisher scoring. For the logit the
link is canonical and the two coincide, which is why that one agrees
exactly. Computing the observed information in R gives 0.4493 and
0.7455, which is what tdaR reports and what Rohwer's own binary reports.
The test asserts against that, not against glm's.

Found only because a first attempt at this comparison used the wrong glm
specification -- Response is 0/1 with Weight a frequency, not a binomial
denominator -- which made the coefficients disagree wildly while the
standard errors happened to agree. Worth remembering: a comparison that
fails is as likely to be the harness as the code.

## boxes where the manual is out of date, settled against Rohwer's binary

Four now, each checked by running the manual's own command file through
a build of Rohwer's unmodified sources and through ours, and getting the
same answer from both:

  6.12.5 Box 2   standard errors 0.5528 / 5.0655, manual 0.5530 / 5.0641
  6.5.2  Box 5   "Duration times limited to: 350", manual 428
  8.4.1  Box 2   19 iterations, 20 evaluations, 35 inclusion evaluations;
                 manual 39, 41, 77
  (and the gdf joint method-1 case above, where the manual was RIGHT and
  the difference was a regression of ours -- which is why each of these
  is checked rather than assumed)

The vignette says so where it shows them, rather than quietly differing.

## the wrapper is more accurate than TDA's own file path

6.15.2.1 Box 2 prints Deviance 2142.4880; the wrapper reports
2142.48772347. Neither is a rounding of the other, and for once the
difference is in our favour.

Reading a data file, TDA stores each variable in the size its nvar
declaration implies -- 4-byte floats by default, which is what the
manual's runs used ("Idx Variable T S" shows S = 4). The wrapper hands
R's doubles over directly, S = 8. The fitted coefficients agree either
way (-143.0269, 3.8990); the deviance does not.

stats::glm on the same data gives 2142.48772347 -- exactly the wrapper's
value, to twelve digits. So the wrapper is right and the command-file
path loses precision in the seventh digit. Asserted in tests against
stats::glm rather than against the manual.

This is why "the manual says otherwise" cannot be the end of a check in
either direction: three boxes today were the manual being out of date,
one was a regression of ours, and this one is the program itself being
less precise than the wrapper.

## gcon's edge layout was unnamed, and 6.2.2's last digits explained

`tda_g_components()` returned V1 to V7 for its edge layout. Named after
the manual's own headers (7.2.3.1, Box 1): component, n, i, j, node_i,
node_j, value. The vignette showed only the node layout; with both, the
box reads 65/65.

6.2.2 Box 12's later rows differ in the last digit or two -- 103.47752
against 103.4777. Not a fault either way: gdf4.cf writes its generated
data to a file at [10.6] and reads it back, so its input is rounded to
six decimals, while the wrapper passes R's doubles straight over. The
same cause as the deviance in 6.15.2.1, and in both cases the wrapper's
value is the more accurate one.

## the glm family boxes printed only their coefficients

Sections 6.15.2.2 to 6.15.2.5 each fit two or three models and showed
the coefficient table alone. The manual's boxes report the degrees of
freedom, the deviance and Pearson statistics, the scale estimates and
the rank as well -- all of them already on the object as $stats. Six
fits now print them; 6.15.2.4 Box 3 reads 20/20, the others gained
between one and four.

What is left in those boxes cannot be matched: the echoed settings (max
iterations, tolerance) and the final scaled parameter change, which is
machine-epsilon noise -- 6.88141e-16 here against the manual's
2.75435e-16, both runs stopping at the tolerance from slightly different
directions. The residual classifier says so now rather than leaving them
in the queue.

## qreg's category table dropped the weighted total

TDA prints one more value on the N line, past the categories: the
weighted number of observations, which its own header calls the
"(Weighted) Observations" column. The parser took the first n values to
match the categories and discarded it.

Kept as the "observations" attribute of $categories -- 581 for the
multinomial logit of qr3.dat, which is what the manual's 6.12.4 Box 2
prints. One value per wave where a run has several.

## gsort's relabelled-edge layout was unnamed; the residual is now explained

`tda_g_toposort()` returned unnamed columns for both its layouts. Named
after the manual's own headers (7.2.2.1, Box 2): i/node_i/label_i for
the node list, and label_i/label_j/node_i/node_j/value for the graph
relabelled. The vignette showed only the first; with both, 57/57.

With that and the category total, the manual comparison has no
unexplained distinctive misses left. What remains in the residual is
named by tools/manual_residual.py: echoed settings, citations, figures
quoted in prose, cross-references, axis labels, memory figures, page
footers, convergence noise, command syntax -- and the boxes settled
against Rohwer's own binary, which the classifier now lists with the
reason rather than leaving them to look outstanding.

## two R CMD check warnings, and examples that never ran

`print..prate_args`: a roxygen block claiming `@exportS3Method
base::print` sat above two plain comment lines and then `.prate_args`,
an internal helper. Roxygen attached it there and registered a print
method for a class named ".prate_args" -- both warnings from one
orphaned block. Removed, with a note so it does not come back.

The undocumented arguments were worse than they looked: `pairs = FALSE`
was still in two signatures, left over from renaming it to `type` --
present, used nowhere, documented nowhere. Removed; `type` documented on
the family page, on tda_seqpd and on tda_seq_info. checkDocFiles,
checkS3methods and codoc all come back empty.

Examples: 217 run, none reports a failure. Five pages had examples that
were never executed. Three are legitimate -- an external data file, a
URL, and lines already runnable above the illustration. Two were not:
tda_arcv and tda_arcvc had only \dontrun{}. They need an archive
description loaded by arcd first, so tda.zad and tda.zoo now ship in
inst/extdata and both examples run for real (four records written).
tda_plot_ps's example drew axes with no ranges, so there was nothing to
render; it now plots five points on a quadratic.

tests/testit/test-examples.R enforces this: no example may report "did
not converge", "cannot handle", "estimated nothing" or an R error, and
none may print nothing unless it draws or is guarded by
requireNamespace(). Both exception lists are named in the test rather
than tolerated silently.

## tda_arcv could never have worked, and its example hid that

`arcv` reads an archive that `arcd` has loaded, and TDA keeps no state
between runs. `tda_arcv()` emitted a bare `arcv` command, so every call
answered "no data archive defined". It now takes `archive=` and issues
the `arcd` itself.

The example was the tell, and I made it worse before I made it better:
asked to give the page a runnable example, I wrote one that drove TDA
with `tda_run(c("arcd = ...", "arcv() = ..."))` -- a working example
that never calls the function it documents, and would have passed the
new test while demonstrating nothing.

Checked the whole package for that shape afterwards: `tda_arcd`'s
example did the same thing for no reason at all -- that function works
standalone. Both now call what they document. `tda_options` is a
documentation-only topic with no function behind the alias, and is
exempt.

Three assertions in test-examples.R: no example may report a failure, be
silent unless it draws or is guarded by requireNamespace(), or fail to
call the function it documents. The failure strings are read out of the
sources rather than guessed -- an earlier version looked for "cannot
handle", which appears nowhere, so it would have passed whatever
happened. Validated against three real failure modes and one healthy
fit.

### the same roxygen trap, twice in one session

Factoring the arcd command into `.arcd_line()` I put the helper between
`tda_arcd`'s roxygen block and `tda_arcd` itself -- so roxygen attached
the documentation to the helper and wrote man/dot-arcd_line.Rd. Exactly
the mistake that produced the `print..prate_args` warnings earlier in
this same session.

The new example test caught it immediately, which is the first time a
test written today has caught a fault made today. The helper now sits
above the block, with a comment saying why it must stay there.

## the contingency measures are formatted in R, not passed through

`print.tda_freq2` echoed TDA's own text for the contingency block. TDA
pads those label/value lines to a fixed column, which is correct in
characters but leaves the column at the mercy of the reader's font: a
face that draws "1" narrower than "0" and centres it in its cell puts
the ink of 1.0000 about 8px right of 0.7071 in a 9.6px cell, and the
column reads as ragged though the characters line up.

`$measures` holds the same numbers, so the method formats them itself:
one width for every value, computed from the values, with the labels
padded to the longest. Every other line TDA printed in that block is
kept -- the cell counts, the degrees of freedom -- selected by label
rather than by shape, after a first attempt at "a label followed by a
number" swallowed "Degrees of freedom: 3".

### and the vignette shows them as a table

Formatting in R fixes the characters but not the rendering: the reader's
font still decides where the ink sits inside each cell, and on a face
that draws "1" narrow and centred a column of 1.0000 and 0.7071 reads as
ragged however the text is padded. Changing the font stack (ui-monospace
first, tabular-nums) did not settle it either.

So the vignette lays $measures out as a table, values right-aligned in
their own cells. Alignment is then layout, not typography, and no font
can undo it.

## the plot folding had eaten two results

5.5.2's running-median box computed the three smoothers with
tda_derive() and printed them, then drew them -- all in one chunk, which
the plot-folding pass then folded whole, hiding the smoothers
themselves. Split: the smoothing and its table in a box of their own,
the drawing in a folded one after it. 6.10.2's npreg box had the same
shape and is split the same way.

Both were found by the check for a folded chunk holding a result rather
than by reading: the earlier version of that check looked only for model
fits, so tda_derive() and tda_npreg() slipped through it. It now covers
the derive/smooth/matrix families too.

Gallant's 30 observations (6.16.1) were typed into the vignette as an
inline read.csv(text = ...). They ship as inst/extdata/exam/gallant.dat
now and the box reads the file, which is why its fingerprint drops from
208 numbers to 55 -- the data is no longer echoed as source.

## seqpd's box showed one sequence structure where the manual shows two

The manual's 3.4.3 runs seq5.cf, which reads seq.d3 and declares TWO
structures: Y0..Y7, one column per time point, and three state/time
pairs read as type 2. Every table it prints carries a column per
structure, and the common time axis runs 0..8 rather than 0..7. The
vignette was running a single structure over seq.d1, so its tables had
one state column and one time point fewer -- the same shape, different
data.

`.seq_desc()` now takes `second =` (and `second_type =`), emitting a
second `seqdef(sn = 2, m = ...)` line for those columns; `tda_seqpd()`
passes them through. seq.d3 ships in inst/extdata/exam.

Checked against the manual's own Box 3: layout 3 gives 1/0/7/1, 1/1/3/1,
... 1/8/-1/-1, and layout 4 gives 1/1/7/3/0/1 through 1/5/7/-1/6/8, both
exactly as printed.

## R substitutes for TDA's own commands

Asked whether the convex hull figures could be reproduced exactly, and
they can: plot17.cf is 100 cases of X = rd, Y = rd at the default seed,
and tda_rng() gives TDA's own draws digit for digit (0.0291038/0.94947,
0.094304/0.700125, 0.891351/0.47206). Our PostScript for that figure has
27 drawing operations, identical to the ones TDA's own run emits.

Reproducing it turned up two places where R was standing in for TDA:

  tda_pl_hull()  computed the hull with grDevices::chull() and drew it
                 as a polyline. plotch IS TDA's convex hull command --
                 "plots a convex hull around the points given by the
                 variables on the right-hand side" -- and the real
                 wrapper for it was sitting under the name
                 tda_pl_curve(). The consequence was not cosmetic:
                 plotch's smoothing (ns=) and the margin it adds to a
                 smoothed hull (ic=) had no way in, so the manual's own
                 Figure 2 (plot18.cf) could not be drawn at all -- ic was
                 rejected as an unknown option. tda_pl_hull() now calls
                 plotch, with smooth= and expand=.

  tda_rd(),      wrapped stats::runif and stats::rnorm while documenting
  tda_rdn()      themselves as nvar()'s own rd() and rdn(). They gave
                 R's numbers, so an example built with them could not
                 reproduce a figure from the book. Both now draw from
                 tda_rng().

A test asserted the second one -- "tda_rd == runif under the same seed"
-- which is how it survived. It now asserts TDA's own first three
uniforms. A second test expected a two-point group to be skipped, which
was the R hull's behaviour; TDA draws a degenerate outline for it, and
the test says so now.

Swept the rest: of 325 commands in tda.hlp, 24 are not mentioned
anywhere in R/ (arcc, cmean, seqdel, seqrd, cblen, mnvar, the macro and
shell commands, spss/stata readers). stats::lowess in tda_pl_regression
stays: xreg's X11 front end is not in this source tree, and that curve
is verified path-for-path against TDA's own scplot lowess.

## npreg's lowess is not scplot's, and the plot wants scplot's

Asked for TDA's lowess rather than R's, which meant finding out why the
plot wrapper had reached for stats::lowess. The reason was a wrong note
in tda_npreg(): opt=4 was recorded as "a genuine TDA syntax error" and
excluded. It is not -- the error was the missing x= parameter, which
every npreg method needs. With x= supplied TDA answers "Method: Lowess:
sig=0.5 ns=2 d=1", and opt=5 gives "Method: Midmeans". Both are
available now.

But npreg's lowess is NOT the curve scplot draws, and cannot be: npreg
resets its band width to 1 when it is below epsilon, BEFORE the lowess
branch, so the delta shortcut is never 0, while scplot calls lowess with
delta hard-coded to 0. Measured on forty points, the two differ by up to
0.47. Setting d=0 through npreg does nothing -- TDA still reports d=1.

So tda_pl_regression's lowess now delegates to tda_pl_scatter, which
draws scplot's own curve. A test that both wrappers trace an identical
PostScript path passes. No stats::lowess call remains in the package.

## 6.5.5 Box 2's second table

The box shows ed1.dat read twice: as single episodes, and as
multi-episode. The second reading needs id= and sn= in the edef block --
the manual's own Box 1 shows them commented out for the first. Without
them a case's later spells count as missing, which is the whole
difference: 3 3 3 3 3 3 4 against 0 0 0 0 0 0 1.

tda_state_dist() now takes id= and spell=, passed to .tda_design() so
the Id and Sn columns exist for edef to name. The box reads 90/90.

Worth noting: the first attempt at this edited tda_episodes() instead --
the two functions share the same two lines, and a string replace found
the wrong one. Reverted. That is the second such slip in this session.

## reading the source instead of probing the binary

The npreg/scplot question was settled by reading t_areg.c and t_plot.c,
after four rounds of trying options against the binary. The answer was
three lines of C. Worth stating as a rule: the sources and a build of
Rohwer's own are both here, and a question about what a command does is
answered by reading it, not by guessing at inputs. The opt=4 "syntax
error" note that caused this was itself a black-box conclusion recorded
as fact.

Applied the same way to plotch, before writing anything: ns= is Akima
smoothing over that many intervals, and TDA smooths only when ns >= 2
AND the hull has more than two vertices. ic= is a margin in millimetres
-- each smoothed point is pushed from the hull's centroid in proportion
to its distance along each axis, via PMIC * UXLen / PXLen. Both are
documented in those words now.

### and the figures that used R's generator

4.4.9's two hull figures drew from set.seed(1) + runif(), with a note
saying TDA's rd() and R's runif() "are different generators, so the
specific points differ". They need not: tda_rng() is TDA's own stream,
and the figures are now the manual's own points. 4.4.x's plot8 and
plot11 likewise used runif() and pnorm() where the command files use
rd(-2,2) and nd(); both now use TDA's, through tda_rd() and
tda_derive(ND = "nd(RDS)"). No set.seed remains in the manual vignette.

## working the shape check: the first four

tools/check_manual_shape.py listed 114 boxes printing fewer numbers than
the manual's. Four cleared:

  6.5.1   ltb.1 holds two blocks -- the counts table and the estimates
          with their standard errors. Only the first was printed;
          $survivor was on the object already.

  7.2.4.1 the box ran two of gep's six output options, and said so in
          the prose. All six now; three boxes cleared at once.

  7.2.6.1 gcyc's option 4 returns 49 records as a ragged list, and
          print.tda_ragged shows ten by default. print(n = Inf).

  7.6.1.2 gio had no output box at all -- the section drew the graph and
          stopped. Added, on gd2b.dat as gd11.cf runs it, and its seven
          columns named after the manual's own header (i, node_i,
          node_j, a_ij, y_ij, c, iterations). gd2b.dat now ships.

Note what these have in common: every one is a box that scored full
marks on check_vignette_manual.py while missing half its content or
more. The membership test cannot see them; the count test finds them
immediately.

## shape check, second batch

  7.2.6.2  gdcyc ran option 1 only; the manual's box shows both.
  7.2.4.3  gtcl ran opt=2 twice, on the undirected and directed graph;
           the manual's box is unvalued AND valued, so both options of
           each.
  6.17.2.1 the rate file opens with the constellation it was computed
           at -- one row per parameter, its coefficient and the
           covariate value used -- which was read past as comment. It
           is the only record of which covariate values a rate belongs
           to, so it is kept as the "constellation" attribute of
           $rates. The box also computed its rates at W = 0 where
           rt1.cf uses COHO3 = 1, W = 1, so every rate in it was the
           wrong constellation's.

That last one is worth noting: the box was not missing content, it was
showing the right table computed from the wrong settings. The membership
check scored it 145/145 because the numbers it printed were real numbers
from a real run -- just not the run the manual describes.

## shape check, third batch

  7.3.1.2  gqap asks for df=q.aux in the command it builds and then
           discarded the file. It is the permuted multigraph, the
           manual's own second output file: each pair, its nodes under
           the permutation, and the flow and distance between them. Now
           $permuted.

  7.2.1.1  gni had no box at all -- tda_g_nodes() exists and the
           vignette never called it. Added on gd1.dat as gd4.cf runs it.

           Its node 8 reads in-degree 0 and out-degree 0 for the first
           graph where the manual prints 1 and 1. TDA gives 0 and 0, and
           so does a build of Rohwer's own sources, so the manual is out
           of date here -- the fifth such box, and the vignette says so
           where it shows it.

That is three findings of three different kinds from the same list: a
second output file discarded (gqap), a command never exercised (gni),
and a manual box that no longer matches the program (gni's node 8).

## shape check, fourth batch

  6.7.2.3  seqm's s= option had no box. print = "sequential" gives the
           distance after each time point (seqm1b.cf) and print = "lcs"
           the longest common subsequence with its positions
           (seqm1c.cf); both were reachable and neither was shown.

           The box also prints TDA's own description of the file it
           writes -- the nvar block naming every column, its width and
           its meaning. dtda= was never requested, so that file was not
           written at all. Asked for now and kept as the "description"
           attribute of the pairs table.

           Three boxes cleared at once.

Worth recording: the box was first placed after the next ### heading, so
it belonged to 6.7.2.4 and the checker went on reporting 6.7.2.3 as
short while the content was right there. Section membership is by the
heading above, and a box written for one section must sit before the
next one starts.

## shape check, fifth batch, and sharpening the tool itself

  6.9.1.3  the box printed the group table and the coefficients but not
           the fit statistics the manual reports beside them -- df, the
           residual sums and norms, R^2 and its adjustment, the rank.
           All were on $stats already.

  6.12.3   qreg writes a description of its predictions file (the nvar
           block naming every column) and the reader used it for the
           column names and dropped it. Kept as the "description"
           attribute. Our widths read [24.16] where the manual's read
           [10.4] -- the wrapper asks for full precision so the
           predictions are not rounded by the print format -- so the box
           is settled rather than short.

Two changes to the tool, because the queue was reporting work that was
not there:

  settled boxes are listed apart, with the reason, instead of sitting at
  the top of the list (9 of them, each checked against the program and
  five against a build of Rohwer's own sources)

  a value the manual prints to more digits than we do is counted as a
  precision difference, not a miss. 6.9.1.3 was 25 short entirely on
  this -- the manual writes 1.045236e+00 where a summary prints 1.0452.

99 boxes left, and they are boxes with genuinely missing content.

## shape check, sixth batch

  5.2.6.6  the manual's box has five change operators -- change, cntch,
           ccntch, lagch, and lag(T,-lagch(S)) -- and the box showed two.
  5.2.6.7  likewise five aggregate operators: gcnt, grec, gsn, gfirst,
           glast. The box showed gcnt and grec, plus gmean and gstd,
           which the manual demonstrates elsewhere.
  8.5.1    the manual's box prints the two input data files beside the
           results, and only the results were shown.

All three were reachable with what the package already had -- no wrapper
changed. They were boxes written to show a couple of examples of a
family where the manual shows the family.

## shape check, seventh batch

  7.2.1.2  gdln has forward and backward links, and tda_g_links() takes
           backward = TRUE. The box showed forward only.
  6.7.2.3  the dp_matrix entry carries seqA and seqB beside D -- the two
           sequences the matrix aligns, which the manual prints above
           it. Only D was shown.

           What is left there is the manual printing the matrix's row
           and column indices 0..8 as headers where ours labels by the
           sequence values. Same matrix; recorded as settled.

Checked before assuming: our D and the manual's agree cell for cell
(0.00 1.00 2.00 3.00 / 1.00 0.30 1.30 2.30 / 2.00 1.30 0.70 1.70) once
the substitution cost matrix seqm3.cf uses is supplied. A first look
without it showed plain integers and would have read as a real
difference.

## shape check: teaching the tool what is not output

The queue went from 95 boxes to 80 without touching the vignette, by
excluding two classes the residual classifier already knew about:

  prose between boxes -- a box's caption runs to the next caption, so
  citations, page numbers and sentences like "convergence is achieved
  with 12 function calls" land inside it. A line of eight or more words
  is prose, not a table row.

  the settings TDA echoes before it estimates -- the algorithm, the
  iteration cap, the tolerances, the covariance type. A vignette box
  prints the fit, not the echo.

6.14.1 Box 3 was 22 short and is almost entirely those two.

Then, on the vignette:

  7.2.9.2  gcset has two layouts, the sets themselves and one row per
           node; the box called it with neither option named.
  7.6.1.1  gbcf likewise -- gd14a.cf asks for opt = 2 and the box showed
           only the default.

Both wrappers already took opt. Checked first that our row counts match
TDA's own run of gd14.cf (16 and 3) before assuming the box was short.

## shape check: the matrix boxes

5.1.4.5, 5.1.4.7, 5.1.4.8, 5.1.4.9 and 5.1.4.13 all failed the same way:
the manual's box prints its input matrix (mpr(A)) and then the result,
and the vignette built the matrix inline inside the call so only the
result was shown. Named and printed first in each. 5.1.4.13 states its
problem as one matrix T -- objective, inequality, equality, one per row
-- which tda_mlp1() takes as three arguments, so T is written out
alongside for the reader to compare.

Four of the five cleared; 5.1.4.5 and 5.1.4.13 are down to single
figures.

## shape check: gtcl's four graph types, and two more tool rules

7.2.4.3's box runs gtcl over all four of gdd's graph types -- gt=1
undirected unvalued, 2 undirected valued, 3 directed unvalued, 4
directed valued -- and the vignette ran opt=1 and opt=2 on two graphs,
which is a different thing entirely. tda_graph()'s `directed` argument
also fights `options = list(gt = ...)`: passing gt alone gives all four
correctly, passing both does not. Worth knowing.

7.6.1.2 Box 6 puts two gio runs side by side, IT 1 and IT 4. gd11.cf
produces the right-hand one exactly; the left-hand one is from a run
that is not among the shipped command files, and mxit = 1 gives 0.505
where it prints 0.5000. Settled rather than chased.

6.2.2 Box 10's rows match exactly -- its "extra row" was a figure's axis
labels, since a box runs to the next caption and a figure in between
lands inside it. The shape tool now skips evenly spaced integer runs, as
the residual classifier already did.

## seq_info takes a type per structure

3.4.2 Box 10 shows seq3.cf declaring two sequence structures and
printing the table after each: structure 1 type 1 over Y0..Y7, time axis
0-7, states 1 3 7; structure 2 type 2 over three state/time pairs, time
axis 0-8, states 1 2 3. `type` was a scalar applied to every structure,
so a run mixing the two could not be described at all. It is recycled
per structure now, and the box prints both.

## the shape check was over-reporting by a third

Its matcher took the FIRST pool value within tolerance. A token the
manual prints without decimals carries a tolerance of 0.5, so a "1"
could consume our 0.5 and leave the real 1 unmatched, cascading through
the rest of the box. 7.2.4.3 read 14 short when a direct count of its
values showed one -- a single zero.

Matching the closest value instead took the queue from 74 boxes to 52.
Found by not trusting the number: the four gtcl tables were checked row
by row against the manual first, and they agreed, which is what made the
14 worth doubting.

6.14.1 cleared on the vignette side: the box printed the fit but not the
convergence detail the manual reports around it.

## int reports more than the value

5.4.1's box prints the rule TDA used, the relative error it worked to,
the number of function calls and whether it succeeded, beside the
approximation. tda_integrate() returned the value alone. All four are on
the result now: 43 calls, 1e-04, method 1 (QNG) -- the manual's own
figures.

(The first version of the parser ran regexpr on the substituted string
and regmatches on the original, so every value came back NA.)

## the shape check, matching in two passes

A coarse match could consume a value an exact match needed: the manual
writes 1.0308e-01 where a fit prints 0.1031, and matching that loosely
first ate a value another token wanted exactly. Exact matches are made
for the whole box first, and only what is left over is retried at four
significant figures.

That plus the closest-match fix took the queue from 74 boxes to 51, and
the largest shortfall from 19 tokens to 7. 6.17.3.1's two boxes cleared
on it alone -- they print everything the manual does, at four decimals
where the manual uses five significant digits.

## the book vignette, first findings from the shape check

  ehe3  TDA prints two tables before the estimates -- the episodes it
        built, and how they fall across the periods -- and the box
        showed neither. Both were on the fit already, as $episodes and
        $periods, and both match the .ref exactly.

  ehd6  likewise $episodes.

  ehd7  the box used Surv(TFP, DES), each spell's own duration, where
        ehd7.cf measures time from the start of the FIRST job: TFC is
        the first spell's starting time carried forward, and the episode
        runs from TStart - TFC to TFin - TFC + 1. The reported episode
        table differed because of it -- TS Min 0 where TDA reports 13,
        TF Max 407 where TDA reports 450.

        The coefficients were NOT affected and never were: for an
        exponential model the exposure is sum(tf - ts), which is the
        same under both clocks. Checked rather than assumed, after
        first writing this up as though the estimates had been wrong.

## the episodes table, across the book vignette

Twelve of the short book boxes have an episode table in their .ref and
eleven did not print it. Fixed for ehi1, ehi2 and ehf4, which needed
their fits naming first -- the boxes called tda_rate()/tda_coxph()
without keeping the result, so there was nothing to ask for $episodes
from.

That table is how TDA reports what it actually modelled: how many
episodes per transition, their weighted count, mean duration, and the
earliest start and latest end. It is the first thing to check when a
model does not look right, and the boxes were printing coefficients
without it.

## the book vignette, episode tables across the model families

ehi7, ehi9, and the two chunks that fit five distribution families at
once now print the episode table their .ref reports. The five-family
chunks fit the same episodes five times over, so one table serves.

ehd7's remaining shortfall was TDA's Sum rows -- a total per spell, and
one for the whole table. Those are sums of the rows themselves:
tapply(episodes, sn, sum) gives 201, 162, 107, 62 and 532, exactly as
TDA prints them. A box showing the table has not lost them, so the check
skips them now.

78 boxes to 70.

## the book vignette, second sweep

ehf6, ehg4, ehe1 given their episode tables; ehf1 and ehf5 now say how
many episodes the split produced -- TDA reports "761 records written"
and "1021 records written", and tda_split() gives exactly those counts,
which the boxes had but did not state.

70 boxes to 62.

## tda_split discarded epdat's episode table

epdat reports the episodes it built before writing them -- one row per
transition with its count, weighted count, mean duration and the
earliest start and latest end, the same table tda_rate() carries as
$episodes. tda_split() returned a bare data frame and read past it. It
is the "episodes" attribute now, and the two split boxes print it.

Also given their episode tables: ehd1, ehc3 (which cleared nine boxes at
once, the life-table sections covering several book boxes each), and
ehi3, whose chunk fits three Cox models to the same split episodes.

62 boxes to 49.

## ehf4, ehf3: a route difference, not a fault

ehf4.cf splits inside edef(split = MarrDate) and so reports the episodes
BEFORE the split -- 600, as 142 censored and 458 events. The box splits
first with tda_split() and reports the split episodes. Both describe
their own route correctly, and the coefficients match the .ref exactly,
all nine rows. Listed as explained rather than short.

Worth recording how nearly that went the other way: a throwaway
diagnostic said 0.5423, 4.7916 and -0.0264 were missing from our box,
which would have meant the estimates were wrong. They are all there. The
script matched the FIRST value within tolerance instead of the closest
-- the same bug that was fixed in check_manual_shape.py, reproduced in
an ad-hoc check. The coefficients were compared row by row against the
.ref before drawing any conclusion.

ehi10 given its episode table. 49 boxes to 43.

## Figure 10.2.1 was drawn wrong, and prate cannot draw it

The first attempt put two flat lines on the page. The point of the
figure is that the rate DECLINES once a gamma mixture is allowed, so it
showed the opposite of what it exists to show.

prate does not apply the mixture. TDA's own run of ehd2.cf with mix = 1
gives the same constant 0.0186 at every time point, and t_prate.c
contains no reference to PMMIX at all -- checked in the source rather
than inferred from the output. So the curve cannot come from prate.

The individual rate is constant; what declines is the rate observed in
the population, because those with the highest rates leave first. That
is h/(1 + sigma^2 H(t)), with sigma^2 the mixture variance -- the fit's
"D Constant" on the exponential scale, -0.8019 here -- and H the
cumulative rate. Computed from the fitted parameters and drawn, with
both curves labelled and the axes titled.

Then drawn WRONG a second time: the dashed line was the mixture model's
own individual rate, 0.0186. The book compares two FITS -- "without
(dotted) and with (solid) a gamma mixture distribution, model
specification according to Box 4.1.5" -- so the dashed line is the rate
from the model fitted without the mixture, 0.0112. The two curves now
cross at about month 80, as the book's does. Axes are the book's, 0 to
0.02. The prose says plainly that prate will not produce it and why.

Also fixed on the way: lty = "dotted" is not a line type TDA takes, so
the two curves were both solid and indistinguishable; tda_pl_text()
takes the label before the position, and the first call passed
coordinates first; and the x-axis title overlapped the tick numbers.

## Figure 10.2.1 removed

Written three times and wrong three times; removed rather than attempted
a fourth. The record above describes the first two errors. The third
fixed the baseline and still did not match the book, so the vignette
carries the two mixture fits (Boxes 10.2.1-10.2.2) and no figure. The
book's Figures 10.2.2 and 10.2.3 were never drawn either.

Nothing else in the vignette depended on it: the prate= arguments added
to the mixture fits for its sake were removed with it, and the QA table
row was withdrawn.

## ehhnew caption tidy

One caption listed its figures twice -- "Figures 7.2.2, 7.3.2, 7.4.2,
7.5.2, 7.6.2 [book pp. 21-206, Fig. 7.2.2, Fig. 7.3.2, Fig. 7.4.2, Fig.
7.5.2, Fig. 7.6.2]" -- and two others repeated a single figure number in
the same way. The title now says what the chunk draws and the bracket
carries the citation.

That "pp. 21-206" was a real error, not a typo in the caption: the QA
table has Fig. 7.6.2 and Box 7.6.3 on "page 21", a digit the scan lost.
Both sit between 7.6.2 on p.210 and 7.6.4 on p.212, so p.211. The table
and both captions are corrected.

Note for next time: shortening "ehh1, ehh2, ehh3, ehh4" to "ehh1 to
ehh4" in a caption broke check_vignette_plots.py, which reads the
command files out of the caption text. Restored.

## mdsn: Kruskal's step adaptation restored; mdsm and mdsn out of "experimental"

`t_mds.c`. The `mdsn` command as shipped carried Kruskal's (1964) step
length adaptation -- the angle, relaxation and good-luck factors --
inside a comment, and stepped with a fixed `slen = 0.2` and `mag = 1.0`
instead; its iteration protocol went to stdout through `tda_out` debug
prints ("PTR:", "HIER MDSXB:", "MDSRN=..."), and `df=` and `pcf=`
printed "N records written" while the writing code was commented out. On
`eurodist` it stopped at stress 0.365 where `MASS::isoMDS` reaches
0.075.

The commented block is now the live update (`slen *= af * rf * gf`, step
`slen * g / (|g| / |x|)`), the debug prints are gone, `df=` writes the
fitted distance matrix and `pcf=` the Shepard triples (dissimilarity,
fitted distance, monotone projection) in ascending order. Stress-1 as
TDA defines it was already right (recomputed from the returned
configuration to three digits). Result: eurodist 0.058 (isoMDS 0.075),
swiss 0.0936 (0.0941), four random problems 0.003 -> 0.0001. Control:
`examples/tests/mdsn.cf`; `examples/coverage/mdsops.ref` regenerated.

`mdsn1` carried the same commented block and gets the same restoration,
plus removal of its live debug prints (`b=======`, `s=... st=...`,
`mag=`). An earlier version of this entry said the adaptation made
`mdsn1` worse; that came from edits that had landed in `mdsn` instead (a
`str.index` on the function name matched its prototype at the top of the
file), so the two builds compared were the same build. Measured properly
(a build with `slen *= af * rf * gf` replaced by `slen = SLen` in
`mdsn1` alone): fixed step 0.06-0.19 on the six embedded problems,
adapted 0-1e-4 and 0.0078/0.058 where isoMDS gets 0.0078/0.075; on the
30 random noisy problems adapted <= fixed in 30/30 and within 0.002 of
isoMDS in 28/30 against 1/30 for the fixed step. Its "(experimental)" is
dropped as well. `mdsn1`'s gradient (`mdsn1_fn`) checks out: with a tiny
step its stress falls monotonically.

`mdsm` was complete as shipped: no debug output, `df=` written, restarts
find distinct minima, raw stress on eurodist 0.072 against cmdscale's
0.090. It fits two dimensions only (`NParm = 2 (n - 1)`), which the R
wrapper now refuses for any other `ndim` instead of passing it through.
The "(experimental)" in the `mdsm`, `mdsn` and `mdsn1` banners is
dropped; `mdsx` keeps its.

(An earlier version of this entry reported `rfit`/`rfit1` hanging on
`rel=3,4,5`. They do not: TDA returns "exceeded maximal number of
solutions" at once when more optimal relations exist than `max=` allows
-- 17 linear orders for 4 nodes at the default 10 -- and my test harness
mistook the resulting non-structured return for a timeout. `tda_rfit()`
now turns that message into an error naming `max_solutions`.)

## mdsn1: the gradient was accumulated across iterations

`t_mds.c`, `mdsn1_fn()`. The gradient array was filled with `+=` and the
caller never cleared it, so from the second iteration on the "gradient"
was the running sum of every gradient since the command began -- across
the random restarts too. With the fixed step this behaved like heavy
momentum and happened to reach reasonable minima; with Kruskal's step
adaptation on top it fell apart (random problem 5: stress 0.0078 ->
0.176, which the previous entry misread as the adaptation being
unsuitable for this command). Assignment instead of accumulation, and
the adaptation restored here as in mdsn: mdsn and mdsn1 now reach the
same stress on every test problem (eurodist 0.058007 both, six random
problems to 6 digits), and mdsn's gradient formula agrees with a
finite-difference gradient of Kruskal stress-1 to 4e-11 (checked in R).
The "(experimental)" in mdsn1's banner stays for now: its
parametrisation fixes the first point at the origin and is
two-dimensional only.

## rxls: the Workbook stream is read through the OLE2 sector table

`t_xls.c`. `rxls` located the BIFF BOF by scanning the raw .xls file and
then read on sequentially. An .xls file is an OLE2 compound document and
its Workbook stream need not be stored in consecutive sectors: readxl's
`datasets.xls` of 1.4.3 has a break in the chain at sector 65 -> 67, and
the read ran into the foreign sector 66 ("Error in length: -31744", the
first two bytes of it as a signed short). `xls_ole_extract()` now walks
the header, DIFAT, FAT and directory (and the mini FAT for a stream
below the 4096-byte cutoff), writes the stream in order to
`tda_rxls.tmp`, and the parser reads that from offset 0; the file is
removed at the end. A file that is not an OLE2 document is read as
before. The loop over BIFF substreams also stops cleanly at the end of
the extracted stream instead of reporting a read error there. Controls:
readxl's datasets.xls of both 1.4.3 (4 sheets, 1253 rows) and 1.4.5 (3
sheets, 1103 rows), clippy, type-me and deaths; test-phase3.R no longer
depends on which of the two datasets.xls is installed. Still unread,
before and after: geometry.xls (no cells found).

## idf opt=2: the self-consistent distribution function, finished

`t_imat.c`, `idf()`. The manual (8.5.2) defines the self-consistent
distribution by the functional equation Pr(S) = (1/n) sum_i Pr(o_i and
S) / Pr(o_i) and says its calculation "is not yet available with the idf
command". The C had the iteration: `opt=2` read, `mxit=` honoured, the
step computed into `AcU` -- and the update `AcZ = AcU`, the convergence
test and every line of output inside a comment, so `opt=2` printed
"Iterative calculation of self-consistent distribution." and nothing
else. The R wrapper had recorded this as "idf does not read opt=".

Finished the way `iddf()` finishes its discrete twin: the iterate
becomes the current distribution, the loop stops when the largest change
is under `tolf` (default 0.001, as iddf; `tolf=` overrides), then
"Convergence [not] reached after N iterations." and a table of the
partition points with the distribution function (export `idf.scdf`).

Validation (test-gap-wrappers.R): TDA's fixed point equals the manual's
equation iterated in R from the same mean-df start to 2e-9; applying the
equation once to TDA's own F moves it by 1e-10; F ends at 1 and is
nondecreasing to 4e-14. On the manual's income example (`id2.dat`,
`examples/tests/idf_sc.cf`) the mass sits on the innermost intervals
(2/11, 4/11, ...) as a Turnbull-type estimate does. `tda_idf()` now
takes `self_consistent=` and `control=` like `tda_iddf()` and returns
`$self_consistent` and `$converged`.
