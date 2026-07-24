## The seq-family variables=/id= interface (tda_seqgc/tda_seqlg/tda_seqsd/
## tda_seqsi/tda_seqen/tda_seqpm) against seq.d1, TDA's example data for
## these commands (examples/exam/seq9.cf runs seqlg, seqgc, seqsd and seqen
## on it together; the .ref values checked below were read from a real run
## of that binary, not guessed).
##
## seq.d1 is ID, Y0..Y7 -- ID 1:5, not a sequence value, and including it as
## an extra time-point column is the exact bug this interface exists to
## rule out (an extra "time point" always present and never -1 shifts ts/tf/
## slen and, for a case with a leading gap, turns it into an internal one).


seqd1 <- utils::read.table(system.file("extdata", "exam", "seq.d1", package = "tdaR"))
names(seqd1) <- c("ID", paste0("Y", seq_len(ncol(seqd1) - 1L)))
yvars <- paste0("Y", 1:8)

## -- variables= selects the sequence columns, ID left alone -------------

lg1 <- tda_seqlg(seqd1, variables = yvars)
lg2 <- tda_seqlg(seqd1[, -1])
# Each result now carries its run (a temp directory, a timestamped
# banner), which necessarily differs between two runs; what this checks
# is that the two ways of naming the sequence columns give the same
# table, so the run is dropped before comparing.
no_run <- function(x) { attr(x, "run") <- NULL; x }
ok("seqlg: variables= and pre-subsetting agree",
   identical(no_run(lg1), no_run(lg2)))

same("seqlg: ts/tf/slen/glen/gmin/gmax/ngap match the manual's Box 1",
     as.matrix(lg1[c("ts", "tf", "slen", "glen", "gmin", "gmax", "ngap")]),
     matrix(c(0, 7, 8, 0, 0, 0, 0,
              0, 7, 8, 2, 2, 2, 1,
              1, 7, 7, 1, 1, 1, 1,
              0, 5, 6, 0, 0, 0, 0,
              0, 7, 8, 0, 0, 0, 0), nrow = 5, byrow = TRUE,
            dimnames = list(NULL, c("ts", "tf", "slen", "glen", "gmin",
                                    "gmax", "ngap"))))

## the bug variables= exists to prevent: ID read as a ninth time point,
## always present and never -1, shifts tf (last valid position) by
## exactly one for every case, but slen = tf - ts + 1, so it moves by
## more wherever ts itself shifts too -- case 3 has a leading gap, which
## an always-present extra column at the front turns from "before the
## sequence starts" into an internal gap, so ts stays 0 (not 1) and ngap
## doubles as well.
lg_bug <- tda_seqlg(seqd1[, c("ID", yvars)], variables = c("ID", yvars))
ok("seqlg: ID included as a sequence column shifts tf by one",
   all(lg_bug$tf == lg1$tf + 1L))
same("seqlg: and turns case 3's leading gap into an internal one",
     c(lg_bug$ts[3L], lg_bug$slen[3L], lg_bug$ngap[3L]), c(0, 9, 2))

## -- id= carries a column through untouched, one value per case ---------

gc1 <- tda_seqgc(seqd1, variables = yvars, id = "ID")
ok("seqgc: id= attaches the requested column", "id" %in% names(gc1))
same("seqgc: with the case's real value, in order",
     gc1$id, seqd1$ID)

## -- seqsd/seqen: no per-case row, id= must error, not silently drop ----

ok("seqsd: id= errors instead of silently accepting and dropping it",
   inherits(try(tda_seqsd(seqd1, variables = yvars, id = "ID"),
                silent = TRUE), "try-error"))
ok("seqen has no id= parameter at all",
   !("id" %in% names(formals(tda_seqen))))

sd1 <- tda_seqsd(seqd1, variables = yvars)
same("seqsd: state distribution matches a real run of seq9.cf's seqsd",
     as.matrix(sd1[c("nst1", "nst3", "nst7", "valid", "nmiss", "total")]),
     matrix(c(1, 1, 2, 4, 1, 5,
              1, 4, 0, 5, 0, 5,
              0, 4, 1, 5, 0, 5,
              3, 1, 0, 4, 1, 5,
              2, 3, 0, 5, 0, 5,
              1, 3, 0, 4, 1, 5,
              1, 0, 2, 3, 2, 5,
              1, 1, 2, 4, 1, 5), nrow = 8, byrow = TRUE,
            dimnames = list(NULL, c("nst1", "nst3", "nst7", "valid",
                                    "nmiss", "total"))))

en1 <- tda_seqen(seqd1, variables = yvars)
same("seqen: entropy matches a real run of seq9.cf's seqen",
     en1$ent, c(1.0397, 0.5004, 0.5004, 0.5623, 0.6730, 0.5623, 0.6365,
               1.0397), 5e-5)

## -- seqsi: per-case-time row, id= carries through the same way ---------

si1 <- tda_seqsi(seqd1, tp = "0(1)7", variables = yvars, id = "ID")
ok("seqsi: id= attaches the requested column", "id" %in% names(si1))
same("seqsi: one row per case", nrow(si1), 5L)

## -- seqpm: id= alongside pattern match counts ---------------------------

pm1 <- tda_seqpm(seqd1, list(c(3, 3)), variables = yvars, id = "ID")
ok("seqpm: id= attaches the requested column, not clashing with `case`",
   all(c("case", "ID") %in% names(pm1$table)))
## seq_search() drops any case with a zero-length sequence or an internal
## gap entirely (cases 2 and 3 here have gaps, per Box 1) -- id= should
## still line up with whichever cases TDA actually kept.
same("seqpm: id= lines up with the cases seq_search() actually kept",
     pm1$table$ID, c(1, 4, 5))

## -- error paths ----------------------------------------------------------

ok("variables/id overlap errors",
   inherits(try(tda_seqgc(seqd1, variables = c("ID", yvars), id = "ID"),
                silent = TRUE), "try-error"))
ok("unknown variables name errors",
   inherits(try(tda_seqgc(seqd1, variables = "NOPE"), silent = TRUE),
            "try-error"))
ok("unknown id name errors",
   inherits(try(tda_seqgc(seqd1, variables = yvars, id = "NOPE"),
                silent = TRUE), "try-error"))

## -- tda_seqm print="sequential" (s=1): replaces DIST, doesn't add to it --
## seqm.d1/seqm1.cf, TDA's example for seqm: 4 cases, real distances
## from Box 4 of the manual are 2,2,2,4,2,2 (checked against a real run of
## seqm1.cf directly, not the manual alone).

seqmd1 <- utils::read.table(system.file("extdata", "exam", "seqm.d1", package = "tdaR"))
names(seqmd1) <- paste0("Y", 1:5)

m0 <- tda_seqm(seqmd1)
same("seqm: default distance matches Box 4 of the manual",
     as.vector(m0), c(2, 2, 4, 2, 2, 2))

m1 <- tda_seqm(seqmd1, print = "sequential")
ok("seqm: print='sequential' still class dist",
   inherits(m1, "dist"))
same("seqm: print='sequential's distance matrix is unchanged from default",
     as.vector(m1), as.vector(m0))
p1 <- attr(m1, "pairs")
ok("seqm: print='sequential' names one dist_t* column per time point",
   all(paste0("dist_t", 1:5) %in% names(p1)))
same("seqm: distance is the *last* dist_t column, not the first",
     p1$distance, p1$dist_t5)
ok("seqm: the first running-cost column is not the final distance",
   !identical(p1$dist_t1, p1$distance))

byt <- attr(m1, "dist_by_time")
ok("seqm: dist_by_time is an n x n x n_time array",
   identical(dim(byt), c(4L, 4L, 5L)))
same("seqm: dist_by_time's names are t1..t5, in order",
     dimnames(byt)[[3L]], paste0("t", 1:5))
same("seqm: dist_by_time's last slice equals the main result",
     byt[, , 5L], as.matrix(m0))
ok("seqm: an earlier slice differs (less of the sequence considered)",
   !identical(byt[, , 1L], as.matrix(m0)))
for (t in 1:5)
    same(paste("seqm: dist_by_time[,,", t, "] matches pairs' dist_t", t),
         byt[, , t][cbind(p1$i, p1$j)], p1[[paste0("dist_t", t)]])

ok("seqm: print='sequential' + compare_with errors instead of mislabeling",
   inherits(try(tda_seqm(seqmd1, print = "sequential", compare_with = 1),
                silent = TRUE), "try-error"))

m2 <- tda_seqm(seqmd1, print = "lcs")
p2 <- attr(m2, "pairs")
ok("seqm: print='lcs' names its columns instead of leaving them V6, V7...",
   all(c("lcs_length", paste0("lcs_t", 1:5)) %in% names(p2)))
same("seqm: print='lcs's distance is unaffected -- adds, doesn't replace",
     p2$distance, p1$distance)

## -- tda_seqm subcost= as a full matrix (TDA's scost=NAME) ----------
## seqm.d3/seqm3.cf, TDA's example for a custom substitution-cost
## matrix; seqm3.d (shipped alongside it, a real run's actual output, not
## a value out of the manual) has the real answer: 3.10.

seqmd3 <- utils::read.table(system.file("extdata", "exam", "seqm.d3", package = "tdaR"),
                            na.strings = "-1")
names(seqmd3) <- paste0("Y", 1:8)
scost <- matrix(0, 10, 10)
for (i in 1:10) for (j in 1:10) scost[i, j] <- 0.1 * abs(i - j)

m3 <- tda_seqm(seqmd3, indel = 1, subcost = scost)
same("seqm: subcost= as a matrix matches a real run of seqm3.cf (seqm3.d)",
     as.vector(m3), 3.1, 1e-6)
same("seqm: states is the sequence's distinct values, ascending",
     attr(m3, "states"), 1:7)

## TDA treats *any* negative raw value as "no state" (t_seq.c's
## seq_sget/get_sdata: `s < 0`, a plain sign test, not `== -1`) --
## states must exclude it whether the gap arrived as R's NA or as a
## literal negative value straight off disk, seqm.d3's convention if
## `sequences` isn't converted to NA first.
seqmd3_raw <- utils::read.table(system.file("extdata", "exam", "seqm.d3", package = "tdaR"))
names(seqmd3_raw) <- paste0("Y", 1:8)
m3_raw <- tda_seqm(seqmd3_raw, indel = 1, subcost = scost)
same("seqm: states excludes a literal -1 the same as it excludes NA",
     attr(m3_raw, "states"), 1:7)
same("seqm: literal -1 and na.strings='-1' give the identical distance",
     as.vector(m3_raw), as.vector(m3))


ok("seqm: subcost= a too-small matrix errors instead of misreading it",
   inherits(try(tda_seqm(seqmd3, subcost = scost[1:3, 1:3]), silent = TRUE),
            "try-error"))
ok("seqm: subcost= with negative values errors",
   inherits(try(tda_seqm(seqmd3, subcost = -scost), silent = TRUE),
            "try-error"))

## -- tda_seqm dp_matrix=TRUE: TDA's tst=2,3/df= alignment dump -----
## Box 12 of the manual, reproducing seqm3.cf's real tst=2,3 output
## exactly (examples/exam/seqm3.tst, checked against a live run above).

m4 <- tda_seqm(seqmd3, indel = 1, subcost = scost, dp_matrix = TRUE)
dp <- attr(m4, "dp_matrix")
## run= is the TDA run: script.cf has tst=2,3/df= in it, and
## dp.tst is the exact file dp_matrix was parsed from -- checked
## directly rather than trusting the parse alone.
run4 <- attr(m4, "run")
cf4 <- run4$commands
ok("seqm: run='s script.cf actually requested tst=2,3 and df=",
   any(grepl("tst\\s*=\\s*2,3", cf4)) && any(grepl("df\\s*=\\s*dp\\.tst", cf4)))
ok("seqm: run$dir's dp.tst is the exact file dp_matrix was parsed from",
   file.exists(file.path(run4$dir, "dp.tst")))
ok("seqm: dp_matrix finds the one pair actually compared",
   identical(names(dp), "2_1"))
d21 <- dp[["2_1"]]$D
same("seqm: dp_matrix's D bottom-right equals the actual distance",
     d21[nrow(d21), ncol(d21)], as.vector(m4))
same("seqm: dp_matrix's D row/column labels match Box 12 exactly",
     dimnames(d21),
     list(c("A", "1", "1", "2", "3", "4", "5", "6"),
         c("B", "4", "5", "6", "7", "7", "7", "7", "7")))
same("seqm: dp_matrix's D matches Box 12 of the manual exactly",
     unname(d21),
     matrix(c(0, 1, 2, 3, 4, 5, 6, 7, 8,
              1, 0.3, 1.3, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3,
              2, 1.3, 0.7, 1.7, 2.7, 3.7, 4.7, 5.7, 6.7,
              3, 2.2, 1.6, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1,
              4, 3.1, 2.4, 1.9, 1.5, 2.5, 3.5, 4.5, 5.5,
              5, 4.0, 3.2, 2.6, 2.2, 1.8, 2.8, 3.8, 4.8,
              6, 5.0, 4.0, 3.3, 2.8, 2.4, 2.0, 3.0, 4.0,
              7, 6.0, 5.0, 4.0, 3.4, 2.9, 2.5, 2.1, 3.1),
            nrow = 8, byrow = TRUE), 1e-6)
ok("seqm: dp_matrix's seqA/seqB are the real sequences, in order",
   identical(dp[["2_1"]]$seqA, c(1, 1, 2, 3, 4, 5, 6)) &&
   identical(dp[["2_1"]]$seqB, c(4, 5, 6, 7, 7, 7, 7, 7)))
ok("seqm: dp_matrix is NULL by default -- opt-in, not a regular output",
   is.null(attr(m0, "dp_matrix")))

## a multi-pair case: every pair's D bottom-right must equal the main
## distance matrix's entry for that pair.
s5 <- data.frame(t1 = c(1, 1, 2), t2 = c(1, 2, 2), t3 = c(2, 2, 3))
m5 <- tda_seqm(s5, dp_matrix = TRUE)
dp5 <- attr(m5, "dp_matrix")
mm5 <- as.matrix(m5)
ok("seqm: dp_matrix finds all three pairs for 3 cases",
   setequal(names(dp5), c("2_1", "3_1", "3_2")))
for (nm in names(dp5)) {
    ij <- as.integer(strsplit(nm, "_")[[1L]])
    dmat <- dp5[[nm]]$D
    same(paste("seqm: dp_matrix pair", nm, "matches the distance matrix"),
         dmat[nrow(dmat), ncol(dmat)], mm5[ij[1L], ij[2L]])
}

## tda_seqev/tda_seqevd/tda_seqmd against seq.d4, TDA's real example
## data for these commands -- ID, six states (Y0..Y5), six more columns
## (S0..S5), two covariates (V1, V2), the realistic shape (a sequence
## alongside other columns in the same data frame), not the states alone.
##
## tda_seqev had a real, confirmed bug: seqev's syntax requires a
## right-hand-side file name (its results are only ever written to a
## file, never printed to the console at all), which an earlier version
## of this function never supplied, then tried to parse console text
## for that never existed. Checked here against a real run of the
## standalone binary with the identical script, not just that it
## returns something.
seqd4 <- utils::read.table(system.file("extdata", "exam", "seq.d4", package = "tdaR"))
names(seqd4) <- c("ID", paste0("Y", 0:5), paste0("S", 0:5), "V1", "V2")
yv4 <- paste0("Y", 0:5)

ev4 <- tda_seqev(yv4, data = seqd4)
ok("seqev: returns a real from/to/count data frame, not text",
   is.data.frame(ev4) && identical(names(ev4), c("from", "to", "count")))
same("seqev: matches a real run of the standalone binary exactly",
     ev4[order(ev4$from, ev4$to), ],
     data.frame(from = c(1L, 2L), to = c(2L, 1L), count = c(4L, 3L)))
ok("seqev: sequence= naming columns not in data= errors clearly",
   inherits(try(tda_seqev("NotAColumn", data = seqd4), silent = TRUE),
            "try-error"))
ok("seqev: a pre-subsetted data frame (no data=) still works, identically",
   isTRUE(all.equal(tda_seqev(seqd4[yv4]),
                    ev4[order(ev4$from, ev4$to), ], check.attributes = FALSE)))

evd4 <- tda_seqevd(yv4, data = seqd4)
## Column names now come from dtda= directly (ev1_2/ev2_1, naming the
## real transition each counts), not a generic events1/events2 --
## confirmed against a real run of the standalone binary (manual's
## Box 4 example), fixed after finding the earlier generic naming threw
## away real information the dtda= description already had.
ok("seqevd: one row per time point, transition-named event columns",
   is.data.frame(evd4) && nrow(evd4) == 5L &&
   all(c("time", "cases", "ev1_2", "ev2_1", "total") %in% names(evd4)))
same("seqevd: matches a direct run of the standalone binary (Box 4)",
     evd4, data.frame(time = 1:5, cases = rep(3, 5),
                      ev1_2 = c(1, 1, 0, 1, 1), ev2_1 = c(0, 2, 0, 1, 0),
                      total = c(1, 3, 0, 2, 1)))

md4 <- tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                 covariates = c("V1", "V2"))
ok("seqmd: sequence=/covariates= as column names in data= both work together",
   is.data.frame(md4) && all(c("id", "time", "event", "V1", "V2") %in%
                             names(md4)))
same("seqmd: matches the same call with a pre-subsetted data frame exactly",
     md4,
     tda_seqmd(seqd4[yv4], event = c(1, 2),
              covariates = seqd4[c("V1", "V2")]))

## summary = TRUE: a separate mode of seqmd (confirmed
## directly, not a different reading of the same output) -- called
## with no right-hand-side file at all, TDA prints a Time/RiskSet/
## Events table to the console instead of writing any per-case
## records. Real, confirmed values from a direct run of the standalone
## binary on this exact data (manual's Box 1 example).
md_sum <- tda_seqmd(yv4, event = c(1, 2), data = seqd4, summary = TRUE)
same("seqmd: summary = TRUE matches a direct run of the standalone binary",
     md_sum,
     data.frame(time = 1:5, riskset = c(2, 1, 2, 2, 2),
               events = c(1, 1, 0, 1, 1)))
ok("seqmd: summary = TRUE rejects covariates=",
   inherits(try(tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                         summary = TRUE, covariates = "V1"),
                silent = TRUE), "try-error"))
ok("seqmd: summary = TRUE rejects event_covariates=",
   inherits(try(tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                         summary = TRUE,
                         event_covariates = list(price = paste0("S", 0:5))),
                silent = TRUE), "try-error"))

## event_covariates= (xe=): different from covariates=/v= --
## one value per time point, contributing whichever value matches the
## event's time, not a fixed value per case -- checked
## against a real run of the standalone binary with the identical
## script (xe=[S0,S1,...,S5], one bracketed group per covariate, not a
## plain comma list the way v= is) from the name alone.
## This is very likely what seq.d4's S0..S5 columns, otherwise
## completely unused by any of this, are actually for.
md_xe <- tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                   covariates = c("V1", "V2"),
                   event_covariates = list(price = paste0("S", 0:5)))
ok("seqmd: event_covariates= adds a real, named, time-varying column",
   "price" %in% names(md_xe))
ok("seqmd: ... without disturbing the id/time/event/period/covariate names",
   identical(names(md_xe)[1:10],
            c("id", "time", "event", paste0("period", 1:5), "V1", "V2")))
ok("seqmd: event_covariates= produces a real, non-missing numeric column",
   is.numeric(md_xe$price) && !anyNA(md_xe$price) &&
   nrow(md_xe) == length(md_xe$price))
ok("seqmd: event_covariates= without a name for each element errors clearly",
   inherits(try(tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                         event_covariates = list(paste0("S", 0:5))),
                silent = TRUE), "try-error"))

## sequences= as a list: several independent sequence data structures
## at once (TDA's seqdef(sn=1), seqdef(sn=2), ...) -- confirmed a
## real, documented pattern, not a hypothetical one: TDA's manual
## builds a second structure from seq.d4's S0..S5 columns this way,
## states {1,3}, entirely separate from the first (Y0..Y5, states
## {1,2}), not a covariate of it -- and sn= to pick which one a call
## analyzes, checked against a direct run of the standalone binary for
## both structures, not just that sn=2 gives *some* different answer.
multi_seq <- list(yv4, paste0("S", 0:5))
ev_sn1 <- tda_seqev(multi_seq, data = seqd4)
ev_sn1_explicit <- tda_seqev(multi_seq, data = seqd4, sn = 1)
ev_sn2 <- tda_seqev(multi_seq, data = seqd4, sn = 2)
ok("seqev: sn= defaults to the first structure, unchanged from before",
   isTRUE(all.equal(no_run(ev_sn1), no_run(ev_sn1_explicit))))
same("seqev: sn=1 (Y0..Y5) matches the direct binary run",
     ev_sn1[order(ev_sn1$from, ev_sn1$to), ],
     data.frame(from = c(1L, 2L), to = c(2L, 1L), count = c(4L, 3L)))
same("seqev: sn=2 (S0..S5) is a different, real structure",
     ev_sn2[order(ev_sn2$from, ev_sn2$to), ],
     data.frame(from = c(1L, 3L), to = c(3L, 1L), count = c(5L, 2L)))
ok("seqev: sn= out of range errors clearly",
   inherits(try(tda_seqev(multi_seq, data = seqd4, sn = 3), silent = TRUE),
            "try-error"))

evd_sn2 <- tda_seqevd(multi_seq, data = seqd4, sn = 2)
ok("seqevd: sn=2's column names reflect its states (1,3), not (1,2)",
   all(c("ev1_3", "ev3_1") %in% names(evd_sn2)))

## select= is a raw string, never translated the way tda_frml's
## definitions= is -- so a real column name in it (Y0, the way the user
## actually named their data) has to be rewritten here to whatever
## internal name TDA was actually given, the same problem and fix as
## Surv()'s arguments in tda_frml. This was a genuine, confirmed
## regression from the sequences=-as-a-list work above: the internal
## names changed from Y1..Y6 to SQ1_1..SQ1_6 for multi-structure
## support, breaking select= silently until fixed here.
sel_ev <- tda_seqev(yv4, data = seqd4, select = "eq(Y0,1)")
sel_ev_unrestricted <- tda_seqev(yv4, data = seqd4)
ok("seqev: select= with the real column name restricts the result",
   !isTRUE(all.equal(sel_ev[order(sel_ev$from, sel_ev$to), ],
                     sel_ev_unrestricted[order(sel_ev_unrestricted$from,
                                              sel_ev_unrestricted$to), ])))
same("seqev: select= gives the exact, real restricted counts",
     sel_ev[order(sel_ev$from, sel_ev$to), ],
     data.frame(from = c(1L, 2L), to = c(2L, 1L), count = c(3L, 2L)))
sel_ev_multi <- tda_seqev(multi_seq, data = seqd4, sn = 2,
                          select = "eq(S0,1)")
ok("seqev: select= translation also works together with sn= on a
   second structure",
   is.data.frame(sel_ev_multi) && nrow(sel_ev_multi) > 0L)

sel_md <- tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                    select = "eq(Y0,1)")
ok("seqmd: select= with the real column name also works",
   nrow(sel_md) == 3L)
sel_md_sum <- tda_seqmd(yv4, event = c(1, 2), data = seqd4,
                        select = "eq(Y0,1)", summary = TRUE)
ok("seqmd: select= works in summary = TRUE mode too",
   is.data.frame(sel_md_sum) && nrow(sel_md_sum) == 5L)

## select= also accepts R's comparison operators (==, and a bare =
## too, since that isn't valid R syntax to pass unquoted as an argument
## and so has to stay reachable through the string form) instead of
## requiring TDA's eq(Y0,1) function-call form up front -- a plain
## text substitution, not R's parser, checked here against the same
## real, restricted counts as the raw TDA-syntax form already verified
## above, and that TDA-syntax form itself still passes through
## untouched.
sel_eq2 <- tda_seqev(yv4, data = seqd4, select = "Y0 == 1")
sel_eq1 <- tda_seqev(yv4, data = seqd4, select = "Y0 = 1")
same("seqev: select= with == matches the raw TDA-syntax eq() form exactly",
     sel_eq2[order(sel_eq2$from, sel_eq2$to), ],
     sel_ev[order(sel_ev$from, sel_ev$to), ])
same("seqev: select= with a bare = also matches",
     sel_eq1[order(sel_eq1$from, sel_eq1$to), ],
     sel_ev[order(sel_ev$from, sel_ev$to), ])
sel_and <- tda_seqev(yv4, data = seqd4, select = "Y0 == 1 & Y1 == 1")
ok("seqev: select= with == combined by & further restricts the result",
   sum(sel_and$count) < sum(sel_eq2$count))
gc_all <- tda_seqgc(seqd4[yv4])
gc_sel <- tda_seqgc(seqd4[yv4], select = "Y0 == 1")
ok("seqgc: select= with a real column name and == also works (.seq_desc path)",
   nrow(gc_sel) < nrow(gc_all))

## tda_seq_info: TDA's "seq;" -- the same Structure/Type/Variables/
## Time-axis/Number-of-States/States table TDA itself prints after
## defining a sequence structure (manual's seq8.cf, Box 3's "upper
## table"). Verified directly against a real run of the standalone
## binary with the identical script, not just that it returns
## something -- and against the same multi-structure data (Y0..Y5 as
## sn=1, S0..S5 as a separate sn=2) already established
## above for tda_seqev's sn= support.
info1 <- tda_seq_info(yv4, data = seqd4)
ok("seq_info: single structure, one row, real fields",
   is.data.frame(info1) && nrow(info1) == 1L &&
   all(c("sn", "type", "variables", "tmin", "tmax", "nstates",
        "states") %in% names(info1)))
same("seq_info: matches a direct run of the standalone binary (single structure)",
     info1[c("sn", "type", "variables", "tmin", "tmax", "nstates")],
     data.frame(sn = 1L, type = 1L, variables = 6L, tmin = 0, tmax = 5,
               nstates = 2L))
ok("seq_info: states is the real state values, not just a count",
   identical(info1$states[[1L]], c(1L, 2L)))

info2 <- tda_seq_info(multi_seq, data = seqd4)
ok("seq_info: several structures at once, one row each",
   nrow(info2) == 2L)
same("seq_info: sn=2's row reflects its real states (1,3), not (1,2)",
     info2$states[[2L]], c(1L, 3L))
ok("seq_info: matches tda_seqev's sn=2 states exactly",
   identical(info2$states[[2L]], sort(unique(c(ev_sn2$from, ev_sn2$to)))))


## -- seqpd: the four layouts -------------------------------------------
##
## Values checked against the manual's section 3.4.3, Boxes 2 and 3,
## which run seqpd over this same seq.d1: case 1's sequence is
## 7 3 3 1 3 3 7 7, and as spells that is 7->3 over 0-1, 3->1 over 1-3,
## 1->3 over 3-4, 3->7 over 4-6, 7->-1 over 6-8.

pd1 <- tda_seqpd(seqd1, variables = yvars, id = "ID")
ok("seqpd: one row per case", nrow(pd1) == nrow(seqd1))
same("seqpd: case 1's sequence, in order",
     unname(unlist(pd1[1L, grep("^y", names(pd1))])),
     c(7, 3, 3, 1, 3, 3, 7, 7))
same("seqpd: id= carries the case's value", pd1$id, seqd1$ID)

pd3 <- tda_seqpd(seqd1, layout = 3, variables = yvars, id = "ID")
ok("seqpd layout 3: one row per case and position",
   nrow(pd3) == nrow(seqd1) * length(yvars))

pd4 <- tda_seqpd(seqd1, layout = 4, variables = yvars, id = "ID")
same("seqpd layout 4: case 1's spells, origin and destination",
     unname(pd4$org[pd4$case == 1]), c(7, 3, 1, 3, 7))
same("seqpd layout 4: and the times they span",
     unname(pd4$ts[pd4$case == 1]), c(0, 1, 3, 4, 6))

ok("seqpd: a layout outside 1:4 is refused, not passed through",
   inherits(try(tda_seqpd(seqd1, layout = 9, variables = yvars),
                silent = TRUE), "try-error"))

## -- seqdef's two data types (manual 3.4.2, Box 6) ----------------------
##
## The same six columns of seq.d2 read as type 1 are six variables over a
## 0-5 time axis with nine states; read as type 2 -- state/time pairs --
## they are three variables over 0-8 with three states. Only type 1 was
## reachable before.

seqd2 <- utils::read.table(system.file("extdata", "exam", "seq.d2", package = "tdaR"),
                           col.names = c("ID", "Y1", "T1", "Y2", "T2",
                                         "Y3", "T3"))
i1 <- tda_seq_info(seqd2[, -1])
i2 <- tda_seq_info(seqd2[, -1], type = 2)
same("seqdef type 1: six variables", i1$variables, 6)
same("seqdef type 2: three", i2$variables, 3)
same("seqdef type 2: the time axis comes from the times", i2$tmax, 8)
same("seqdef type 2: and three states, not nine", i2$nstates, 3)
