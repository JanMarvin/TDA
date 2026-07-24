if (have_examples) {

## rdataframe -- the data matrix built from R's columns rather than written
## out as text and parsed back.  The numbers have to be the ones TDA gets
## from its file, so ehd2 is the check: same log likelihood, same
## coefficients.
rdf <- tda_rrdat()
rdfit <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, rdf)
ok("rdataframe: the command is what runs",
   any(grepl("^rdataframe;$", rdfit$run$commands)))
same("rdataframe: ehd2 through memory", as.numeric(stats::logLik(rdfit)),
     -2465.9873, 5e-4)
same("rdataframe: and its coefficients",
     unname(round(coef(rdfit), 4)),
     c(-4.4894, 0.0773, 0.6080, 0.6108, -0.0032, 0.0596, -0.0280), 5e-5)

## NA becomes the code rstata() and rspss1() use, ctx->PMMSYS, and says so.
rdna <- data.frame(t = c(10, 20, NA, 40, 50, 60), s = c(1L, 1L, 1L, 0L, 1L, 1L),
                   x = c(1, 2, 3, NA, 5, 6))
rna <- tda_run(c(tda_nvar(rdna),
                 tda_block("edef", ts = 0, tf = "t", org = 0, des = "s")),
               data = rdna)
ok("rdataframe: missing values are substituted and reported",
   any(grepl("2 missing values were substituted by: -5", rna$output)))

## NA must never reach TDA as a number.  It is written as the point marker,
## which nvar() counts and mpnt= gives a value to -- the same value the
## in-memory path substitutes, so the two agree.
nad <- data.frame(Dur = c(10, 20, NA, 40, NaN, 60),
                  Des = c(1L, 1L, 1L, 0L, 1L, NA),
                  X = c(1, 2, 3, Inf, 5, 6))
naf <- tempfile()
tda_write_data(nad, naf)
ok("write_data: every kind of NA becomes the point marker",
   sum(vapply(strsplit(readLines(naf), " "), function(r) sum(r == "."),
              integer(1))) == 4L)
nar <- tda_run(c(tda_nvar(nad, file = "data.dat"),
                 tda_block("edef", ts = 0, tf = "Dur", org = 0, des = "Des")),
               data = nad)
ok("nvar: they are counted as missing, not read as data",
   any(grepl("^Point +4 +-5", nar$output)))

## The Excl column of the episode table: "*" marks episodes belonging to no
## transition, which TDA leaves out of every estimate and warns about.  An
## origin state with no way out of it produces one.
exd <- tda_rrdat()
exd$ORG <- ifelse(exd$NOJ > 3L, 1L, 0L)
exd$DES2 <- ifelse(exd$ORG == 1L, 1L, exd$DES)
exr <- tda_run(c(tda_nvar(exd),
                 tda_block("edef", ts = 0, tf = "TFP", org = "ORG",
                           des = "DES2")), data = exd)
extab <- tdaR:::.tda_episode_table(exr)
ok("episodes: the Excl column comes back",
   "excluded" %in% names(extab))
ok("episodes: and marks the ones TDA drops",
   identical(extab$excluded, c(FALSE, FALSE, TRUE)))

## The spatial round trip: sdgen writes the objects, sdnvar reads them back.
## SDID, SDTyp and SDN are columns of the object record; only SDPtr is the
## dummy, and `rd` is a uniform random deviate rather than "read".
sdd <- tda_spatial(data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1)))
sdi <- tda_sd_info(sdd)
ok("sdnvar: the objects come back",
   any(grepl("^Number of points: 4", sdi$run$output)))
ok("sdnvar: with the right bounding box",
   any(grepl("XMax: *1\\.0+ *YMax: *1\\.0+", sdi$run$output)))
for (sdnm in c("intersect", "lines", "points", "data")) {
    ok(sprintf("spatial: tda_sd_%s runs", sdnm),
       !inherits(try(get(paste0("tda_sd_", sdnm))(sdd), silent = TRUE),
                 "try-error"))
}
ok("spatial: tda_sd_clip takes the rectangle as numbers",
   !inherits(try(tda_sd_clip(sdd, c(0, 0, 1, 1)), silent = TRUE), "try-error"))
ok("spatial: and refuses to guess one",
   inherits(try(tda_sd_clip(sdd), silent = TRUE), "try-error"))
## The three-dimensional plots need psetup3 -- pxa, pya, pza and no pylen --
## and either a marker or a zvar, or TDA draws nothing at all by design.
sd3 <- tda_spatial(data.frame(id = 1:8, x = seq(0, 1, length.out = 8),
                              y = seq(1, 0, length.out = 8),
                              z = seq(0, 1, length.out = 8)),
                   attributes = "z")
ok("spatial: tda_sd_plot3 draws in three dimensions",
   !inherits(try(tda_sd_plot3(sd3, z = c(0, 1), symbol = 1, size = 3),
                 silent = TRUE), "try-error"))
ok("spatial: and says so rather than drawing nothing",
   inherits(try(tda_sd_plot3(sd3, z = c(0, 1), symbol = NULL), silent = TRUE),
            "try-error"))
ok("spatial: z has no default to fall back on",
   inherits(try(tda_sd_plot3(sd3), silent = TRUE), "try-error"))

ok("sdvd: a Voronoi diagram of them",
   any(grepl("Number of Delaunay triangles: 2", tda_sd_voronoi(sdd)$run$output)))

## The rest of the graph family.  Each wrapper is only the command name and
## the table TDA writes, so what is worth checking is that every one runs on
## a graph of the kind it needs -- TDA refuses the others rather than giving
## a wrong answer, which is how the requirements were found.
ge <- data.frame(from = c(1, 1, 2, 2, 3, 4), to = c(2, 3, 3, 4, 4, 5),
                 value = 1)
gdd <- tda_graph(ge, directed = TRUE)
gdu2 <- tda_graph(ge, directed = FALSE)
gwants <- c(aggregate = 1, symmetric = 1, toposort = 1, nodes = 1, dot = 1,
            links = 1, dcycles = 1, dblocks = 1, random = 1, ownership = 1,
            backward = 1, paths = 1, flowcontrol = 1, reachable = 1,
            compact = 0, gcliques = 0, centred = 0, spantrees = 0)
for (gnm in names(gwants)) {
    gg <- if (gwants[[gnm]]) gdd else gdu2
    ok(sprintf("graph: tda_g_%s runs", gnm),
       !inherits(try(get(paste0("tda_g_", gnm))(gg), silent = TRUE),
                 "try-error"))
}
## gcni is the exception: it needs an unvalued graph, and tda_graph() always
## builds a valued one, so the wrapper sets gt itself.
ok("graph: tda_g_neighbourhoods makes the graph unvalued for gcni",
   any(grepl("undirected, unvalued",
             tda_g_neighbourhoods(
                 tda_graph(ge[c("from", "to")], directed = FALSE))$run$output)))

## TDA's manual, read from the copy in inst/extdata.  help() paginates
## on stdin, which would hang a run under R, so hlp_inp() returns "next
## page" when the package is built.
hlp <- utils::capture.output(hres <- tda_help("rate"))
hlp <- hres
ok("help: the manual entry comes back",
   any(grepl("The rate command estimates transition rate models", hlp)))
ok("help: without TDA's banner", !any(grepl("^TDA\\. Analysis", hlp)))
ok("help: an empty topic lists the keywords",
   sum(lengths(strsplit(
       utils::capture.output(htop <- tda_help("")), "[[:space:]]{2,}"))) > 150L)
ok("time: reports a time", any(grepl("^Current time:", utils::capture.output(tda_time()))))
ok("mem: reports bytes", any(grepl("bytes", utils::capture.output(tda_memory()))))

## sdshp, against a real ESRI shapefile if one has been provided.  nc.shp
## from the sf package is the usual one: 100 North Carolina counties, of
## which 6 have more than one part.
if (have_ext("nc.shp", "nc.shx", "nc.dbf")) {
    ncs <- tda_read_shapefile(file.path(EXT, "nc.shp"))
    ok("sdshp: reads all hundred shapes",
       any(grepl("Number of shapes: 100", ncs$run$output)))
    ok("sdshp: and the multi-part ones become separate polygons",
       any(grepl("Object-description records: 108", ncs$run$output)))
    ok("sdshp: the result is a spatial structure the rest can use",
       any(grepl("Number of polygons: 108", tda_sd_info(ncs)$run$output)))
    ok("sdshp: refuses a stem without its index and attributes",
       inherits(try(tda_read_shapefile(file.path(EXT, "nowhere.shp")),
                    silent = TRUE), "try-error"))

    ## A map: psetupg builds the projection, sdpmap draws the objects,
    ## sdpgrat the graticule and sdpgeo the marked points.
    ok("map: all hundred and eight polygons are projected and drawn",
       !inherits(try(tda_map(ncs, view = c(-80, 35), region = c(6, 3),
                             graticule = list(lon = "-84(2)-76",
                                              lat = "34(1)37"),
                             points = data.frame(lon = -78.64, lat = 35.78)),
                     silent = TRUE), "try-error"))
    ok("plot3: polygons draw flat without a zvar",
       !inherits(try(tda_sd_plot3_polygons(ncs, z = c(0, 1),
                                           xlim = c(-84.4, -75.4),
                                           ylim = c(33.8, 36.6)),
                     silent = TRUE), "try-error"))
    ok("plot3: and a zvar raises them into prisms",
       !inherits(try(tda_sd_plot3_polygons(ncs, z = c(0, 0.25),
                                           xlim = c(-84.4, -75.4),
                                           ylim = c(33.8, 36.6),
                                           zvar = "AREA", gs = 0.85),
                     silent = TRUE), "try-error"))
    ok("map: region has no default to fall back on",
       inherits(try(tda_map(ncs, view = c(-80, 35)), silent = TRUE),
                "try-error"))
}

## sdcpol: the four sides of a square, reassembled into the square.
sqseg <- data.frame(x1 = c(0, 1, 1, 0), y1 = c(0, 0, 1, 1),
                    x2 = c(1, 1, 0, 0), y2 = c(0, 1, 1, 0))
sqpol <- tda_polygons(sqseg)
ok("sdcpol: the segments become one polygon",
   any(grepl("Number of polygons: 1", tda_sd_info(sqpol)$run$output)))
ok("sdcpol: option 1 gives the node list instead",
   nrow(tda_polygons(sqseg, option = 1)$table) == 4L)
ok("sdcpol: four columns are required",
   inherits(try(tda_polygons(sqseg[1:2]), silent = TRUE), "try-error"))

## Projection 20 (orthographic).  The map viewport is a rectangle in
## projected coordinates, and the projected x of the region edge uses
## sin(lonD), which folds back past 90 degrees: TDA used to accept a
## longitude half-width in (90, 180] and either draw the map for
## 180 - lonD instead (150 drew 30's map) or emit coordinates around
## 1e18 (region 180,80).  lonD > 90 is refused now; corners of a
## rectangle that spills off the projection disc (60,60) were always
## refused, correctly.
pdf(NULL)
p20 <- function(reg) try(tda_map(sqpol, view = c(0, 0), region = reg,
                                 projection = 20,
                                 graticule = list(lon = 30, lat = 30)),
                         silent = TRUE)
ok("proj 20: a region inside the disc draws",
   !inherits(p20(c(45, 45)), "try-error"))
r <- p20(c(150, 40))
ok("proj 20: the fold-back band is refused, not silently remapped",
   inherits(r, "try-error") &&
   grepl("region parameter", attr(r, "condition")$message))
ok("proj 20: a half-width of 180 is refused",
   inherits(p20(c(180, 80)), "try-error"))
ok("proj 20: a rectangle whose corners leave the disc is refused",
   inherits(p20(c(60, 60)), "try-error"))
dev.off()

## The other readers have no test data, so all that can be checked is that
## they refuse a file that is not there rather than failing later.
for (rd in c("tda_read_e00", "tda_read_gshhs", "tda_read_dcw")) {
    ok(sprintf("%s: refuses a missing file", rd),
       inherits(try(get(rd)("nowhere.dat"), silent = TRUE), "try-error"))
}

## lsreg1 names its independent variables on the right-hand side, with the
## response given by yl= -- leaving the response in the varlist regresses it
## on itself, which fits perfectly and returns 1 on the response and 0 on
## everything else.  And it writes its estimates to ppar= rather than
## printing them, so the wrapper asks for that file and reads it back.
set.seed(1)
cd <- data.frame(x = 1:40, z = rnorm(40))
cd$y <- 2 + 0.5 * cd$x + cd$z
cd$cen <- 0L
cd$yc <- cd$y
same("lsreg1 with nothing censored reproduces least squares",
     unname(round(coef(tda_lsreg1(yc ~ x + z, cd, censor = "cen")), 6)),
     unname(round(coef(tda_lsreg(y ~ x + z, cd)), 6)), 1e-6)
ok("lsreg1 returns its coefficients without being asked",
   length(coef(tda_lsreg1(yc ~ x + z, cd, censor = "cen"))) == 3L)

## print() and summary() have to survive a fit with no log likelihood and a
## coefficient table with no standard errors -- a regression has neither, and
## both methods used to fail on it: `if (!is.na(NULL))` is an error, and
## printCoefmat wants an estimate, a standard error and a test statistic.
lsf <- tda_lsreg(y ~ x + z, cd)
ok("print works on a fit with no logLik",
   !inherits(try(utils::capture.output(print(lsf)), silent = TRUE),
             "try-error"))
ok("summary works on a table with no standard errors",
   !inherits(try(utils::capture.output(print(summary(lsf))), silent = TRUE),
             "try-error"))
ok("summary still prints the full table when there is one",
   any(grepl("C/Error",
             utils::capture.output(print(summary(
                 tda_rate(Surv(TFP, DES) ~ EDU, tda_rrdat())))))))
## "---" is TDA's "no value" and becomes NA in the parsed table, including
## when a whole column is "---" -- an exact fit leaves no residual variance,
## so every standard error is one.  Before, the column stayed character and
## the fallback in vcov() did "---"^2.
cd$yp <- 2 + 0.5 * cd$x + cd$z
lsx <- tda_lsreg(yp ~ x + z, cd)
ok("--- becomes NA even when the whole column is ---",
   all(is.na(lsx$estimates$Error)) && is.numeric(lsx$estimates$Error))
ok("and vcov() survives it",
   !inherits(try(suppressWarnings(vcov(lsx)), silent = TRUE), "try-error"))
## The rest needs a fit with residual variance: cd's response is an exact
## function of its predictors, which is the case the two checks above are
## about.
cd$yn <- cd$y + rnorm(nrow(cd))
lsn <- tda_lsreg(yn ~ x + z, cd)

## summary() also keeps what TDA reports about the fit itself, which used to
## be thrown away: R squared, the F statistic, the residual variance and the
## degrees of freedom, so the report reads like lm's.
lso <- utils::capture.output(print(summary(lsn)))
ok("summary reports the residual standard error and df",
   any(grepl("Residual standard error: .* on 37 degrees of freedom", lso)))
ok("summary reports R squared and the F statistic",
   any(grepl("Multiple R-squared", lso)) && any(grepl("F-statistic", lso)))
ok("a model summary is unaffected: it has no such lines",
   !any(grepl("R-squared", utils::capture.output(print(summary(
       tda_rate(Surv(TFP, DES) ~ EDU, tda_rrdat())))))))

ok("lsreg's Coeff/E is the same column as the models' C/Error",
   "C/Error" %in% names(lsn$estimates))
ok("vcov() is the matrix TDA wrote, not a reconstructed diagonal",
   { v <- vcov(lsn); is.matrix(v) && nrow(v) == 3L &&
     any(abs(v[upper.tri(v)]) > 0) })

## Censoring, which is easy to get backwards: TDA's cen= marks the cases it
## can use -- gdf_dcheck() calls a case right censored when the indicator is
## zero -- while an R user writes 1 for a censored case. tda_lsreg1()
## translates, so the count TDA reports must match the count marked.
ce <- data.frame(Yy = as.numeric(1:10), Xx = as.numeric(1:10),
                 Cc = as.integer((1:10) > 7))
cef <- tda_lsreg1(Yy ~ Xx, ce, censor = "Cc")
ok("lsreg1: three marked censored are three censored to TDA",
   any(grepl("Number of right censored cases: 3", cef$run$output)) &&
   any(grepl("Number of exact cases: 7", cef$run$output)))
ok("lsreg1: and none marked means none censored",
   any(grepl("Number of right censored cases: 0",
             tda_lsreg1(Yy ~ Xx, transform(ce, Cc = 0L),
                        censor = "Cc")$run$output)))

## A character `data` is a file for TDA to parse, not something R reads first.
rfnames <- c("ID", "NOJ", "TStart", "TFin", "SEX", "TI", "TB", "TE", "TMAR",
             "PRES", "PRESN", "EDU")
rfd <- stats::setNames(as.data.frame(matrix(0L, 1L, 12L)), rfnames)
# noc comes from the frame's rows, and this one is only a carrier for the
# column names, so the file's case count is given explicitly.
rfile <- tda_run(c(tda_nvar(rfd, file = "rrdat.1", noc = 600), "mem;"),
                 data = file.path(EX, "ehhnew", "rrdat.1"))
ok("a data file is left for TDA to read",
   any(grepl("Number of cases: 600", rfile$output)))

# stringsAsFactors: the package fixes the default once in its own
# namespace rather than at 250 call sites, so it behaves the same on any
# R version.  R 4.0 changed data.frame()'s default from TRUE to FALSE,
# and nearly every data.frame() here relied on it -- on R 3.x each
# character column would quietly have become a factor.
ok("data.frame inside the package defaults to character",
   is.character(get("data.frame",
                    envir = asNamespace("tdaR"))(a = c("x", "y"))$a))
ok("as.data.frame inside the package defaults to character",
   is.character(get("as.data.frame",
                    envir = asNamespace("tdaR"))(list(a = c("x", "y")))$a))
# and an explicit value still wins -- passing it twice was the first
# version of this shim's bug ("matched by multiple actual arguments")
ok("an explicit stringsAsFactors is still honoured",
   is.factor(get("as.data.frame", envir = asNamespace("tdaR"))(
       list(a = c("x", "y")), stringsAsFactors = TRUE)$a))

## A character or factor column is now CARRIED, as a TDA string
## variable: it cannot go through rdataframe (that path hands TDA the
## columns already in memory and has no type for one), so tda_nvar()
## falls back to a file and declares it with `= str(n,m)`, TDA's only
## way of making one.  A bare "rdataframe;" therefore still cannot take
## it -- there is no file for the string to live in.
# rdataframe hands TDA the columns already in memory and save_var() has
# no string branch, so it refuses -- in TDA's output, not as an R
# condition, because the caller wrote the command by hand.
ok("rdataframe alone still cannot carry a factor column",
   any(grepl("column f is not numeric",
             tda_run("rdataframe;",
                     data = data.frame(f = factor("a")))$output)))
local({
    sd_ <- data.frame(ID = 1:3, S = c("pear", "fig", "date"),
                      X = c(1.5, 2.5, 3.5), stringsAsFactors = FALSE)
    ok("tda_nvar declares a character column with str(n,m)",
       any(grepl("S = str\\([0-9]+,[0-9]+\\)", tda_nvar(sd_))))
    r <- tda_run(c(tda_nvar(sd_), "pdata() = out.txt;"), data = sd_)
    ok("a character column round-trips as character",
       identical(tda_strings(r)$S, sd_$S))
    # a factor is its labels, the reading .tda_columns() already took
    fd_ <- sd_
    fd_$S <- factor(fd_$S)
    rf <- tda_run(c(tda_nvar(fd_), "pdata() = out.txt;"), data = fd_)
    ok("a factor column comes back as its labels, not its codes",
       identical(tda_strings(rf)$S, as.character(fd_$S)))
})
# Strings are written LAST in the data file whatever order the caller
# used: TDA cannot reach a `= cK` field lying after a string one, and
# fails with "Can't read c3 in data file" and no data matrix at all.
# Declaration order is independent of file order, so the caller's
# column order survives.
# tda_help reads the section from tda.hlp rather than re-splitting
# TDA's printed output.  TDA's help printer JOINS logical lines and pads
# them apart -- one line of the `rate` entry arrives 816 characters long
# holding a dozen real lines -- and no re-splitting recovers what each
# was indented by, because the padding replaced it.
local({
    h <- tda_help("rate")
    ok("help: the option table keeps its indentation",
       any(grepl("^        tp=", h)))
    ok("help: a continuation note stays at its deeper indent",
       any(grepl("^                        in addition", h)))
    ok("help: the options after that note are NOT flushed left",
       any(grepl("^        mplog=", h)))
    ok("help: the closing bracket sits where it was written",
       any(grepl("^    \\) = model_number;", h)))
    ok("help: no line is a joined run of several",
       max(nchar(h)) < 200)
})
ok("help: every section of the file resolves from the file",
   local({
       f <- system.file("extdata", "tda.hlp", package = "tdaR")
       nm <- trimws(sub("^##", "", grep("^##", readLines(f, warn = FALSE),
                                        value = TRUE)))
       all(vapply(nm, function(t) !is.null(tdaR:::.hlp_section(f, t)), NA))
   }))
# section names carry trailing spaces in the file ("##nvar   "), which
# is what made nvar miss and fall through to the joined-line fallback
ok("help: a padded section heading still matches",
   any(grepl("^       vdef,", tda_help("nvar"))))
ok("help: an unknown topic still reports through TDA itself",
   any(grepl("No matching entry", tda_help("no-such-topic-xyz"))))

# The four string operators, each against its R equivalent.  extra=
# takes a { } block of plain R, so they read as R rather than TDA text.
local({
    sd2 <- data.frame(ID = 1:4, S = c("pear", "fig", "date", "kiwi"),
                      N = c("0012", "0345", "0007", "1234"),
                      stringsAsFactors = FALSE)
    dd <- tempfile("strops")
    r <- tda_run(c(tda_nvar(sd2, extra = { Len = strlen(S)
                                           Rank = strsp(S)
                                           Num = strv(N)
                                           Sub = strvp(N, 2, 3) }),
                   "pdata(fmt=24.16) = out.txt;"), data = sd2, dir = dd)
    t <- utils::read.table(file.path(dd, "out.txt"), fill = TRUE)
    names(t) <- c("ID", "S", "N", "Len", "Rank", "Num", "Sub")
    ok("strsp() is the alphabetical position, as R's rank()",
       identical(as.numeric(t$Rank), as.numeric(rank(sd2$S))))
    ok("strv() converts a digit string as as.numeric() does",
       identical(as.numeric(t$Num), as.numeric(sd2$N)))
    ok("strvp(N,2,3) is substr(N,2,3) as a number",
       identical(as.numeric(t$Sub), as.numeric(substr(sd2$N, 2, 3))))
    ok("strlen() is the STORAGE width, not the visible characters",
       all(t$Len == 4L) && min(nchar(sd2$S)) == 3L)
})

ok("a numeric column after a string one still works",
   local({
       d2 <- data.frame(A = 1:3, S = c("x", "yy", "zzz"), B = c(9, 8, 7),
                        stringsAsFactors = FALSE)
       r2 <- tda_run(c(tda_nvar(d2), "pdata() = out.txt;"), data = d2)
       any(grepl("3 records with 3 variables", r2$output))
   }))


# A rate model driven entirely from R: the data arrives as a data frame rather
# than a file, and the reference numbers come from examples/ehhnew/ehg1.ref.
rr <- read.table(file.path(EX, "ehhnew/rrdat.1"))
names(rr) <- c("ID", "NOJ", "TStart", "TFin", "SEX", "TI", "TB", "TE",
               "TMAR", "PRES", "PRESN", "EDU")

cmds <- c(
    tda_nvar(rr,
             fmt = c("3.0", "2.0", "3.0", "3.0", "2.0", "3.0", "3.0", "3.0",
                     "3.0", "3.0", "3.0", "2.0"),
             extra = c("DES [1.0] = if eq(TFin,TI) then 0 else 1",
                       "TFP [3.0] = TFin - TStart + 1")),
    tda_block("edef", ts = 0, tf = "TFP", org = 0, des = "DES"),
    tda_block("rate", "prate(tab=0(1)300) = rates.prs", rhs = 6))

r <- tda_run(cmds, data = rr)

ok("gompertz: run succeeded", r$status == 0)
est <- tda_estimates(r)
ok("gompertz: two estimates", is.data.frame(est) && nrow(est) == 2)
# tda_estimates() now rebuilds these from the coefficient export, so
# they are the doubles TDA computed rather than its 4-decimal print.
# eq()'s tolerance is RELATIVE, which a small coefficient fails on a
# difference the printed reference cannot even express: -0.0066922988
# against -0.0067 is 7.7e-6 absolute but 1.2e-3 relative.  near()
# applies the window the printed value actually justifies.
near("gompertz: baseline constant", est$Coeff[1], -4.0729, decimals = 4L)
near("gompertz: shape constant", est$Coeff[2], -0.0067, decimals = 4L)
near("gompertz: standard error", est$Error[1], 0.0633, decimals = 4L)
eq("gompertz: log likelihood", as.numeric(logLik(r)), -2474.5059, 1e-4)
ok("gompertz: coef names disambiguated",
   identical(names(coef(r)), c("B Constant", "C Constant")))

rates <- tda_file(r, "rates.prs")
ok("prate: 301 rows", nrow(rates) == 301)
ok("prate: columns named from the comment header",
   identical(names(rates), c("ID", "Time", "Surv.F", "Density", "Rate")))
eq("prate: survivor starts at 1", rates$Surv.F[1], 1.0, 1e-6)
ok("prate: survivor decreasing", all(diff(rates$Surv.F) <= 0))


# The same model run from the shipped command file must give the same answer,
# so the data-frame path is not quietly doing something else.
d <- file.path(tempdir(), "cf")
dir.create(d, showWarnings = FALSE)
invisible(file.copy(list.files(file.path(EX, "ehhnew"), full.names = TRUE), d))
rf <- tda_run_cf(file.path(d, "ehg1.cf"))
eq("data frame path matches command file path",
   coef(r), coef(rf), 1e-12)


# Parsing the three estimate table layouts TDA emits.
lay <- file.path(tempdir(), "lay")
dir.create(lay, showWarnings = FALSE)
invisible(file.copy(list.files(file.path(EX, "exam"), full.names = TRUE), lay))

q <- tda_run_cf(file.path(lay, "qr7.cf"))
eq2 <- tda_estimates(q)
ok("qreg: Cat/Term layout parsed", is.data.frame(eq2) && nrow(eq2) == 4)
ok("qreg: name with spaces kept whole", eq2$Variable[2] == "Sigma 2, 1")
ok("qreg: fixed parameter is NA", is.na(eq2$Signif[3]))

f <- tda_run_cf(file.path(lay, "cd2.cf"))
ef <- tda_estimates(f)
ok("fml: Parameter/Value layout parsed", is.data.frame(ef) && nrow(ef) == 11)
ok("fml: Value renamed to Coeff", "Coeff" %in% names(ef))
eq("fml: log likelihood", as.numeric(logLik(f)), -67.8354, 1e-4)


# Repeated runs in one process must not drift: the state that used to leak
# between them is what the context threading and tda_reset_globals() fixed.
c1 <- coef(tda_run(cmds, data = rr))
c2 <- coef(tda_run(cmds, data = rr))
eq("re-entrancy: third run equals first", c1, coef(r), 1e-12)
eq("re-entrancy: fourth run equals third", c2, c1, 1e-12)


# A TDA error must come back as an R condition, not take the session down.
bad <- try(tda_run("nonsense_command(;"), silent = TRUE)
ok("malformed input survives", inherits(bad, "try-error") ||
                               inherits(bad, "tda_result"))
ok("session alive afterwards", is.data.frame(tda_estimates(tda_run(cmds, data = rr))))

}


## tda_rng()/tda_runif()/tda_rnorm(): bit-exact against random1()/normal()
## themselves, compiled from this package's t_rand.c and run directly
## (a standalone C driver, RD1Seed=13421773, RD1I=1, s_normal_dflg=0 --
## TDA's real defaults, confirmed in tda_context_init.c) -- not
## inferred from the algorithm description alone. The reference values
## below came from that driver, printed to 15 decimal digits; the
## comparisons allow 1e-14 for exactly that reason, not to paper over a
## real mismatch (bumping the driver to full double precision and
## rerunning showed agreement to the last bit).

ref_x <- c(0.029103830456734, 0.462388113141060, 0.962853565812111,
          0.278095975518227, 0.049923494458199)
ref_e <- c(-0.461167522868052, -1.599802987375197, 1.058071537208986,
          0.076249701213566, -1.351163065961058)

gen <- tda_rng()
x5 <- e5 <- numeric(5)
for (i in 1:5) {
    x5[i] <- gen$rd()
    e5[i] <- gen$rdn()
}
same("tda_rng: interleaved $rd()/$rdn() matches random1()/normal() exactly",
     x5, ref_x, 1e-14)
same("tda_rng: ... and the normal deviates too, including the cached pair",
     e5, ref_e, 1e-14)

same("tda_runif: matches a fresh tda_rng()$rd() run of the same length",
     tda_runif(5), tda_rng()$rd(5))
gen2 <- tda_rng()
same("tda_rnorm: matches a fresh tda_rng()$rdn() run of the same length",
     tda_rnorm(5), gen2$rdn(5))

ok("tda_rng: different seeds give different streams",
   !identical(tda_rng(1)$rd(5), tda_rng(2)$rd(5)))
ok("tda_rng: same seed gives the same stream, deterministically",
   identical(tda_rng(99)$rd(10), tda_rng(99)$rd(10)))

## plot.tda_graph(): TDA's pltree refuses a directed graph
## ("graph must be undirected"), which is a real constraint of pltree
## itself, not something this package adds -- but relaying that raw
## message verbatim gives no hint that the fix is tda_graph(...,
## directed = FALSE), found only by asking why a parameter that only
## errors is worth having at all. A directed tree-shaped graph (built
## with the default directed = TRUE, which is right for an ordinary
## graph but not a tree) now gets a specific, actionable message
## instead of TDA's bare one.
ge_tree <- data.frame(from = c(1, 1, 2), to = c(2, 3, 4))
g_directed_tree <- tda_graph(ge_tree, directed = TRUE)
err_tree <- tryCatch({
    pdf(NULL); on.exit(dev.off())
    plot(g_directed_tree, layout = "tree")
    NULL
}, error = function(e) conditionMessage(e))
ok("plot.tda_graph: a directed graph forced into layout=\"tree\" gets a specific, actionable error",
   !is.null(err_tree) && grepl("directed = FALSE", err_tree, fixed = TRUE))
g_undirected_tree <- tda_graph(ge_tree, directed = FALSE)
ok("plot.tda_graph: the same tree with directed = FALSE draws without error",
   isTRUE(tryCatch({
       pdf(NULL); on.exit(dev.off())
       plot(g_undirected_tree, layout = "tree")
       TRUE
   }, error = function(e) FALSE)))


## the argument grammar ---------------------------------------------------
#
# ?tdaR promises that a handful of argument names mean the same thing
# wherever they appear.  Enforce the vocabulary: the synonyms that used to
# exist (cen, w, group-for-sel=) must not come back in exported signatures.
grammar_fns <- Filter(function(n) is.function(getExportedValue("tdaR", n)),
                      getNamespaceExports("tdaR"))
sig <- lapply(grammar_fns, function(n) names(formals(getExportedValue("tdaR", n))))
names(sig) <- grammar_fns
uses <- function(a) grammar_fns[vapply(sig, function(s) a %in% s, NA)]
assert("argument grammar: no 'cen' where the vocabulary is 'censor'",
       length(uses("cen")) == 0L)
assert("argument grammar: 'w' only where the weights ARE the object",
       all(uses("w") %in% c("tda_cwt_norm")))
assert("argument grammar: sel= is reachable as 'select' on the pl side too",
       "select" %in% sig$tda_pl_lines, !"group" %in% sig$tda_pl_lines)
assert("argument grammar: fits lead with formula, data",
       all(vapply(c("tda_rate", "tda_lsreg", "tda_glm", "tda_qreg",
                    "tda_ple", "tda_ltb", "tda_loglin"),
                  function(n) identical(sig[[n]][1:2], c("formula", "data")),
                  NA)))
assert("argument grammar: every options= is a list() escape hatch",
       all(vapply(uses("options"), function(n)
           identical(formals(getExportedValue("tdaR", n))$options,
                     quote(list())), NA)))
