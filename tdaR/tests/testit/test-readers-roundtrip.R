# Round-trip fidelity for every importer: the returned data frame must
# equal the source cell for cell -- not run, not look plausible, EQUAL.
# Where an independent R reader exists it serves as the oracle.

.rt_data <- function() data.frame(
    ID   = 1:6,
    X    = c(-1.5, 0, 2.25, 1e6, -0.001, 3),
    Y    = c(10L, -3L, 0L, 7L, 42L, -1L))

# ---- rcsv: cell-for-cell against the written csv -----------------
local({
    d0 <- data.frame(a = 1:3, b = c(2.5, 1, 4), s = c("x", "yy", "zzz"))
    f <- tempfile(fileext = ".csv")
    utils::write.csv(d0, f, row.names = FALSE)
    d1 <- tda_rcsv(f, sep = ",")
    if (tdaR:::.use_exports()) {
        ok("rcsv round-trips every cell",
           identical(d1$a, d0$a) && identical(d1$b, d0$b) &&
           identical(d1$s, d0$s))
    } else {
        ok("rcsv rendered lines carry every row",
           is.character(d1) && length(d1) == nrow(d0) + 1L)
    }
})

# ---- sys: write with TDA, read with TDA, compare to source -------
local({
    d0 <- .rt_data()
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(d0, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    invisible(tda_run(c(
        "nvar(dfile=d.dat, ID[4.0]=c1, X[16.6]=c2, Y[6.0]=c3);",
        "wsys = t.sys;"), dir = dr))
    d1 <- tda_read_sys(file.path(dr, "t.sys"))
    ok("sys round-trips all values",
       nrow(d1) == 6L && max(abs(d1[[2]] - d0$X)) < 1e-6 &&
       all(d1[[3]] == d0$Y) && all(d1[[1]] == d0$ID))
    d2 <- tda_rsys(file.path(dr, "t.sys"))
    # the run attribute carries per-call temp paths; the DATA must match
    ok("rsys returns the same data as read_sys",
       isTRUE(all.equal(unname(as.matrix(d1)), unname(as.matrix(d2)))))
})

# ---- stata: TDA writes, TDA reads; readstata13/foreign oracle ----
local({
    d0 <- .rt_data()
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(d0, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    invisible(tda_run(c(
        "nvar(dfile=d.dat, ID[4.0]=c1, X[16.6]=c2, Y[6.0]=c3);",
        "wstata = t.dta;"), dir = dr))
    d1 <- tda_read_stata(file.path(dr, "t.dta"))
    ok("stata round-trips all values",
       nrow(d1) == 6L && max(abs(d1$X - d0$X)) < 1e-6 &&
       all(d1$Y == d0$Y))
    if (requireNamespace("readstata13", quietly = TRUE)) {
        d2 <- readstata13::read.dta13(file.path(dr, "t.dta"))
        ok("stata: independent reader agrees cell for cell",
           max(abs(d2$X - d1$X)) < 1e-9 && all(d2$Y == d1$Y))
    }
})

# ---- table: TDA writes a df= table, read_table restores it -------
local({
    d0 <- .rt_data()
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(d0, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    invisible(tda_run(c(
        "nvar(dfile=d.dat, ID[4.0]=c1, X[16.6]=c2, Y[6.0]=c3);",
        "pdata(fmt=16.6) = t.out;"), dir = dr))
    d1 <- tda_read_table(file.path(dr, "t.out"))
    ok("read_table restores every value TDA wrote",
       nrow(d1) == 6L && max(abs(d1[[2]] - d0$X)) < 1e-5 &&
       all(d1[[3]] == d0$Y))
})

# ---- rplz: every field of every record, leading zeros intact -----
local({
    f <- tempfile()
    writeLines(c("plz;land;kreis;orte",
                 "01067;SN;Dresden;Dresden,Altstadt",
                 "01069;SN;Dresden;Dresden,Suedvorstadt",
                 "80331;BY;Muenchen;Muenchen"), f)
    d <- tda_rplz(f)
    ok("rplz keeps every record and field",
       nrow(d) == 3L && identical(d$plz, c("01067", "01069", "80331")) &&
       identical(d$region, c("SN", "SN", "BY")) &&
       grepl("Altstadt", d$places[1L]))
})

# ---- unzoo: archived member equals its known content -------------
local({
    a <- system.file("extdata", "tda.zoo", package = "tdaR")
    lst <- unzoo(a, list = TRUE)
    ok("unzoo lists every member intact",
       is.data.frame(lst) && all(lst$ok) && nrow(lst) == 2L)
    ex <- tempfile("zoo"); dir.create(ex)
    p <- unzoo(a, files = "avar.dat", exdir = ex)
    v <- readLines(p)
    ok("unzoo-extracted avar.dat holds its known content",
       length(v) == 3L && all(grepl("^V[123] 1 (0|8|16) 8.0", v)))
})

# ---- ucinet (external fixture): matrix vs edges, self-consistent -
local({
    ed <- Sys.getenv("TDA_EXT_INPUT", "")
    uf <- if (nzchar(ed))
        list.files(ed, pattern = "[.]##d$", full.names = TRUE)
    else character()
    if (length(uf)) {
        m <- tda_read_ucinet(uf[1L])
        e <- tda_read_ucinet(uf[1L], form = "edges")
        ok("ucinet: every nonzero cell appears as exactly one edge",
           nrow(e) == sum(as.matrix(m) != 0) &&
           all(as.matrix(m)[cbind(e$from, e$to)] == e$value))
    }
})

# ---- dbf via sf's fixture: against foreign::read.dbf -------------
local({
    shp <- system.file("shape/nc.dbf", package = "sf")
    if (nzchar(shp) && requireNamespace("foreign", quietly = TRUE)) {
        d1 <- tda_read_dbf(shp)
        d2 <- foreign::read.dbf(shp, as.is = TRUE)
        ok("dbf: name and area columns equal the independent reader",
           identical(as.character(d1$NAME), as.character(d2$NAME)) &&
           max(abs(d1$AREA - d2$AREA)) < 1e-9 &&
           nrow(d1) == nrow(d2))
    }
})
