# zoo()/unzoo() -- these need no TDA run at all, so unlike the other test
# files here they do not depend on TDA_EXAMPLES / have_examples.

## The LZH decoder is checked against deha1.zoo, the archive from the TDA
## teaching pages (every compressed member uses LZH, packing method 2),
## kept outside the package as an external fixture -- 67 files, not the
## phantom 68th a stray terminal directory record used to add.
if (have_ext("deha1.zoo")) {
archive <- file.path(EXT, "deha1.zoo")
lst <- unzoo(archive, list = TRUE)
ok("unzoo: lists the real member count, not one more", nrow(lst) == 67L)
ok("unzoo: every member decodes correctly", all(lst$ok))
ok("unzoo: methods are only ones Zoo defines", all(lst$method %in% 0:2))
ok("unzoo: sizes are all recorded and positive",
   all(is.finite(lst$size) & lst$size > 0))

## Extracting one member and checking its content against what the archive
## is known to contain.
ex <- tempfile("unzoo")
p <- unzoo(archive, files = "eha1.dat", exdir = ex)
ok("unzoo: returns the path it wrote, like unzip()", file.exists(p))
eha1 <- read.table(p, col.names = c("ID", "DUR", "CEN"))
same("unzoo: LZH-decoded content is exactly right", eha1$DUR,
     c(17, 5, 22, 13, 2, 9, 12, 15))
same("unzoo: including the censoring indicator", eha1$CEN,
     c(1, 0, 1, 1, 0, 1, 0, 1))

## The command file for the same example, to check LZH text content too
## (eha1.dat above is numeric; this exercises a different byte pattern).
cf <- unzoo(archive, files = "eha1.cf", exdir = ex)
ok("unzoo: LZH-decoded text is intact",
   any(grepl("dfile = eha1.dat", readLines(cf), fixed = TRUE)))

## The one stored (uncompressed, method 0) member in this archive.
ok("unzoo: the archive's one stored member decodes too",
   lst$ok[lst$name == "eha5.dat" & lst$method == 0L])
}

## Round trip on the package's own LZD archive: extract, write back out,
## and read it back. zoo() writes CRC-16 checksums; a bad one would still
## list, but this checks the reader agrees, not just that nothing errored.
archive <- system.file("extdata", "tda.zoo", package = "tdaR")
ex <- tempfile("unzoo")
p <- unzoo(archive, files = "avar.dat", exdir = ex)
allpaths <- unzoo(archive, exdir = ex)
rt <- tempfile(fileext = ".zoo")
zoo(rt, allpaths)
rtlst <- unzoo(rt, list = TRUE)
ok("zoo: round trip keeps every file", nrow(rtlst) == length(allpaths))
ok("zoo: round trip decodes cleanly", all(rtlst$ok))
ok("zoo: round-tripped names match", setequal(rtlst$name, basename(allpaths)))

ex2 <- tempfile("unzoo2")
back <- unzoo(rt, files = "avar.dat", exdir = ex2)
same("zoo: round-tripped content is byte-identical",
     readBin(back, "raw", file.size(back)),
     readBin(p, "raw", file.size(p)))

## Long names and directories (Zoo 2.1's type-2 directory entries).
## long.zoo was written by the reference zoo tool (github.com/troglobit/zoo,
## "zoo a") from TDA's own example data under long names; the fields of
## every entry zoo() writes were compared byte for byte with that tool's
## by tools/check_zoo_ref.R.
ref <- system.file("extdata", "long.zoo", package = "tdaR")
rl <- unzoo(ref, list = TRUE)
same("unzoo: long names and directories from the reference tool's archive",
     sort(rl$name),
     c("data/avar_description.dat", "data/sub/adata_first_wave_2026.dat"))
ok("unzoo: the reference tool's LZD members decode", all(rl$ok))
ex3 <- tempfile("unzoo3")
got <- unzoo(ref, exdir = ex3)
ok("unzoo: creates the member's directory under exdir",
   file.exists(file.path(ex3, "data", "sub", "adata_first_wave_2026.dat")))
same("unzoo: the long-named member is the example data",
     readBin(file.path(ex3, "data", "sub", "adata_first_wave_2026.dat"),
             "raw", 480L),
     readBin(file.path(ex, "adata.dat"), "raw", 480L))

## TDA itself finds the members by their full stored names.
dr <- tempfile("tda"); dir.create(dr)
file.copy(ref, file.path(dr, "long.zoo"))
file.copy(system.file("extdata", "long.zad", package = "tdaR"),
          file.path(dr, "long.zad"))
r <- tda_run(c("arcd = long.zad;", "arcc;",
               "nvar(V1<5>[8.0] = A:V1, V2<5>[8.0] = A:V2, V3<5>[8.0] = A:V3);",
               "dstat;"), dir = dr)
ok("arcd: TDA finds a member by its long name inside a directory",
   any(grepl("Using archive data file: data/sub/adata_first_wave_2026.dat",
             r$output, fixed = TRUE)) && r$errors == 0)
same("arcd: the data read through the long name are the example data",
     r$exports$dstat.stats[, 3L], c(10.5, 12.3, 8.9))

## zoo() writes the same: a long name, a directory, a name with several
## dots, a dot-file, and an 8.3 name that needs no long form.
src <- tempfile("zoosrc"); dir.create(file.path(src, "results", "sub"),
                                      recursive = TRUE)
nm <- c("short.dat", "a_much_longer_file_name_than_dos_allows.txt",
        "results/sub/nested-name.with.dots.csv", ".hidden")
for (x in nm) writeLines(rep(paste("line of", x), 50), file.path(src, x))
lz <- tempfile(fileext = ".zoo")
zoo(lz, file.path(src, nm), names = nm)
ll <- unzoo(lz, list = TRUE)
same("zoo: long names and directories are stored as given", ll$name, nm)
ex4 <- tempfile("unzoo4")
unzoo(lz, exdir = ex4)
for (x in nm)
    same(paste("zoo: round trip of", x),
         readLines(file.path(ex4, x)), readLines(file.path(src, x)))
## ... and TDA's arcd reads what zoo() wrote, by the long names.
dr2 <- tempfile("tda"); dir.create(dr2)
zoo(file.path(dr2, "mine.zoo"),
    c(file.path(ex, "adata.dat"), file.path(ex, "avar.dat")),
    names = c("wave-one/adata_first_wave_2026.dat", "avar_description.dat"))
writeLines(c("mine.zoo", "1 wave-one/adata_first_wave_2026.dat 1 24 20 3",
             "2 avar_description.dat 2 40 3 0"), file.path(dr2, "mine.zad"))
r2 <- tda_run(c("arcd = mine.zad;",
                "nvar(V1<5>[8.0] = A:V1, V2<5>[8.0] = A:V2, V3<5>[8.0] = A:V3);",
                "dstat;"), dir = dr2)
same("arcd: TDA reads an archive zoo() wrote, through its long names",
     r2$exports$dstat.stats[, 3L], c(10.5, 12.3, 8.9))

ok("zoo: a member name may not climb out of the archive",
   inherits(tryCatch(zoo(lz, file.path(src, nm[1L]), names = "../x"),
                     error = function(e) e), "error"))
ok("zoo: nor be absolute",
   inherits(tryCatch(zoo(lz, file.path(src, nm[1L]), names = "/x"),
                     error = function(e) e), "error"))

unlink(c(ex, ex2, ex3, ex4, rt, lz, src, dr, dr2), recursive = TRUE)
