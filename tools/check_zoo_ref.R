## Cross-check zoo()/unzoo() against the reference zoo tool
## (https://github.com/troglobit/zoo), in both directions:
##
##   1. an archive zoo() wrote is listed and extracted by the reference
##      tool, which verifies each member's CRC-16 and its directory
##      entry's CRC as it goes; every extracted file must be byte-identical
##      to the input;
##   2. an archive the reference tool wrote ("zoo a", which stores
##      directories) is read by unzoo() with the same names and bytes;
##   3. the archive header, every directory entry's name fields (the 8.3
##      short name, the long name, the directory) and version bytes, the
##      file leader before each member and the terminal record are
##      identical in the two archives.
##
## Usage: Rscript tools/check_zoo_ref.R /path/to/zoo

args <- commandArgs(TRUE)
if (length(args) != 1L || !file.exists(args[1L]))
    stop("usage: Rscript tools/check_zoo_ref.R /path/to/zoo")
ZOO <- normalizePath(args[1L])
library(tdaR)

wd <- tempfile("zooref"); dir.create(wd)
old <- setwd(wd); on.exit(setwd(old))
dir.create("in/results/sub", recursive = TRUE)
dir.create("in/results/deeper/still", recursive = TRUE)
nm <- c("short.dat", "a_much_longer_file_name_than_dos_allows.txt",
        "results/sub/nested-name.with.dots.csv", "README", ".hidden",
        "UPPER_CASE_NAME.DATA", "x", "trailing.dot.", "eight_ch",
        "abcdefghijkl", "data.json", "a.b.c.d", "my file (1)+~.txt",
        "results/deeper/still/x.y")
for (x in nm)
    writeLines(rep(paste("line of", x), 40), file.path("in", x))

## 1. ours, read by the reference
zoo("ours.zoo", file.path("in", nm), names = nm)
dir.create("ref")
setwd("ref")
st <- system2(ZOO, c("x//", "../ours.zoo"), stdout = TRUE, stderr = TRUE)
setwd(wd)
if (!all(grepl("-- extracted", grep("^Zoo:", st, value = TRUE))))
    stop("reference zoo did not extract every member cleanly:\n",
         paste(st, collapse = "\n"))
for (x in nm)
    if (!identical(readBin(file.path("in", x), "raw", 1e6),
                   readBin(file.path("ref", x), "raw", 1e6)))
        stop("reference extraction of ", x, " differs from the input")

## 2. theirs, read by us
setwd("in")
system2(ZOO, c("a", "../theirs.zoo", shQuote(nm)), stdout = FALSE)
setwd(wd)
lst <- unzoo("theirs.zoo", list = TRUE)
if (!setequal(lst$name, nm) || !all(lst$ok))
    stop("unzoo() does not list the reference archive's members as written")
unzoo("theirs.zoo", exdir = "back")
for (x in nm)
    if (!identical(readBin(file.path("in", x), "raw", 1e6),
                   readBin(file.path("back", x), "raw", 1e6)))
        stop("unzoo() extraction of ", x, " differs from the input")

## 3. directory-entry name fields, both archives
entries <- function(path) {
    d <- readBin(path, "raw", file.size(path))
    u16 <- function(i) as.integer(d[i + 1L]) + 256L * as.integer(d[i + 2L])
    u32 <- function(i) u16(i) + 65536 * u16(i + 2L)
    off <- u32(24L)
    out <- list()
    header <- d[1:42]
    repeat {
        nxt <- u32(off + 6L)
        if (nxt == 0) {
            terminal <- d[off + seq_len(length(d) - off)]
            break
        }
        ofs <- u32(off + 10L)
        leader <- rawToChar(d[ofs - 4:0][1:4])
        fn <- rawToChar(d[off + 39L:51L][seq_len(
            which(d[off + 39L:51L] == as.raw(0))[1L] - 1L)])
        namlen <- as.integer(d[off + 57L]); dirlen <- as.integer(d[off + 58L])
        lf <- if (namlen > 1L) rawToChar(d[off + 59L + seq_len(namlen - 1L) - 1L]) else ""
        dn <- if (dirlen > 1L) rawToChar(d[off + 59L + namlen + seq_len(dirlen - 1L) - 1L]) else ""
        out[[length(out) + 1L]] <- data.frame(type = as.integer(d[off + 5L]),
                                              fname = fn, lfname = lf, dir = dn,
                                              vmaj = as.integer(d[off + 29L]),
                                              vmin = as.integer(d[off + 30L]),
                                              leader = leader)
        off <- nxt
    }
    o <- do.call(rbind, out)
    list(header = header, entries = o[order(o$fname, o$lfname, o$dir), ],
         terminal = terminal)
}
a <- entries("ours.zoo"); b <- entries("theirs.zoo")
rownames(a$entries) <- rownames(b$entries) <- NULL
if (!identical(a$entries, b$entries)) {
    print(a$entries); print(b$entries)
    stop("directory-entry fields differ from the reference tool's")
}
if (!identical(a$header, b$header))
    stop("archive header differs from the reference tool's: ",
         paste(a$header, collapse = ""), " vs ", paste(b$header, collapse = ""))
if (!identical(a$terminal, b$terminal))
    stop("terminal directory record differs from the reference tool's")
cat("zoo()/unzoo() agree with the reference zoo on", length(nm), "members\n")
