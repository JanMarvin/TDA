# ---- zoo archives -----------------------------------------------------

#' Read and write Zoo archives
#'
#' \code{unzoo} extracts a \code{.zoo} archive, and \code{zoo} creates one --
#' \code{\link[utils]{unzip}} and \code{\link[utils]{zip}} for Zoo rather
#' than Zip. This is tdaR's code, not part of TDA (TDA can only load an
#' archive's contents into its data matrix, via a command that needs a
#' companion \code{.zad} description file, and has no general extract or
#' create capability at all).
#'
#' \code{unzoo} decodes both of Zoo's packing methods, LZD and LZH, ported
#' from TDA's reference implementation -- the same
#' code path TDA itself has used since the 1990s to read \code{.zoo}
#' archives, so archives TDA itself could read, this reads too. Both
#' methods are tested against a historical TDA archive (\code{deha1.zoo}
#' from the TDA teaching pages, kept outside the package as an external
#' test fixture): every decoded byte matches the CRC-16 the original
#' archive itself stored, and separately, every file
#' \code{zoo} writes here was extracted correctly, with a passing CRC-16,
#' by the independent, unrelated \code{zoo}/\code{unzoo} reference tools
#' (\url{https://github.com/troglobit/zoo}) -- not just read back
#' correctly by this same code, which would only prove self-consistency.
#' \code{zoo} writes LZD-compressed entries by default, or stored
#' (uncompressed) ones with \code{method = 0}, each with its CRC-16
#' checksum, readable by any Zoo implementation including TDA's.
#'
#' Member names are Zoo 2.1's: a long file name and a directory name are
#' stored in the directory entry beside the DOS-style 8.3 short name, and
#' both \code{unzoo} and \code{zoo} use them, so \code{sub/long_name.dat}
#' round-trips as itself. \code{unzoo} creates the directories it needs
#' under \code{exdir}. An archive \code{zoo} writes is, field for field,
#' what the reference \code{zoo} tool writes for the same members (its
#' header, every directory entry's names and version bytes, the file
#' leaders, the terminal record; only the members' dates differ, which
#' \code{zoo} leaves unset), and TDA's own \code{arcd} finds a member by
#' its full stored name. TDA's archive description file separates its
#' fields by blanks, so a member TDA is to read may not have a blank in
#' its name.
#'
#' Zoo can, in principle, keep more than one stored version of a file under
#' the same name (its version-history feature) -- \code{unzoo} does not
#' collapse or deduplicate by name, so if an archive has several entries
#' called the same thing, all of them come back as separate rows/files
#' rather than only the latest.
#'
#' @param zoofile path to a \code{.zoo} file, or an \code{http(s)://} URL,
#'   which is downloaded to a temporary file first.
#' @param files which members to extract; \code{NULL} (the default) means
#'   all of them. Members inside a directory are named with it,
#'   \code{"sub/x.dat"}.
#' @param method packing method for the members: \code{1} (LZD, the
#'   default) or \code{0} (stored, uncompressed). TDA's \code{arcd}
#'   reads all of 0, 1 and 2, but other readers of TDA-written archives
#'   may not, and LZD is what \code{zoo} itself has always produced.
#' @param list if \code{TRUE}, return a data frame describing the archive's
#'   contents instead of extracting anything.
#' @param exdir directory to extract into, created if it does not exist.
#' @param overwrite whether to replace files already present in \code{exdir}.
#' @return For \code{list = TRUE}, a data frame with one row per member:
#'   \code{name} (with its directory, if the entry stores one), \code{method} (0 stored, 1 LZD, 2 LZH), \code{size}
#'   (original, uncompressed), and \code{ok} (whether this reads correctly --
#'   always \code{TRUE} for methods 0/1/2, since those are the only methods
#'   Zoo defines). Otherwise, the paths of the files written, invisibly the
#'   same as \code{\link[utils]{unzip}}.
#' @family zoo archives
#' @examples
#' archive <- system.file("extdata", "tda.zoo", package = "tdaR")
#' unzoo(archive, list = TRUE)
#'
#' ex <- tempfile()
#' paths <- unzoo(archive, files = "avar.dat", exdir = ex)
#' readLines(paths)
#'
#' # round-trip: write what was just extracted back out, and read it again
#' out <- tempfile(fileext = ".zoo")
#' zoo(out, paths)
#' unzoo(out, list = TRUE)
#'
#' # names longer than 8.3 and a directory are kept as given
#' zoo(out, paths, names = "results/avar_first_run.dat")
#' unzoo(out, list = TRUE)
#' ex2 <- tempfile()
#' unzoo(out, exdir = ex2)
#' list.files(ex2, recursive = TRUE)
#' \dontrun{
#' unzoo("https://example.com/archive.zoo", exdir = tempfile())
#' }
#' @export
unzoo <- function(zoofile, files = NULL, list = FALSE, exdir = ".",
                  overwrite = TRUE) {
    if (grepl("^https?://", zoofile)) {
        tmp <- tempfile(fileext = ".zoo")
        utils::download.file(zoofile, tmp, mode = "wb", quiet = TRUE)
        on.exit(unlink(tmp), add = TRUE)
        zoofile <- tmp
    }
    zoofile <- normalizePath(zoofile, mustWork = TRUE)
    r <- .Call(C_tda_unzoo, zoofile)
    d <- data.frame(name = r$name, method = r$method, size = r$size,
                    ok = r$ok, stringsAsFactors = FALSE)
    if (list)
        return(d)

    keep <- if (is.null(files)) rep(TRUE, nrow(d)) else d$name %in% files
    if (any(!d$ok[keep]))
        warning("could not decode: ",
                paste(d$name[keep & !d$ok], collapse = ", "))
    if (!dir.exists(exdir))
        dir.create(exdir, recursive = TRUE)

    bad <- grepl("^/|^[.][.]/|/[.][.]/|/[.][.]$|^[.][.]$", d$name[keep])
    if (any(bad))
        stop("member name(s) leave the extraction directory: ",
             paste(d$name[keep][bad], collapse = ", "), call. = FALSE)
    paths <- character(0)
    for (i in which(keep & d$ok)) {
        p <- file.path(exdir, d$name[i])
        if (!overwrite && file.exists(p))
            next
        if (!dir.exists(dirname(p)))
            dir.create(dirname(p), recursive = TRUE)
        writeBin(r$data[[i]], p)
        paths <- c(paths, p)
    }
    invisible(paths)
}

#' @param files paths of the files to add.
#' @param names the names the members get in the archive, one per file:
#'   by default the files' base names. A name may carry a directory,
#'   \code{"sub/x.dat"}, with \code{/} as the separator; it must be
#'   relative and may not contain \code{..}.
#' @rdname unzoo
#' @export
zoo <- function(zoofile, files, method = 1L, names = basename(files)) {
    # method 1 (LZD) by default: TDA's arcd refuses any packing method
    # outside 1..2, so a stored archive is the one archive TDA cannot open.
    # method = 0 stores, for a reader that wants no decompression at all.
    if (!method %in% c(0L, 1L))
        stop("`method` must be 0 (stored) or 1 (LZD)", call. = FALSE)
    missing_f <- files[!file.exists(files)]
    if (length(missing_f))
        stop("file(s) not found: ", paste(missing_f, collapse = ", "))
    if (length(names) != length(files) || anyNA(names))
        stop("`names` must name every file", call. = FALSE)
    names <- gsub("\\\\", "/", names)
    names <- sub("^[.]/", "", names)
    bad <- !nzchar(basename(names)) |
        grepl("^/|^[.][.]/|/[.][.]/|/[.][.]$|^[.][.]$|//", names)
    if (any(bad))
        stop("member name(s) must be relative paths without `..`: ",
             paste(names[bad], collapse = ", "), call. = FALSE)
    if (any(nchar(basename(names), "bytes") > 255L |
            nchar(dirname(names), "bytes") > 255L))
        stop("member name or directory longer than 255 bytes", call. = FALSE)
    contents <- lapply(files, function(f) readBin(f, "raw", n = file.size(f)))
    .Call(C_tda_zoo, path.expand(zoofile), names, contents,
          as.integer(method))
    invisible(zoofile)
}
