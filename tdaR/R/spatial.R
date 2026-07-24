# Wrapping TDA's spatial commands.
#
# These work on a spatial data structure rather than the ordinary data matrix:
# sdgen writes one from a data frame of coordinates, sdnvar reads one back in,
# and the analysis and plotting commands then operate on it.  So the wrapper
# is an object holding that structure, built once and passed to the commands,
# in the same way tda_graph works for the graph family.
#
# sf does all of this better.  These exist because this is meant to be TDA in
# R rather than a selection from it.

#' Spatial data for TDA's spatial commands
#'
#' Builds TDA's spatial data structure from a data frame of coordinates. Each
#' object -- a point, a line or a polygon -- is a group of rows sharing an id.
#'
#' @param data a data frame with an id column and two coordinate columns.
#' @param id,x,y the column names.
#' @param attributes further columns to carry along, TDA's \code{av=}.
#' @param type type of spatial object, TDA's \code{typ=}: 1 for points.
#' @param options a named list of further \code{sdgen} options.
#' @param dir working directory; the spatial file is written there and stays
#'   for the commands that read it.
#' @return A \code{tda_spatial}.
#' @family spatial analysis
#' @examples
#' d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
#' s <- tda_spatial(d)
#' tda_sd_info(s)
#' @export
tda_spatial <- function(data, id = "id", x = "x", y = "y",
                        attributes = NULL, type = 1, options = list(),
                        dir = tempfile("tda")) {
    d <- as.data.frame(data)
    for (v in c(id, x, y))
        if (!v %in% names(d))
            stop("no column called '", v, "' in `data`")
    keep <- c(id, x, y, attributes)
    d <- d[, keep, drop = FALSE]
    lab <- names(d)
    names(d) <- .tda_names(lab)

    opts <- c(list(id = names(d)[1L],
                   xyv = paste(names(d)[2:3], collapse = ","),
                   typ = type), .tda_extra(options))
    if (!is.null(attributes))
        opts$av <- paste(names(d)[-(1:3)], collapse = ",")

    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "sdgen"), opts,
                                          list(rhs = "spatial.sd")))),
                   data = d, dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not build the spatial data: ", err[1L], call. = FALSE)
    # the copy the caller sees keeps the caller's column names; the
    # .tda_names() spellings existed only for TDA's benefit
    names(d) <- lab
    structure(list(dir = dir, file = "spatial.sd", data = d, xlab = lab,
                   n = nrow(d), run = res), class = "tda_spatial")
}

# sdnvar needs the variable definitions as well as the file: given only dfile
# it reports "No variables defined".  check_sd() in t_gdat.c names the four
# it insists on -- SDID, SDTyp, SDN and SDPtr, by those exact names, with
# SDPtr required to be type 3 where the others may be 2 or 3.
#
# The declaration below is the one that works, and each piece is load-
# bearing: SDID/SDTyp/SDN are columns 1-3 of the object description record
# that sdgen writes ("id type count", then count coordinate pairs); SDPtr is
# the dummy TDA's manual calls it and consumes NO column -- it is declared
# as `rd`, a uniform random deviate (not "read", which is what it looks
# like), with an explicit <5> width to satisfy the type check.  Declaring
# SDPtr = 3 instead fails check_sd(); declaring all four as rd fills the
# data matrix with random numbers and sdinf rejects it ("SDN contains a
# value less than 1" -- the error an earlier, wrong version of this block
# hit, preserved here because it is the symptom you see whenever the
# column mapping drifts).  With this shape the whole reading chain works:
# sdinf reports the right object counts and bounds, sdvd computes Voronoi
# diagrams and Delaunay triangulations, and polygon-type files (typ=3)
# read back with their segment statistics -- pinned in
# tests/testit/test-print.R and test-models.R.
# The object description record carries any attributes after the object
# fields, in the order the file was written, and where they start depends on
# who wrote it: sdgen puts them straight after SDID, SDTyp and SDN, so at
# c4, while sdshp writes the shape number and the part number first, so at
# c6.  SDPtr is `rd`, a random deviate, and consumes no column at all.  The
# structure carries its offset rather than the caller counting fields.
# The x/y range to use when none was given: from $data directly for a
# tda_spatial() structure built from an R data frame, or parsed from
# tda_sd_info()'s reported bounds for anything else (a shapefile or
# other file read straight by TDA, with no R data frame ever involved).
# Used by plot.tda_spatial(), tda_sd_ps3(), and the legacy .sd_plot3() --
# previously duplicated in the first two, and used to default to c(0, 1)
# in the third, drawing everything outside the visible window for any
# real-world coordinate data.
.sd_bounds <- function(s) {
    if (!is.null(s$data))
        return(list(xlim = range(s$data[[2L]], na.rm = TRUE),
                    ylim = range(s$data[[3L]], na.rm = TRUE)))
    info <- tda_sd_info(s)$text
    xr <- as.numeric(regmatches(info, regexpr("(?<=XMin:)\\s*-?[0-9.]+",
                                               info, perl = TRUE)))
    xr <- c(xr, as.numeric(regmatches(info,
                regexpr("(?<=XMax:)\\s*-?[0-9.]+", info, perl = TRUE))))
    yr <- as.numeric(regmatches(info, regexpr("(?<=YMin:)\\s*-?[0-9.]+",
                                               info, perl = TRUE)))
    yr <- c(yr, as.numeric(regmatches(info,
                regexpr("(?<=YMax:)\\s*-?[0-9.]+", info, perl = TRUE))))
    if (!length(xr) || !length(yr))
        stop("no coordinate range available: pass xlim/ylim explicitly")
    list(xlim = range(xr, na.rm = TRUE), ylim = range(yr, na.rm = TRUE))
}

.sd_open <- function(s, attributes = NULL) {
    # TDA variable names must start with an uppercase letter (see
    # .tda_names() and the "Variable names" section on ?tda_write_data) --
    # attribute names come straight from the R data frame's column
    # names here, which are not guaranteed to start uppercase, and were
    # never run through .tda_names() before. TDA's error for this
    # ("Syntax error: x=c4,)") does not mention case at all, so this is
    # easy to rediscover the hard way; see the "Variable names" section.
    att <- .tda_names(attributes)
    sprintf(paste("sdnvar(",
                  "    dfile = %s,",
                  "    SDID     = c1,",
                  "    SDTyp    = c2,",
                  "    SDN      = c3,",
                  "    SDPtr<5> = rd,",
                  if (length(att))
                      paste(sprintf("    %s = c%d,", att,
                                    (s$attr_col %||% 4L) - 1L +
                                        seq_along(att)),
                            collapse = "\n"),
                  ");", sep = "\n"), s$file)
}

#' Run a spatial command
#'
#' \code{tda_sd} reaches any of TDA's spatial commands; the named functions
#' are the ones with an obvious meaning. The spatial data is opened with
#' \code{sdnvar} first, which is what the commands operate on.
#'
#' @param s a \code{\link{tda_spatial}}.
#' @param cmd the command name, e.g. \code{"sdinf"}, \code{"sdvd"}.
#' @param options a named list of options for the command.
#' @param rhs the command's right-hand side, when it takes one.
#' @param setup extra commands to run before it, such as a PostScript setup.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying the run and, where the command writes one, a
#'   \code{table}.
#' @family spatial analysis
#' @examples
#' d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
#' s <- tda_spatial(d)
#' r <- tda_sd(s, "sdinf")
#' cat(r$run$output, sep = "\n")
#' @export
tda_sd <- function(s, cmd, options = list(), rhs = NULL, setup = character(),
                   ...) {
    stopifnot(inherits(s, "tda_spatial"))
    res <- tda_run(c(.sd_open(s), setup,
                     do.call(tda_block, c(list(name = cmd), .tda_extra(options),
                                          list(rhs = rhs)))),
                   dir = s$dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ", cmd, ": ", err[1L], call. = FALSE)
    out <- if (!is.null(rhs) && grepl("\\.txt$", rhs)) {
        .read_g_output(res, rhs)
    } else {
        list(table = NULL, text = NULL)
    }
    # A command with its own producer overlays the parsed frame.  The
    # spatial files are written at fmt=, default 10.4, which is fine for
    # a unit square and destroys real geographic coordinates: reading
    # zimbabwe.pnt back gave 30.2718 for a longitude of
    # 30.271785736084.  SDVarX/Y are double, so the export has the value
    # the file rounded away.
    if (!is.null(out$table))
        out$table <- .overlay_num(out$table,
                                  res$exports[[paste0(cmd, ".table")]])
    if (is.null(out$table) && is.null(out$text)) {
        # Either no file was requested (sdinf), or one was and TDA wrote
        # nothing real to it ("0 records written") because its actual
        # answer -- sdrel's adjacency comparisons, say -- went to the
        # console instead. Either way, the block of text between the last
        # two "----" rules brackets one command's printed output in
        # TDA's log, and is the only place left to find it. Previously
        # this fallback only ran when no file was requested at all, so a
        # command that asked for one and got nothing there still lost its
        # real answer -- $table and $text both NULL despite TDA having
        # printed it in plain sight.
        out$text <- .last_output_block(res$output)
    }
    # sdinf answers in prose; the counts and the bounding box are its
    # whole content, so they come back as a one-row table, numbers from
    # the console tap
    if (cmd == "sdinf" && is.null(out$table))
        out$table <- .sdinf_table(res)
    structure(list(call = sys.call(-1L), run = res, command = cmd, n = s$n,
                   table = out$table, text = out$text, spatial = s, rhs = rhs),
              class = c(paste0("tda_", cmd), "tda_table"))
}

.sdinf_table <- function(res) {
    out <- res$output
    one <- function(pat, which = 1L) {
        i <- grep(pat, out)
        if (!length(i)) return(NA_real_)
        v <- .tap_line(res, i[1L])
        if (!is.null(v) && length(v) >= which) return(v[which])
        m <- regmatches(out[i[1L]], gregexpr("-?[0-9]+\\.?[0-9]*(e[+-]?[0-9]+)?", out[i[1L]]))[[1L]]
        if (length(m) >= which) as.numeric(m[which]) else NA_real_
    }
    d <- data.frame(points = one("^Number of points:"),
                    lines = one("^Number of lines:"),
                    polygons = one("^Number of polygons:"),
                    xmin = one("^XMin:", 1L), ymin = one("^XMin:", 2L),
                    xmax = one("^XMax:", 1L), ymax = one("^XMax:", 2L))
    if (all(is.na(d))) NULL else d
}

# The text strictly between the last two "----...----" rule lines in a run's
# output -- which brackets the most recently run command's printed
# answer, for commands that write nothing to a file.
.last_output_block <- function(lines) {
    rules <- grep("^-{10,}$", lines)
    if (length(rules) < 2L)
        return(NULL)
    n <- length(rules)
    block <- lines[seq.int(rules[n - 1L] + 1L, rules[n] - 1L)]
    block <- block[nzchar(trimws(block))]
    if (length(block)) block else NULL
}

#' Spatial analyses
#'
#' One function per spatial command, each a thin wrapper: the command name,
#' the options it takes and the table TDA writes. They all work on a
#' \code{\link{tda_spatial}} structure, however it was built -- from a data
#' frame with \code{tda_spatial}, from a shapefile with
#' \code{\link{tda_read_shapefile}}, or from line segments with
#' \code{\link{tda_polygons}}.
#'
#' \code{tda_sd_info} reports what a structure holds: how many points, lines
#' and polygons, and the bounding box. \code{tda_sd_voronoi} builds a
#' Voronoi diagram and the Delaunay triangulation of the points;
#' \code{tda_sd_enclosing} the enclosing rectangle or the convex hull,
#' \code{tda_sd_neighbours} nearest-neighbour lists, and
#' \code{tda_sd_relations} the relations between objects.
#' \code{tda_sd_select} picks out objects in a region -- like
#' \code{tda_sd_clip}, it needs \code{rec} (as an option:
#' \code{tda_sd_select(s, rec = "0,0,1,1")}), which TDA will not guess. The
#' polygon operations -- \code{tda_sd_intersect}, \code{tda_sd_lines},
#' \code{tda_sd_points} and \code{tda_sd_clip} -- each write a record file
#' that comes back as a table.
#'
#' @param s a \code{\link{tda_spatial}}.
#' @param shape for \code{tda_sd_enclosing}, \code{"rectangle"} (default)
#'   or \code{"hull"} (the convex hull); \code{sdencl}'s \code{opt=}.
#' @param ... options for the command.
#' @return An object carrying the run and, where there is one, a table.
#' @family spatial analysis
#' @examples
#' # four points, as a spatial structure
#' d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
#' s <- tda_spatial(d)
#'
#' # what is in it
#' tda_sd_info(s)
#'
#' # the Voronoi diagram of those points: what sdvd writes is itself a
#' # spatial data file (points for the Voronoi vertices, lines for the
#' # edges), so plot() reads it back with the input's plot machinery
#' v <- tda_sd_voronoi(s, opt = 2)
#' plot(v, xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5))
#'
#' # clip to a rectangle; rec is required, TDA will not guess one
#' tda_sd_clip(s, c(0, 0, 0.5, 0.5))
#'
#' # a select needs one too, the region to select within
#' tda_sd_select(s, rec = "0,0,1,1")$table
#'
#' # relations and nearest neighbours, on the same points
#' tda_sd_relations(s)
#' tda_sd_neighbours(s)
#' tda_sd_enclosing(s, shape = "hull")$table
#'
#' # polygon data: typ = 3 groups consecutive rows sharing an id into one
#' # polygon each -- here two squares
#' poly <- data.frame(id = c(1, 1, 1, 1, 2, 2, 2, 2),
#'                    x  = c(0, 2, 2, 0, 3, 5, 5, 3),
#'                    y  = c(0, 0, 2, 2, 0, 0, 2, 2))
#' sp <- tda_spatial(poly, type = 3)
#' tda_sd_data(sp)$table   # the polygons' own vertex coordinates back out
#'
#' # tda_sd_intersect/tda_sd_lines/tda_sd_points compare polygons against a
#' # *different* object type (points or lines) already in the same spatial
#' # structure -- with only the two polygons above and nothing else, there
#' # is nothing for them to find, which is what they correctly report
#' tda_sd_lines(sp)$table
#' tda_sd_intersect(sp)$table
#' tda_sd_points(sp)$table
#'
#' # a real shapefile, if one is available: sf's example data, North
#' # Carolina counties -- everything above works on synthetic points or two
#' # squares; this is what these look like on real, complicated polygons
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
#'   print(tda_sd_info(nc))                 # 108 counties, 2421 vertices
#'
#'   # relations' real answer is printed, not written to a file -- $text
#'   # holds the boundary comparisons TDA's console output shows
#'   head(tda_sd_relations(nc)$text, 4)
#'
#'   # a neighbour network needs line data; a set of polygons alone
#'   # legitimately has none to report -- 0 nodes, 0 edges, not a bug
#'   tda_sd_neighbours(nc)$text
#'
#'   tda_sd_enclosing(nc)$table             # the whole state's bounding box
#'
#'   # select/clip write one record per vertex, not all the same length --
#'   # a polygon's header row, then that many coordinate pairs, the
#'   # same ragged shape as sdvd's edge list above -- a proper rectangular
#'   # table is what tda_sd_data() gives instead
#'   tda_sd_select(nc, rec = "-79,35,-77,36")$table   # counties in a region
#'   tda_sd_clip(nc, c(-79, 35, -77, 36))$table       # the same, cut to it
#'   head(tda_sd_data(nc)$table)            # every county's vertices,
#'                                          # as an ordinary data frame
#' }
#' @name tda_sd_analyses
NULL

#' @rdname tda_sd_analyses
#' @export
tda_sd_info <- function(s, ...) tda_sd(s, "sdinf", list(...))

#' @rdname tda_sd_analyses
#' @export
tda_sd_relations <- function(s, ...)
    tda_sd(s, "sdrel", list(...), rhs = "out.txt")

#' @rdname tda_sd_analyses
#' @export
tda_sd_neighbours <- function(s, ...)
    tda_sd(s, "sdnl", list(...), rhs = "out.txt")

#' @rdname tda_sd_analyses
#' @export
tda_sd_enclosing <- function(s, shape = c("rectangle", "hull"), ...) {
    opts <- list(...)
    if (!missing(shape))
        opts$opt <- match(match.arg(shape), c("rectangle", "hull"))
    tda_sd(s, "sdencl", opts, rhs = "out.txt")
}

#' @rdname tda_sd_analyses
#' @export
tda_sd_voronoi <- function(s, ...)
    tda_sd(s, "sdvd", list(...), rhs = "out.txt")

# sdvd's output file (with opt = 1 or 2 -- the Voronoi diagram, not the
# Delaunay triangulation) is itself a valid spatial-data file: TDA's
# per-object header (SDID, SDType, SDN, SDPtr) that sdnvar/sdplot already
# read, holding the Voronoi vertices as points and the Voronoi edges as
# lines. Rather than parse that format by hand a second time, this points
# the *same* spatial structure's plot machinery (sdplot) at the new
# file and lets it do the reading, exactly as tda_sd_voronoi()'s
# `dir = s$dir` guarantees the file is sitting right where it was written.
#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::plot
plot.tda_sdvd <- function(x, ...) {
    vs <- x$spatial
    vs$file <- x$rhs
    plot(vs, ...)
}

#' Read a shapefile into a TDA spatial structure
#'
#' \code{sdshp} reads an ESRI shapefile -- the \code{.shp}, its \code{.shx}
#' index and the \code{.dbf} of attributes -- and writes TDA's spatial
#' data file, which is what the rest of the spatial commands work on.
#'
#' @param file path to the shapefile, with or without the \code{.shp}
#'   extension. The \code{.shx} and \code{.dbf} beside it are read too.
#' @param dir working directory for the run.
#' @param ... further options for \code{sdshp}.
#' @return A \code{tda_spatial} object, as \code{\link{tda_spatial}}
#'   returns.
#' @family spatial analysis
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
#'   print(tda_sd_info(nc))
#'   plot(nc)
#' }
#' @export
tda_read_shapefile <- function(file, dir = tempfile("tda"), ...) {
    base <- sub("\\.shp$", "", file, ignore.case = TRUE)
    have <- vapply(c("shp", "shx", "dbf"),
                   function(e) file.exists(paste0(base, ".", e)), NA)
    if (!all(have))
        stop("missing ", paste(names(have)[!have], collapse = ", "),
             " beside ", base)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    # TDA opens the three files itself from one stem, so all three are copied
    # in under a name it can find.
    for (e in c("shp", "shx", "dbf"))
        file.copy(paste0(base, ".", e), file.path(dir, paste0("shape.", e)),
                  overwrite = TRUE)
    # shape.sd is OUR OWN intermediate, and TDA's default print format
    # for it (12.6) costs six digits of every coordinate.  .sd_read()
    # was given fmt = "24.16" for exactly this reason when a real DCW
    # file exposed it; this reader builds its run and so never got
    # the fix -- nc.shp came back at 4 decimals (-84.3239 for
    # -84.32385254) until it did.  Overridable with fmt=.
    o <- .tda_extra(list(...))
    if (is.null(o$fmt))
        o$fmt <- "24.16"
    res <- tda_run(do.call(tda_block,
                           c(list(name = "sdshp"), o,
                             list(df = "shape.sd", rhs = "shape"))),
                   dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not read the shapefile: ", err[1L], call. = FALSE)
    n <- suppressWarnings(as.integer(regmatches(res$output,
        regexpr("(?<=Object-description records: )[0-9]+", res$output,
               perl = TRUE))))
    # The .dbf beside the geometry holds the attribute table -- one row
    # per shape, and for real data that is where the names live (nc.shp's
    # county names, for instance).  sdshp reads the geometry only, so the
    # attributes are read with TDA's rdbf and attached rather than
    # left for the caller to find.  A file with no readable attributes
    # simply gets NULL.
    att <- tryCatch(tda_read_dbf(file.path(dir, "shape.dbf")),
                    error = function(e) NULL)
    structure(list(call = sys.call(), run = res, dir = dir,
                   file = "shape.sd", data = NULL, attr_col = 6L,
                   attributes = att,
                   n = if (length(n)) n[1L] else NA_integer_),
              class = "tda_spatial")
}

#' Three-dimensional spatial plots
#'
#' \code{sdplot31}, \code{sdplot32} and \code{sdplot33} draw spatial objects
#' in a three-dimensional coordinate system rather than the flat one
#' \code{plot()} uses. That system comes from \code{psetup3}, so the ranges
#' of all three axes and the direction the plot is viewed from are set up
#' once with \code{tda_sd_ps3()}, the 3-D counterpart of \code{\link{tda_ps}}
#' -- \code{tda_sd_pl3_points()}, \code{tda_sd_pl3_lines()} and
#' \code{tda_sd_pl3_polygons()} then add drawing commands to it the way
#' \code{\link{tda_pl_lines}} and its siblings add to a \code{tda_ps}
#' session, and \code{plot()} draws the result. \code{\link{tda_pl_text3}}
#' and \code{\link{tda_pl_points3}} (TDA's 3-D-aware commands, not tied
#' to a spatial file) work on a \code{tda_sd_ps3} session too, which is how
#' a title or an extra marker gets added -- TDA's ordinary \code{plabel()}
#' (\code{\link{tda_pl_labels}}) is 2-D only and refuses outright with
#' \dQuote{current coordinate system is 3-dimensional} if tried here.
#' \code{\link{tda_pl_lines3}}, by contrast, does \strong{not} work on a
#' \code{tda_sd_ps3} session: it folds its coordinates into a regular
#' \code{nvar} variable, which needs an ordinary \code{nvar} declaration
#' that a spatial session, declared with \code{sdnvar} instead, does not
#' have -- \dQuote{Syntax error or undefined variables}. Use
#' \code{tda_pl_points3} instead, which takes literal coordinates and
#' needs no such declaration.
#'
#' \strong{A north-up view of geographic data} (longitude/latitude, the
#' way an ordinary map looks) is \code{view = c(-90, 90)}, not TDA's
#' own default of \code{c(30, 30)}. This is a property of TDA's
#' azimuth/elevation convention
#' itself, not of any particular dataset, so it is the same for any
#' longitude/latitude data.
#'
#' \strong{Each of the three drawing commands only draws one kind of
#' object}, silently drawing nothing at all for any other kind:
#' \code{tda_sd_pl3_points} draws points, \code{tda_sd_pl3_lines} draws
#' lines, \code{tda_sd_pl3_polygons} draws polygons. A shapefile of
#' counties, say, is polygons -- points on it draws nothing, not because
#' anything is wrong, but because it is looking for points that are not
#' there. Check \code{\link{tda_sd_info}} first if it is not obvious which
#' kind a structure holds.
#'
#' \code{tda_sd_plot3}, \code{tda_sd_plot3_lines} and
#' \code{tda_sd_plot3_polygons} are one-shot equivalents -- set up a
#' session, add one drawing command, plot it -- kept for the simple case
#' where nothing else needs adding to the same plot.
#'
#' @param s a \code{\link{tda_spatial}} structure.
#' @param p a \code{tda_sd_ps3} session.
#' @param cmd the TDA command name.
#' @param z,zlim the range of the z axis, as two numbers. Unlike x and y
#'   there is nothing in the spatial file to take it from, so it has no
#'   default.
#' @param xlim,ylim ranges of the x and y axes, taken from the structure
#'   when not given.
#' @param view direction of the projection, as longitude and latitude in
#'   degrees; TDA's default is 30,30.
#' @param width horizontal size of the plot in mm.
#' @param symbol marker symbol to draw at each point, as TDA's \code{s=}.
#'   It is named \code{symbol} rather than \code{s} because \code{s} is
#'   the spatial object.
#' @param size marker size in mm.
#' @param zvar an attribute carried in from the source file, used as each
#'   point's or polygon's height instead of a flat \code{zval} --
#'   for the polygon commands, each polygon is raised into a prism of
#'   that height rather than drawn flat. Name it as it appears in the
#'   object record -- for a shapefile, a column of the \code{.dbf}.
#' @param zval a single height for every point or polygon, when they do
#'   not have their own via \code{zvar}; TDA's default is 0.
#' @param newpage start a new grid page.
#' @param ... further options for the command.
#' @return \code{tda_sd_ps3()} and the \code{tda_sd_pl3_*} functions return
#'   the session, with the command added. \code{plot()} on a session, and
#'   the one-shot functions, draw the plot and return the parsed drawing
#'   operations invisibly.
#' @family spatial analysis
#' @examples
#' set.seed(1)
#' d <- data.frame(id = 1:8, x = runif(8), y = runif(8), z = runif(8))
#' s <- tda_spatial(d, attributes = "z")
#' pdf(NULL)
#'
#' # composable form: build a session, add commands, plot it -- a title
#' # (via TDA's 3-D-aware text command, not the 2-D-only plabel) and an
#' # extra marker layered onto the same plot as the data
#' p <- tda_sd_ps3(s, zlim = c(0, 1))
#' p <- tda_sd_pl3_points(p, symbol = 1, size = 3)
#' p <- tda_pl_text3(p, "eight points", at = c(0.5, 0.5, 0.9))
#' plot(p)
#'
#' # the one-shot form, for when nothing else needs adding
#' tda_sd_plot3(s, z = c(0, 1), symbol = 1, size = 3)
#'
#' # plot3_lines wants line-type spatial data, not scattered points --
#' # a single connected line here, four vertices sharing one id
#' ln <- data.frame(id = 1, x = c(0, 0.3, 0.6, 1), y = c(0, 0.4, 0.5, 1),
#'                  z = c(0, 0.3, 0.6, 1))
#' sl <- tda_spatial(ln, type = 2, attributes = "z")
#' tda_sd_plot3_lines(sl, z = c(0, 1), symbol = 1)
#'
#' # plot3_polygons draws flat polygons at a single height in the z range,
#' # or, given zvar, raises each one into a prism of that height -- h here
#' poly <- data.frame(id = 1, x = c(0, 2, 2, 0), y = c(0, 0, 2, 2),
#'                    h = c(2, 2, 2, 2))
#' sp <- tda_spatial(poly, type = 3, attributes = "h")
#' tda_sd_plot3_polygons(sp, z = c(0, 3), zvar = "h")
#'
#' # a real shapefile, if one is available: sf's example data, North
#' # Carolina counties, in a proper north-up view, with city markers and
#' # labels layered on via TDA's 3-D-aware point/text commands
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
#'   cities <- data.frame(
#'     name = c("Raleigh", "Charlotte"),
#'     lon = c(-78.6382, -80.8431), lat = c(35.7796, 35.2271))
#'
#'   pnc <- tda_sd_ps3(nc, zlim = c(0, 1), view = c(-90, 90))
#'   pnc <- tda_sd_pl3_polygons(pnc)
#'   # fs= is the marker/font size in mm, independent for each command --
#'   # the default marker size (2mm) reads as a large blob at this map's
#'   # own small default width (100mm), so it is turned down here, and
#'   # the default text size turned up, rather than using both defaults
#'   pnc <- tda_pl_points3(pnc, cities$lon, cities$lat, rep(0.5, 2),
#'                         symbol = 5, size = 0.8, lty = 0)
#'   for (i in seq_len(nrow(cities)))
#'     pnc <- tda_pl_text3(pnc, cities$name[i],
#'                         at = c(cities$lon[i] + 0.35, cities$lat[i], 0.5),
#'                         fs = 4)
#'   plot(pnc)
#' }
#' dev.off()
#' @export
tda_sd_ps3 <- function(s, xlim = NULL, ylim = NULL, zlim, view = c(30, 30),
                       width = 100, ...) {
    stopifnot(inherits(s, "tda_spatial"))
    if (missing(zlim) || length(zlim) != 2L)
        stop("`zlim` is required: the range of the z axis, as two numbers")
    if (is.null(xlim) || is.null(ylim)) {
        b <- .sd_bounds(s)
        if (is.null(xlim)) xlim <- b$xlim
        if (is.null(ylim)) ylim <- b$ylim
    }
    structure(list(s = s, xlim = xlim, ylim = ylim, zlim = zlim, view = view,
                   width = width, setup_opts = list(...), zvar = character(),
                   file = "spatial.ps", cmds = character(), data = NULL,
                   xlab = character(), xname = character(), run = NULL),
              class = "tda_sd_ps3")
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_pl3 <- function(p, cmd, ..., zvar = NULL) {
    stopifnot(inherits(p, "tda_sd_ps3"))
    opts <- list(...)
    # a zvar has to be named in the sdnvar block before it can be used --
    # and TDA variable names must start uppercase (see .sd_open() and
    # ?tda_write_data's "Variable names" section), which the declaration
    # and this reference both have to agree on or one silently points at
    # a variable that was never actually declared.
    if (!is.null(zvar)) {
        zvar <- .tda_names(as.character(zvar))
        opts$zvar <- zvar
        p$zvar <- union(p$zvar, zvar)
    }
    p$cmds <- c(p$cmds, do.call(tda_block, c(list(name = cmd), .tda_extra(opts))))
    p
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_pl3_points <- function(p, symbol = 1, size = NULL, zvar = NULL,
                              zval = NULL, ...) {
    opts <- list(...)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(zval)) opts$zval <- zval
    if (is.null(opts$s))
        stop("give `symbol`; without one sdplot31 draws nothing")
    do.call(tda_sd_pl3, c(list(p, "sdplot31"), opts, list(zvar = zvar)))
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_pl3_lines <- function(p, symbol = NULL, size = NULL, zvar = NULL,
                             zval = NULL, ...) {
    opts <- list(...)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(zval)) opts$zval <- zval
    if (is.null(opts$s))
        stop("give `symbol`; without one sdplot32 draws nothing")
    do.call(tda_sd_pl3, c(list(p, "sdplot32"), opts, list(zvar = zvar)))
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_pl3_polygons <- function(p, symbol = NULL, size = NULL, zvar = NULL,
                                zval = NULL, ...) {
    opts <- list(...)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(zval)) opts$zval <- zval
    do.call(tda_sd_pl3, c(list(p, "sdplot33"), opts, list(zvar = zvar)))
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::plot
plot.tda_sd_ps3 <- function(x, newpage = TRUE, ...) {
    res <- tda_run(c(.sd_open(x$s, x$zvar),
                     sprintf("psfile = %s;", x$file),
                     do.call(tda_block, c(list(name = "psetup3"), x$setup_opts,
                                          list(view = paste(x$view, collapse = ","),
                                              pxa = paste(x$xlim, collapse = ","),
                                              pya = paste(x$ylim, collapse = ","),
                                              pza = paste(x$zlim, collapse = ","),
                                              pxlen = x$width))),
                     x$cmds,
                     "psclose;"),
                   dir = x$s$dir)
    err <- grep("^Error|^Syntax", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not draw this: ", err[1L], call. = FALSE)
    tda_plot_ps(tda_read_ps(res, which = x$file), newpage = newpage)
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_plot3 <- function(s, z, xlim = NULL, ylim = NULL,
                         view = c(30, 30), width = 100, newpage = TRUE,
                         symbol = 1, size = NULL, zvar = NULL, zval = NULL,
                         ...) {
    p <- tda_sd_ps3(s, xlim = xlim, ylim = ylim, zlim = z, view = view,
                    width = width)
    p <- tda_sd_pl3_points(p, symbol = symbol, size = size, zvar = zvar,
                           zval = zval, ...)
    plot(p, newpage = newpage)
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_plot3_lines <- function(s, z, xlim = NULL, ylim = NULL,
                               view = c(30, 30), width = 100, newpage = TRUE,
                               symbol = NULL, size = NULL, zvar = NULL,
                               zval = NULL, ...) {
    p <- tda_sd_ps3(s, xlim = xlim, ylim = ylim, zlim = z, view = view,
                    width = width)
    p <- tda_sd_pl3_lines(p, symbol = symbol, size = size, zvar = zvar,
                          zval = zval, ...)
    plot(p, newpage = newpage)
}

#' @rdname tda_sd_ps3
#' @export
tda_sd_plot3_polygons <- function(s, z, xlim = NULL, ylim = NULL,
                                  view = c(30, 30), width = 100,
                                  newpage = TRUE, symbol = NULL,
                                  size = NULL, zvar = NULL, zval = NULL,
                                  ...) {
    p <- tda_sd_ps3(s, xlim = xlim, ylim = ylim, zlim = z, view = view,
                    width = width)
    p <- tda_sd_pl3_polygons(p, symbol = symbol, size = size, zvar = zvar,
                             zval = zval, ...)
    plot(p, newpage = newpage)
}

#' Maps
#'
#' \code{tda_map} draws spatial objects in a geographical coordinate system
#' -- a projection, rather than the plain x and y of \code{plot()}. TDA
#' builds that system with \code{psetupg} and then draws with \code{sdpmap};
#' \code{graticule} adds meridians and parallels with \code{sdpgrat}, and
#' \code{points} adds coordinates given directly with \code{sdpgeo}.
#'
#' \code{region} is required and has no default: it is the half-width and
#' half-height of the mapped area in degrees, measured from \code{view}. So
#' \code{view = c(-80, 35), region = c(6, 3)} maps 6 degrees of longitude
#' and 3 of latitude either side of that point.
#'
#' @param s a \code{\link{tda_spatial}} structure.
#' @param view centre of the projection, as longitude and latitude.
#' @param region half-width and half-height of the mapped area, in degrees.
#' @param projection 10 cylindrical equidistant, 11 cylindrical equal-area,
#'   12 Mercator, 20 azimuthal orthographic. For projection 20 the
#'   longitude half-width is limited to 90 degrees, and the whole region
#'   rectangle must fit on the projection disc, so both halves of
#'   \code{region} need \code{sin(lon)^2 + sin(lat)^2 <= 1} at the view
#'   latitude.
#' @param width width of the plot in mm.
#' @param graticule a list of options for the meridians and parallels, or
#'   \code{NULL} for none. \code{lon} and \code{lat} take TDA's sequence
#'   notation, \code{"-84(2)-76"} meaning every 2 degrees from -84 to -76.
#' @param points a data frame of \code{lon} and \code{lat} to mark, or
#'   \code{NULL}.
#' @param symbol marker symbol for \code{points}; TDA's \code{s=}, which
#'   cannot be spelled that way here because \code{s} is the structure. 4 is
#'   an open circle, 5 a filled one, 8 a square, 15 a diamond.
#' @param newpage start a new grid page.
#' @param ... further options for \code{sdpmap}.
#' @return The map, drawn; the session is returned invisibly.
#' @family spatial analysis
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
#'   tda_map(nc, view = c(-80, 35), region = c(6, 3),
#'           graticule = list(lon = "-84(2)-76", lat = "34(1)37",
#'                            fsx = 2, fsy = 2))
#' }
#' @export
tda_map <- function(s, view, region, projection = 10, width = 150,
                    graticule = list(), points = NULL, symbol = 1,
                    newpage = TRUE, ...) {
    if (missing(view) || length(view) != 2L)
        stop("`view` is required: the centre, as longitude and latitude")
    if (missing(region) || length(region) != 2L)
        stop("`region` is required and has no default: how far the map ",
             "reaches either side of `view`, in degrees")
    cmds <- c(.sd_open(s), "psfile = map.ps;",
              tda_block("psetupg", proj = projection,
                        view = paste(view, collapse = ","),
                        region = paste(region, collapse = ","),
                        pxlen = width))
    if (!is.null(graticule)) {
        g <- .tda_extra(graticule)
        if (is.null(g$lw))
            g$lw <- 0.1
        cmds <- c(cmds, do.call(tda_block, c(list(name = "sdpgrat"), g)))
    }
    cmds <- c(cmds, do.call(tda_block,
                            c(list(name = "sdpmap"), .tda_extra(list(...)))))
    if (!is.null(points)) {
        p <- as.data.frame(points)
        # sdpgeo takes the coordinates on the right-hand side, in pairs
        cmds <- c(cmds, tda_block("sdpgeo", s = symbol,
                                  rhs = paste(as.vector(t(as.matrix(
                                      p[, 1:2, drop = FALSE]))),
                                      collapse = ",")))
    }
    res <- tda_run(c(cmds, "psclose;"), dir = s$dir)
    err <- grep("^Error|^Syntax", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not draw this map: ", err[1L], call. = FALSE)
    tda_plot_ps(tda_read_ps(res, which = "map.ps"), newpage = newpage)
}

#' Build polygons from line segments
#'
#' \code{sdcpol} takes a set of line segments, finds where they meet and
#' assembles the polygons they enclose. It works on a plain table of
#' segments rather than on a \code{\link{tda_spatial}} structure, and with
#' \code{option = 3} it writes a spatial data file the rest of the family
#' can read.
#'
#' @param segments a data frame of four columns: the two ends of each
#'   segment, as x1, y1, x2, y2.
#' @param option 1 prints the node list of the graph the segments form, 2
#'   its edge list, 3 constructs the polygons.
#' @param dir working directory for the run.
#' @param ... further options for \code{sdcpol}; \code{tol} is the
#'   tolerance for treating two points as the same, default 1e-4.
#' @return For \code{option = 3} a \code{tda_spatial} structure; otherwise
#'   the run, with the table in \code{$table}.
#' @family spatial analysis
#' @examples
#' sq <- data.frame(x1 = c(0, 1, 1, 0), y1 = c(0, 0, 1, 1),
#'                  x2 = c(1, 1, 0, 0), y2 = c(0, 1, 1, 0))
#' tda_polygons(sq)
#' @export
tda_polygons <- function(segments, option = 3, dir = tempfile("tda"), ...) {
    d <- as.data.frame(segments)
    if (ncol(d) < 4L)
        stop("`segments` needs four columns: x1, y1, x2, y2")
    # sdcpol names the columns on its right-hand side, and TDA's variable
    # names are case-sensitive, so they are set here rather than taken from
    # whatever the caller happened to call them.
    d <- d[, 1:4]
    names(d) <- c("X1", "Y1", "X2", "Y2")
    out <- if (option == 3) "poly.sd" else "poly.out"
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block,
                             c(list(name = "sdcpol"), list(opt = option),
                               .tda_extra(list(...)),
                               list(df = out, rhs = "X1,Y1,X2,Y2")))),
                   data = d, dir = dir)
    err <- grep("^Error|^Syntax", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not build the polygons: ", err[1L], call. = FALSE)
    if (option != 3)
        return(structure(list(call = sys.call(), run = res, command = "sdcpol",
                              table = tryCatch(tda_file(res, out),
                                               error = function(e) NULL)),
                         class = c("tda_sdcpol", "tda_table")))
    structure(list(call = sys.call(), run = res, dir = dir, file = out,
                   data = NULL), class = "tda_spatial")
}

# The remaining format readers.  Each takes its input file on the
# right-hand side and writes TDA's spatial data file with df=, exactly
# as sdshp does, so one helper covers all of them.
.sd_read <- function(cmd, file, ext, dir, opts) {
    if (!file.exists(file))
        stop("no such file: ", file)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- paste0("input", ext)
    file.copy(file, file.path(dir, stem), overwrite = TRUE)
    # shape.sd is OUR OWN intermediate, not something the caller asked
    # for, and TDA's default print format for it (12.6) silently costs
    # six digits of every coordinate: reading a real DCW outline back
    # gave 30.271786 for a longitude of 30.271785736084, and every
    # spatial command downstream then worked from the rounded value.
    # A wide format costs nothing here -- the file is written once and
    # read once -- so it is the default, overridable with fmt=.
    o <- .tda_extra(opts)
    if (is.null(o$fmt))
        o$fmt <- "24.16"
    res <- tda_run(do.call(tda_block,
                           c(list(name = cmd), o,
                             list(df = "shape.sd", rhs = stem))),
                   dir = dir)
    err <- grep("^Error|^Syntax", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not read this with ", cmd, ": ", err[1L],
             call. = FALSE)
    # Object count, for $n -- phrased differently by each reader (sdshp's
    # own "Object-description records: N" vs sde00/sdgshhs/sddcwp's
    # "Number of points/lines/polygons: N", sometimes several of the
    # latter for a mixed file), so every phrasing actually seen is tried
    # and, for the "Number of ..." case, summed across geometry types.
    n <- suppressWarnings(as.integer(regmatches(res$output,
        regexpr("(?<=Object-description records: )[0-9]+", res$output,
               perl = TRUE))))
    if (!length(n)) {
        m <- regmatches(res$output,
                        regexpr("Number of (points|lines|polygons): [0-9]+",
                               res$output, perl = TRUE))
        if (length(m))
            n <- sum(as.integer(sub(".*: ", "", m)), na.rm = TRUE)
    }
    structure(list(call = sys.call(-1L), run = res, dir = dir,
                   file = "shape.sd", data = NULL, attr_col = 6L,
                   n = if (length(n)) n[1L] else NA_integer_),
              class = "tda_spatial")
}

#' Read other spatial formats
#'
#' The same shape as \code{\link{tda_read_shapefile}}: TDA reads the file
#' and writes its spatial data file, which the rest of the family works
#' on. \code{tda_read_e00} reads an uncompressed ARC/INFO export,
#' \code{tda_read_gshhs} a GSHHS coastline file, and \code{tda_read_dcw} a
#' Digital Chart of the World polygon point file.
#'
#' All three target older formats from GIS's pre-shapefile era, and only
#' shapefile really remains in wide, current use -- e00 and DCW were both
#' effectively superseded by it, and DCW's distribution site has been
#' gone for years, which is worth knowing going in:
#' \itemize{
#'   \item \code{tda_read_e00} works, tested against a real ARC/INFO
#'     export (a point coverage).
#'   \item \code{tda_read_dcw} works, tested against a \code{.pnt} file
#'     rebuilt from a real DCW coverage read via the \pkg{sf} package
#'     (\code{st_read()} on the original coverage directory, coordinates
#'     regrouped into rings and written out in the plain text form
#'     \code{sddcwp} expects). Coordinates round-trip exactly: the
#'     intermediate this reader writes uses \code{fmt = "24.16"}, not
#'     TDA's default, which would cost six digits of each one.
#'   \item \code{tda_read_gshhs} reads both GSHHS header layouts: the
#'     pre-2.2 one and GSHHG 2.x, which is what you would download
#'     today. The layout is detected from the file itself and the
#'     version reported in the run's output.
#'
#'     \strong{Straight lines across a world map} would come from the
#'     dateline; that is handled for you, and there is nothing to
#'     configure. GSHHS stores
#'     the polygons that span the dateline with longitudes running past
#'     180 (Eurasia reaches 190). Those are cut at the dateline and come
#'     back as separate objects, so each piece is drawn where it belongs
#'     and nothing is drawn between them; a polygon that does not span
#'     it is untouched. The Antarctic coastline, which runs from one
#'     edge of the map to the other, is returned as a line rather than a
#'     closed area for the same reason -- closing it would draw a
#'     straight segment back across the world.
#'
#'     This is what the GSHHG maintainers do in their own shapefile
#'     distribution, where the dateline-straddling polygons -- the
#'     Antarctic cap chief among them -- are split into east and west
#'     components. The native binary files are left unsplit for the
#'     reader to handle.
#'
#'     \code{centre = "pacific"} cuts the map at Greenwich instead,
#'     giving longitudes on 0..360; the splitting then happens at that
#'     seam rather than at the dateline. Either way no polygon is drawn
#'     across the map.
#' }
#'
#' @section Further options:
#' Anything else these commands accept can be passed through \code{...}
#' under TDA's option name. \code{tda_help("sdshp")},
#' \code{tda_help("sde00")}, \code{tda_help("sddcwp")} and
#' \code{tda_help("sdgshhs")} print the full lists. One worth knowing is
#' \code{nc}, which tells the DCW reader that the input has no comment
#' lines.
#'
#' \code{df=}, \code{dtda=} and \code{fmt=} are set by the reader
#' itself: the first two name the intermediate files it writes and reads
#' back, and \code{fmt} is fixed at \code{"24.16"} so that no
#' coordinate is rounded on the way through.
#'
#' @param centre for \code{tda_read_gshhs}: where the world map is cut.
#'   \code{"greenwich"} (default) gives the familiar view, longitudes
#'   -180..180, cut at the dateline. \code{"pacific"} cuts at the
#'   Greenwich meridian instead, longitudes 0..360. Polygons spanning
#'   whichever seam is in force are split either way.
#' @param file path to the input file.
#' @param level for \code{tda_read_gshhs}: which hierarchy levels to
#'   keep -- 1 land, 2 lake, 3 island in lake, 4 pond in island in lake,
#'   and since GSHHG 2.3.0 5 the Antarctic ice front and 6 the Antarctic
#'   grounding line. Several may be given, \code{level = c(1, 2)} for
#'   land and lakes. \code{NULL}, the default, keeps all of them; level
#'   1 alone is much the fastest way to draw a coastline.
#'
#'   \strong{Antarctica appears twice} in a 2.3.x file, once as level 5
#'   and once as level 6, and they are about six degrees of latitude
#'   apart near the dateline. Keeping all levels draws both, which is
#'   why a world map shows two near-parallel southern coastlines; GSHHG
#'   intends you to pick one. Neither is an artefact -- the ice front
#'   runs along roughly 78 degrees south there (the Ross Ice Shelf edge,
#'   near-straight) and the grounding line along roughly 84.
#' @param dir working directory for the run.
#' @param ... further options for the command. \code{sde00} takes
#'   \code{attr} for how many attribute variables to keep; \code{sdgshhs}
#'   takes \code{level} to select polygons; \code{sddcwp} takes \code{nc}
#'   to drop centroids.
#' @return A \code{tda_spatial} structure.
#' @family spatial analysis
#' @examples
#' # nc.shp ships with sf, so this one runs where sf is installed.
#' # The .dbf beside the geometry holds the attribute table -- county
#' # names and areas -- and is attached as $attributes.
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
#'   head(nc$attributes[, c("NAME", "AREA")])
#' }
#'
#' \dontrun{
#' # the other readers need files that do not ship with any package
#' tda_read_e00("world.e00")
#' tda_read_gshhs("gshhs_f.b", level = 1)  # either header layout
#' tda_read_dcw("polygon.pnt", nc = TRUE)
#' }
#' @name tda_read_spatial
NULL

#' @rdname tda_read_spatial
#' @export
tda_read_e00 <- function(file, dir = tempfile("tda"), ...)
    .sd_read("sde00", file, ".e00", dir, list(...))

#' @rdname tda_read_spatial
#' @export
tda_read_gshhs <- function(file, level = NULL, centre = c("greenwich", "pacific"),
                           dir = tempfile("tda"), ...) {
    o <- list(...)
    if (!is.null(level)) {
        if (!is.numeric(level) || anyNA(level) ||
            any(level < 1) || any(level > 6) || any(level != trunc(level)))
            stop("`level` must be whole numbers between 1 and 6", call. = FALSE)
        level <- sort(unique(as.integer(level)))
        # TDA's level= takes one level.  A set is passed as a bitmask
        # offset by 64 -- bit 0 for level 1, and so on -- which sdgshhs
        # unpacks; a single level still goes through as itself, so the
        # command file reads level=1 for the common case rather than
        # level=65.
        o$level <- if (length(level) == 1L) level
                   else 64L + sum(bitwShiftL(1L, level - 1L))
    }
    # TDA spells this opt=, 1 "standard coding of longitudes" and 2
    # "longitudes unchanged", neither of which says what it does to a
    # map.  What it picks is where the map is cut: at the dateline
    # (longitudes -180..180, the familiar view) or at the Greenwich
    # meridian (0..360, centred on the Pacific).  Polygons spanning
    # whichever seam is in force are split by sdgshhs either way.
    if (!is.null(o$opt))
        stop("use centre = \"greenwich\" or \"pacific\" rather than TDA's opt=",
             call. = FALSE)
    o$opt <- if (match.arg(centre) == "greenwich") 1L else 2L
    .sd_read("sdgshhs", file, ".b", dir, o)
}

#' @rdname tda_read_spatial
#' @export
tda_read_dcw <- function(file, dir = tempfile("tda"), ...)
    .sd_read("sddcwp", file, ".pnt", dir, list(...))

#' @rdname tda_sd_analyses
#' @export
tda_sd_select <- function(s, ...) tda_sd(s, "sdsel", list(...), rhs = "out.txt")

# The polygon operations.  Each writes a record file, so the table comes
# back the same way the others do.

#' @rdname tda_sd_analyses
#' @export
tda_sd_intersect <- function(s, ...)
    tda_sd(s, "sdipol", list(...), rhs = "out.txt")

#' @rdname tda_sd_analyses
#' @export
tda_sd_lines <- function(s, ...)
    tda_sd(s, "sdlpol", list(...), rhs = "out.txt")

#' @rdname tda_sd_analyses
#' @export
tda_sd_points <- function(s, ...)
    tda_sd(s, "sdppol", list(...), rhs = "out.txt")

#' @rdname tda_sd_analyses
#' @export
tda_sd_data <- function(s, ...)
    tda_sd(s, "sdpdata", list(...), rhs = "out.txt")

# sdclip needs the rectangle to clip to; TDA refuses without it, so it is a
# named argument rather than something to discover from an error message.
#' @rdname tda_sd_analyses
#' @param rec the clipping rectangle, \code{xmin,ymin,xmax,ymax}.
#' @export
tda_sd_clip <- function(s, rec, ...) {
    if (missing(rec))
        stop("`rec` is required: the rectangle to clip to, as ",
             "xmin,ymin,xmax,ymax")
    if (is.numeric(rec))
        rec <- paste(rec, collapse = ",")
    tda_sd(s, "sdclip", c(list(rec = rec), list(...)), rhs = "out.txt")
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_spatial <- function(x, ...) {
    # Two kinds of object share this class: one built from an in-memory
    # data frame (tda_spatial(), which keeps $data), and one built by a
    # file reader (shapefile/e00/DCW/GSHHS, where $data is NULL and the
    # information lives in the reader's console output).  A header
    # alone told the user nothing about what was actually read; show
    # the coordinates' shape either way.
    cat("TDA spatial data (", x$file, ")\n", sep = "")
    if (is.data.frame(x$data) && nrow(x$data)) {
        d <- x$data
        idc <- d[[1L]]
        xs <- d[[2L]]
        ys <- d[[3L]]
        cat("  ", nrow(d), " coordinate rows, ",
            length(unique(idc)), " object",
            if (length(unique(idc)) != 1L) "s",
            ", x in [", format(min(xs)), ", ", format(max(xs)),
            "], y in [", format(min(ys)), ", ", format(max(ys)), "]\n",
            sep = "")
        cat("  columns:", paste(x$xlab, collapse = ", "), "\n\n")
        if (length(x$xlab) == ncol(d))
            names(d) <- x$xlab
        print(utils::head(d, 6L))
        if (nrow(d) > 6L)
            cat("  ...", nrow(d) - 6L, "more rows\n")
    } else {
        if (!is.null(x$n) && !is.na(x$n))
            cat("  ", x$n, " object-description records\n", sep = "")
        # the reader's summary: object counts and the bounding box
        info <- grep(paste0("^(Number of (points|lines|polygons|arcs|",
                            "records|objects)|Object-description|",
                            "XMin|XMax|YMin|YMax)"),
                     x$run$output, value = TRUE)
        if (length(info))
            cat(paste0("  ", info), sep = "\n")
        else
            cat("  (no summary in the reader's output; see $run$output)\n")
    }
    invisible(x)
}

#' Draw spatial data
#'
#' \code{sdplot} needs a PostScript coordinate system as well as the spatial
#' data, so the ranges have to be given.
#'
#' Points need a marker or nothing is drawn: give \code{symbol} (or another
#' TDA symbol number). This is the same trap noted for
#' \code{\link{tda_sd_plot3}} -- without a marker or a \code{zvar} attribute,
#' \code{sdplot} draws nothing at all rather than erroring.
#'
#' @param x a \code{\link{tda_spatial}}.
#' @param xlim,ylim ranges of the coordinate system.
#' @param width,height size of the plotting area in millimetres; the
#'   height follows the data's aspect ratio unless given.
#' @param newpage start a new page before drawing.
#' @param what what to draw: \code{"objects"} (default), the convex hull
#'   over every point, or the convex hull of each polygon separately;
#'   \code{sdplot}'s \code{opt=}.
#' @param symbol,size,lty,lw,gray marker symbol and size, line type, line
#'   width in mm, and grey level from 0 (black) to 1 (white); \code{sdplot}'s
#'   own \code{s=}/\code{fs=}/\code{lt=}/\code{lw=}/\code{gs=}.
#' @param ... further options for \code{sdplot}.
#' @return The parsed drawing operations, invisibly.
#' @family spatial analysis
#' @examples
#' d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
#' s <- tda_spatial(d)
#' plot(s, symbol = 1)  # a marker is needed; points draw nothing without one
#' @exportS3Method base::plot
plot.tda_spatial <- function(x, xlim = NULL, ylim = NULL, width = 100,
                             height = NULL, newpage = TRUE,
                             what = c("objects", "hull", "hull_per_polygon"),
                             symbol = NULL, size = NULL, lty = NULL,
                             lw = NULL, gray = NULL, ...) {
    # tda_spatial()-built structures carry the original coordinates in
    # $data, so the range can be read off them directly -- but a structure
    # from tda_read_shapefile() (or any other route that reads a spatial
    # file TDA already has, without going through R first) has $data =
    # NULL: there was never an R data frame to keep. range(NULL) silently
    # returns Inf/-Inf with a warning rather than an error, which produced
    # a degenerate xlim/ylim here and "nothing to draw" with no indication
    # why. .sd_bounds() falls back to tda_sd_info()'s report either way.
    if (is.null(xlim) || is.null(ylim)) {
        b <- .sd_bounds(x)
        if (is.null(xlim)) xlim <- b$xlim
        if (is.null(ylim)) ylim <- b$ylim
    }
    # TDA's psetup() scales x and y independently, so a map drawn into a
    # fixed 100 x 80 mm box is stretched to it.  Left unset, the height
    # follows the data: equal scales in both directions, within reason.
    if (is.null(height)) {
        r <- diff(ylim) / diff(xlim)
        height <- width * max(0.25, min(2, if (is.finite(r) && r > 0) r else 0.8))
    }
    opts <- list(...)
    if (!missing(what))
        opts$opt <- match(match.arg(what),
                          c("objects", "hull", "hull_per_polygon"))
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    res <- tda_run(c(.sd_open(x),
                     "psfile = spatial.ps;",
                     tda_block("psetup", pxlen = width, pylen = height,
                               pxa = paste(xlim, collapse = ","),
                               pya = paste(ylim, collapse = ",")),
                     do.call(tda_block, c(list(name = "sdplot"),
                                          .tda_extra(opts))),
                     "psclose;"),
                   dir = x$dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not draw this: ", err[1L], call. = FALSE)
    tda_plot_ps(tda_read_ps(res, which = "spatial.ps"), newpage = newpage)
}

#' Read a dBase (.dbf) attribute table
#'
#' \code{tda_read_dbf} reads a dBase file through TDA's \code{rdbf} and
#' returns it as a data frame. Character fields come back as character,
#' numeric ones as numeric, taken from the field types the file itself
#' declares.
#'
#' This is the attribute half of a shapefile. \code{\link{tda_read_shapefile}}
#' reads the geometry; the \code{.dbf} beside it holds one row per shape
#' with whatever columns the data has, which for real data is usually
#' where the names live. \code{\link{tda_read_shapefile}} attaches it as
#' \code{attr(x, "attributes")}.
#'
#' @param file path to a \code{.dbf} file.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame; \code{attr(x, "run")} carries the run.
#' @family spatial analysis
#' @examples
#' # a shapefile's .dbf, which sf ships alongside nc.shp
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   f <- sub("[.]shp$", ".dbf", system.file("shape/nc.shp", package = "sf"))
#'   head(tda_read_dbf(f)[, c("NAME", "AREA")])
#' }
#' @export
tda_read_dbf <- function(file, options = list(), dir = tempfile("tda"), ...) {
    if (!file.exists(file))
        stop("no such file: ", file)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- basename(file)
    file.copy(file, file.path(dir, stem), overwrite = TRUE)
    res <- tda_run(do.call(tda_block,
                           c(list(name = "rdbf"), .tda_extra(options),
                             list(df = "out.txt", rhs = stem))),
                   dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not read this with rdbf: ", err[1L], call. = FALSE)
    tab <- .rdbf_frame(res)
    if (is.null(tab))
        tab <- .rdbf_frame_file(res, dir)
    if (is.null(tab))
        stop("rdbf produced no data")
    attr(tab, "run") <- res
    tab
}

.rdbf_frame_file <- function(res, dir) {
    # flag-off fallback: the df file is strictly fixed-width -- each
    # field rendered at the Length from the protocol's field listing,
    # one space between fields -- so it is sliced by that map rather
    # than split on whitespace (a numeric-to-string boundary has only
    # a single space).
    f <- file.path(dir, "out.txt")
    if (!file.exists(f))
        return(NULL)
    ln <- readLines(f)
    ln <- ln[nzchar(trimws(ln))]
    if (!length(ln))
        return(NULL)
    i <- grep("^Idx +Type +Length", res$output)
    if (!length(i))
        return(NULL)
    body <- res$output[(i[1L] + 2L):length(res$output)]
    body <- body[seq_len(max(0L,
                 which(!grepl("^ *[0-9]+ ", body))[1L] - 1L))]
    nm  <- sub(".* (\\S+) *$", "\\1", body)
    ty  <- sub("^ *[0-9]+ +([0-9]+).*", "\\1", body)
    len <- as.integer(sub("^ *[0-9]+ +[0-9]+ +([0-9]+).*", "\\1", body))
    if (anyNA(len))
        return(NULL)
    pos <- cumsum(c(1L, len + 1L))
    d <- as.data.frame(lapply(seq_along(len), function(k)
        trimws(substr(ln, pos[k], pos[k] + len[k] - 1L))),
        stringsAsFactors = FALSE)
    names(d) <- nm
    for (k in seq_along(d))
        if (ty[k] != "1")
            d[[k]] <- as.numeric(d[[k]])
    d
}

# The attribute table from rdbf's exports: rdbf.fields is "<name>|<type>"
# per field (1 = character, 0 = numeric, 2 = anything else), rdbf.values
# every field of every record as text, in row-major order.  The values
# travel as text because a dbf holds character fields as well as numbers;
# the conversion happens here, per field, from the type the file declares
# rather than by guessing from the content.
.rdbf_frame <- function(res) {
    if (!.use_exports())
        return(NULL)
    f <- res$exports[["rdbf.fields"]]
    v <- res$exports[["rdbf.values"]]
    if (!is.character(f) || !is.character(v) || !length(f) || !length(v))
        return(NULL)
    parts <- strsplit(f, "|", fixed = TRUE)
    if (any(lengths(parts) != 2L))
        return(NULL)
    nm <- trimws(vapply(parts, `[`, "", 1L))
    ty <- vapply(parts, `[`, "", 2L)
    if (length(v) %% length(nm) != 0L)
        return(NULL)
    m <- matrix(v, ncol = length(nm), byrow = TRUE)
    cols <- lapply(seq_along(nm), function(j) {
        x <- trimws(m[, j])
        if (ty[j] == "0") suppressWarnings(as.numeric(x)) else x
    })
    names(cols) <- nm
    as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
}

