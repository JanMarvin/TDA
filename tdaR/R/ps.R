# Replaying TDA's PostScript output through grid.
#
# TDA draws with a small, fixed subset of PostScript: it defines its
# operators in a prolog -- m, l, rm, rl, the symbol procedures, center and
# adjust -- and then emits a stream using those plus stroke, fill, clip,
# gsave/grestore, setlinewidth, setdash, setgray and show.  That is a few
# dozen operators in total, not a general PostScript interpreter's worth, so
# reading it back is tractable.
#
# This is a stepping stone rather than the destination.  The goal is for each
# plot command to build grid output directly; replaying the PostScript proves
# the geometry is right and gives something to check that against, because the
# .ps file is what TDA itself considers correct.

# TDA's 17 numbered plot symbols (ps_sym() in t_plot.c): each is a shape,
# filled white ("invers"), filled black, or stroked as an outline only, and
# some additionally have a cross or x mark drawn on top. Not just "circle"
# vs "square": a filled-white shape with no mark on top (types 1/2 without
# their mark would be, though none actually go that far) is genuinely
# invisible against a white background, which is what happens if this
# table is not used and the shape and its mark are not both drawn.
.SYMBOL_TABLE <- list(
    `1`  = list(shape = "circle", fill = "white", mark = "cross"),
    `2`  = list(shape = "circle", fill = "white", mark = "xsym"),
    `3`  = list(shape = "circle", fill = "white", mark = "both"),
    `4`  = list(shape = "circle", fill = "white", outline = TRUE),
    `5`  = list(shape = "circle", fill = "black"),
    `6`  = list(shape = "circle", fill = "white", outline = TRUE, mark = "cross"),
    `7`  = list(shape = "circle", fill = "white", outline = TRUE, mark = "xsym"),
    `8`  = list(shape = "square", outline = TRUE),
    `9`  = list(shape = "square", fill = "black"),
    `10` = list(shape = "trian1", outline = TRUE),
    `11` = list(shape = "trian1", fill = "black"),
    `12` = list(shape = "trian2", outline = TRUE),
    `13` = list(shape = "trian2", fill = "black"),
    `14` = list(shape = "square", outline = TRUE, mark = "cross"),
    `15` = list(shape = "rhomb", outline = TRUE),
    `16` = list(shape = "rhomb", fill = "black"),
    `17` = list(shape = "rhomb", outline = TRUE, mark = "cross")
)

# Everything before the setup is the prolog: procedure definitions, not
# drawing.  "setlinejoin" appears there as well as in the setup, so the last
# occurrence is the start of the drawing section, not the first.
.ps_body <- function(txt) {
    # %%BeginProlog/%%EndProlog are standard PostScript structuring
    # comments bracketing the prolog (font/symbol procedure
    # definitions) -- exactly one pair for an ordinary, single-panel
    # plot, but dplot() (combining several already-created files into
    # a grid, examples/exam/plot7.cf) keeps each combined file's own,
    # complete, separate prolog rather than sharing one: confirmed
    # directly, a real four-panel combined file has four full
    # %%BeginProlog...%%EndProlog pairs, wrapped by dplot()'s
    # "gsave % added" / "grestore % added" isolation around each one
    # (which precede/follow the pair, not sit inside it -- that
    # ordering is why the first version of this fix, taking only the
    # first %%EndProlog as the single boundary, still broke: it
    # stripped the very first panel's "gsave % added" along with
    # its prolog, leaving that gsave's matching "grestore % added"
    # unmatched later and the panel's coordinate-system state
    # never actually isolated from the next panel's own, so a scale
    # factor compounded across panels, 0.43 becoming 0.43^2
    # for every panel after the first -- checked by tracing
    # the parser's transform stack). Removing every
    # %%BeginProlog...%%EndProlog pair found, not just content before
    # the first %%EndProlog, keeps each panel's surrounding
    # "gsave % added"/"grestore % added" markers intact while still
    # discarding all four panels' own, purely-boilerplate prolog text.
    begins <- grep("%%BeginProlog", txt)
    ends <- grep("%%EndProlog", txt)
    if (length(begins) && length(ends) &&
        length(begins) == length(ends)) {
        drop <- unlist(Map(function(b, e) seq.int(b, e), begins, ends))
        body <- txt[-drop]
    } else {
        # Fallback for a file with neither marker (not expected from
        # TDA itself, but defensive): the last "setlinejoin" is at
        # least closer to the true boundary than nothing at all.
        start <- grep("setlinejoin", txt)
        body <- txt[seq.int(if (length(start)) start[length(start)] + 1L
                            else 1L, length(txt))]
    }
    .ps_strip_procs(body)
}

# dplot() (combining several already-created PostScript files into a grid,
# examples/exam/plot7.cf) re-emits each combined file's full set of
# procedure definitions ("/sclear {...} def", "/center {...} def", and so
# on -- the same block %%EndProlog is meant to be the single, one-time
# boundary for) once per panel, not just at the very start of the file:
# checked, "/sclear {...} def" appears four separate times in a
# real four-panel combined file, evenly spaced, one per panel. A procedure
# *definition* is never executed line by line the way the rest of the
# drawing stream is -- but this parser has no notion of PostScript's
# "/name { ... } def" syntax and, without this, read every gsave/grestore/
# setgray/stroke *inside* such a block as if it were a real, sequential
# drawing command. That desynchronised this parser's gsave/grestore
# stack from the file's true nesting (checked: a real
# four-panel file came back with 24 gsave against only 20 grestore within
# what should have been one, self-contained panel), leaving every panel
# after the first stuck inside a colour state (white) that belonged to an
# unrelated, never-executed procedure body, not the one actually active in
# the real drawing stream at that point -- confirmed against Ghostscript's
# own rendering of the identical, unmodified file, which shows all four
# panels correctly. Stripped here by brace-depth matching (a nested
# "{...} if", as adjust's body has, needs more than a same-line
# open/close check) rather than assumed to appear only where the prolog
# itself does.
.ps_strip_procs <- function(body) {
    depth <- 0L
    in_def <- FALSE
    keep <- rep(TRUE, length(body))
    # One vectorised test instead of a grepl() per line.  A world map's
    # PostScript is ~38000 lines and holds about thirty procedure
    # definitions, so all but a handful of those calls found nothing;
    # Rprof put this function at a quarter of tda_read_ps()'s time.
    # The loop still has to run -- the brace depth is sequential state --
    # but it now only indexes a logical vector.
    is_start <- grepl("^/[A-Za-z_][A-Za-z0-9_]*\\s*\\{", body)
    for (i in seq_along(body)) {
        ln <- body[[i]]
        if (!in_def && is_start[i]) {
            in_def <- TRUE
            depth <- 0L
        }
        if (in_def) {
            keep[i] <- FALSE
            depth <- depth + lengths(regmatches(ln, gregexpr("\\{", ln))) -
                lengths(regmatches(ln, gregexpr("\\}", ln)))
            if (depth <= 0L)
                in_def <- FALSE
        }
    }
    body[keep]
}


# The graphics state carries the current transformation matrix, not just an
# offset.  TDA moves, turns and scales the coordinate system rather than
# emitting absolute coordinates in several places: pltext does
# "gsave / x y translate / 0 0 m / show", and an arrowhead is
# "x y m / a rotate / sx sy scale / rl rl rl / fill" -- so a parser that
# handles only translate draws every arrowhead pointing the same way.
#
# ctm is (a, b, c, d, e, f) as PostScript orders it: a point (x, y) maps to
# (a*x + c*y + e, b*x + d*y + f).
.ps_state <- function()
    list(lwd = 1, lty = "solid", col = "black", fontsize = 8.5,
         ssiz = 1, lpt = NULL, clip = NULL, ctm = c(1, 0, 0, 1, 0, 0),
         font = "FT")

# Apply the current transform to a point in user space.
.ps_map <- function(m, x, y)
    c(m[1L] * x + m[3L] * y + m[5L], m[2L] * x + m[4L] * y + m[6L])

# Concatenate a new transform onto the current one, as PostScript does:
# the new one applies first, then the old.
.ps_concat <- function(m, n)
    c(n[1L] * m[1L] + n[2L] * m[3L],
      n[1L] * m[2L] + n[2L] * m[4L],
      n[3L] * m[1L] + n[4L] * m[3L],
      n[3L] * m[2L] + n[4L] * m[4L],
      n[5L] * m[1L] + n[6L] * m[3L] + m[5L],
      n[5L] * m[2L] + n[6L] * m[4L] + m[6L])

.ps_num <- function(x) suppressWarnings(as.numeric(x))

# A PostScript string literal is delimited by parentheses, so a parenthesis
# inside one is written \( or \), and a backslash \\.  Reading a label back
# has to undo that or the escapes are drawn as part of the text -- ehe2p's
# "Piecewise Constant Exponential (I)" comes back with its brackets
# backslashed.  \ddd is octal and the rest are the usual control characters.
.ps_unescape <- function(x, font = "FT", spanvec_table = NULL) {
    esc <- c("(" = "(", ")" = ")", "\\" = "\\", "n" = "\n", "r" = "\r",
             "t" = "\t", "b" = "\b", "f" = "\f")
    vapply(x, function(s) {
        if (!grepl("\\\\", s)) return(s)
        ch <- strsplit(s, "")[[1L]]
        out <- character(0)
        i <- 1L
        while (i <= length(ch)) {
            if (ch[i] != "\\" || i == length(ch)) {
                out <- c(out, ch[i])
                i <- i + 1L
                next
            }
            nx <- ch[i + 1L]
            oct <- ch[seq(i + 1L, min(i + 3L, length(ch)))]
            oct <- oct[cumprod(grepl("^[0-7]$", oct)) == 1L]
            if (length(oct)) {
                # Which table an octal escape maps through depends on
                # which font is active when it is drawn -- see the
                # .PS_SYMBOL_FONT/.PS_AGL/.ps_spanvec() comments for
                # why neither is the plain Latin-1 byte value this
                # used to decode to regardless of font.
                key <- paste(oct, collapse = "")
                mapped <- if (identical(font, "FS"))
                    .PS_SYMBOL_FONT[[key]]
                else if (!is.null(spanvec_table) && key %in% names(spanvec_table))
                    spanvec_table[[key]]
                else
                    NULL
                if (is.null(mapped) || !nzchar(mapped)) {
                    mapped <- rawToChar(as.raw(strtoi(key, 8L)))
                    Encoding(mapped) <- "latin1"
                    mapped <- enc2utf8(mapped)
                }
                out <- c(out, mapped)
                i <- i + 1L + length(oct)
            } else if (nx %in% names(esc)) {
                out <- c(out, esc[[nx]])
                i <- i + 2L
            } else {
                # a backslash before anything else is dropped, as PostScript
                # drops it
                out <- c(out, nx)
                i <- i + 2L
            }
        }
        paste(out, collapse = "")
    }, character(1), USE.NAMES = FALSE)
}

# TDA writes special characters two different ways, both by number, in a
# text string: an octal escape ("\245") drawn in the current font, or an
# "@nnn" TDA-source form (which reaches this parser already turned into a
# raw octal escape, with the font switched to "FS" first, by the time the
# PostScript is written) -- checked by comparing a real .cf
# file's "@141" against its output, "/fsiz ... FS (\141) show").
# Which table an octal code maps through depends entirely on which font
# is active: "FT" (Times-Roman, re-encoded via spanvec, embedded in every
# TDA-generated file's prolog) or "FS" (the standard, fixed Adobe
# Symbol font, used for Greek letters and math symbols, not TDA-specific
# at all). Decoding every octal escape through .ps_unescape()'s plain
# Latin-1 byte value regardless of which font was active -- what this
# parser did before -- is wrong for both: spanvec's codes are not
# Latin-1 (its own 8#245 is a bullet, not the Latin-1 byte 0o245), and
# Symbol's codes were never Latin-1 to begin with (its own 8#141 is
# lower-case alpha, not the ASCII letter "a" that byte value would
# otherwise be). Checked: examples/exam/plots.cf's three
# lines of specially-encoded text rendered as raw, undecoded ASCII
# instead of "alpha = 2.5 (bullet) beta", "Sigma = { theta | theta
# (elementof) Theta }" and "Sm(oslash)rebr(oslash)d" before this.
#
# The standard Adobe Symbol font encoding, reached via "FS" -- fixed
# across every PostScript implementation that ships this font, not
# TDA-specific, and not the full 256-entry table, only the range TDA's
# own "@nnn" text actually reaches in practice (lower-case and
# upper-case Greek letters at the same octal positions as their Latin
# look-alikes, plus a handful of common math symbols).
.PS_SYMBOL_FONT <- c(
    # Row 4 (0100-0117): upper-case Greek, plus congruent at 0100.
    `100` = "\u2245", `101` = "\u0391", `102` = "\u0392", `103` = "\u03a7",
    `104` = "\u0394", `105` = "\u0395", `106` = "\u03a6", `107` = "\u0393",
    `110` = "\u0397", `111` = "\u0399", `112` = "\u03d1", `113` = "\u039a",
    `114` = "\u039b", `115` = "\u039c", `116` = "\u039d", `117` = "\u039f",
    # Row 5 (0120-0137): more upper-case Greek, brackets, therefore/perp.
    `120` = "\u03a0", `121` = "\u0398", `122` = "\u03a1", `123` = "\u03a3",
    `124` = "\u03a4", `125` = "\u03a5", `126` = "\u03c2", `127` = "\u03a9",
    `130` = "\u039e", `131` = "\u03a8", `132` = "\u0396", `133` = "[",
    `134` = "\u2234", `135` = "]", `136` = "\u22a5", `137` = "_",
    # Row 6 (0140-0157): lower-case Greek.
    `140` = "\u203e", `141` = "\u03b1", `142` = "\u03b2", `143` = "\u03c7",
    `144` = "\u03b4", `145` = "\u03b5", `146` = "\u03c6", `147` = "\u03b3",
    `150` = "\u03b7", `151` = "\u03b9", `152` = "\u03d5", `153` = "\u03ba",
    `154` = "\u03bb", `155` = "\u03bc", `156` = "\u03bd", `157` = "\u03bf",
    # Row 7 (0160-0177): more lower-case Greek, braces, tilde, CR arrow.
    `160` = "\u03c0", `161` = "\u03b8", `162` = "\u03c1", `163` = "\u03c3",
    `164` = "\u03c4", `165` = "\u03c5", `166` = "\u03d6", `167` = "\u03c9",
    `170` = "\u03be", `171` = "\u03c8", `172` = "\u03b6", `173` = "{",
    `174` = "|", `175` = "}", `176` = "\u223c", `177` = "\u21b5",
    # Row A (0240-0257): currency, comparisons, arrows.
    `240` = "\u20ac", `241` = "\u03d2", `242` = "\u2032", `243` = "\u2264",
    `244` = "\u2044", `245` = "\u221e", `246` = "\u0192", `247` = "\u2663",
    `250` = "\u2666", `251` = "\u2665", `252` = "\u2660", `253` = "\u2194",
    `254` = "\u2190", `255` = "\u2191", `256` = "\u2192", `257` = "\u2193",
    # Row B (0260-0277): degree, comparisons, calculus.
    `260` = "\u00b0", `261` = "\u00b1", `262` = "\u2033", `263` = "\u2265",
    `264` = "\u00d7", `265` = "\u221d", `266` = "\u2202", `267` = "\u2022",
    `270` = "\u00f7", `271` = "\u2260", `272` = "\u2261", `273` = "\u2248",
    `274` = "\u2026", `276` = "\u21b5", `277` = "\u21b5",
    # Row C (0300-0317): set theory.
    `300` = "\u2135", `301` = "\u2111", `302` = "\u211c", `303` = "\u2118",
    `304` = "\u2297", `305` = "\u2295", `306` = "\u2205", `307` = "\u2229",
    `310` = "\u222a", `311` = "\u2283", `312` = "\u2287", `313` = "\u2284",
    `314` = "\u2282", `315` = "\u2286", `316` = "\u2208", `317` = "\u2209",
    # Row D (0320-0337): angle, gradient, logic, double arrows.
    `320` = "\u2220", `321` = "\u2207", `322` = "\u00ae", `323` = "\u00a9",
    `324` = "\u2122", `325` = "\u220f", `326` = "\u221a", `327` = "\u22c5",
    `330` = "\u00ac", `331` = "\u2227", `332` = "\u2228", `333` = "\u21d4",
    `334` = "\u21d0", `335` = "\u21d1", `336` = "\u21d2", `337` = "\u21d3",
    # Row E (0340-0357): lozenge, angle brackets, integral, big-paren
    # pieces (rendered as their nearest plain-Unicode equivalent, since
    # these are technically font-specific glyph fragments meant to be
    # stacked, not standalone characters).
    `340` = "\u25ca", `341` = "\u3008", `342` = "\u00ae", `343` = "\u00a9",
    `344` = "\u2122", `345` = "\u2211", `346` = "\u239b", `347` = "\u239c",
    `350` = "\u239d", `351` = "\u23a1", `352` = "\u23a2", `353` = "\u23a3",
    `354` = "\u23a7", `355` = "\u23a8", `356` = "\u23a9", `357` = "\u23aa",
    # Row F (0360-0377): apple, angle bracket, integral pieces.
    `361` = "\u3009", `362` = "\u222b", `363` = "\u2320", `364` = "\u23ae",
    `365` = "\u2321", `366` = "\u239e", `367` = "\u239f", `370` = "\u23a0",
    `371` = "\u23a4", `372` = "\u23a5", `373` = "\u23a6", `374` = "\u23ab",
    `375` = "\u23ac", `376` = "\u23ad")

# TDA's Times-Roman re-encoding, reached via "FT" (the default text
# font, active whenever "FS" has not been switched to): the Adobe Glyph
# List names spanvec's table actually uses, mapped to Unicode. spanvec
# itself is extracted fresh from each file's prolog (see
# .ps_spanvec()), not hardcoded here, since only the glyph *names* are
# fixed, standard Adobe names -- what octal code each one sits at is
# whatever this specific file's prolog says.
.PS_AGL <- c(
    NUL = "", HT = "\t", LF = "\n", CR = "\r", DLE = "", DC = "", RS = "",
    US = "",
    Eth = "\u00d0", eth = "\u00f0", Lslash = "\u0141", lslash = "\u0142",
    Scaron = "\u0160", scaron = "\u0161", Yacute = "\u00dd",
    yacute = "\u00fd", Thorn = "\u00de", thorn = "\u00fe",
    Zcaron = "\u017d", zcaron = "\u017e", onehalf = "\u00bd",
    onequarter = "\u00bc", onesuperior = "\u00b9",
    threequarters = "\u00be", threesuperior = "\u00b3",
    twosuperior = "\u00b2", brokenbar = "\u00a6", minus = "\u2212",
    multiply = "\u00d7", plusminus = "\u00b1", divide = "\u00f7",
    degree = "\u00b0", approxequal = "\u2248", greaterequal = "\u2265",
    lessequal = "\u2264", notequal = "\u2260", infinity = "\u221e",
    integral = "\u222b", summation = "\u2211", product = "\u220f",
    radical = "\u221a", logicalnot = "\u00ac", mu = "\u03bc", pi = "\u03c0",
    fraction = "\u2044", florin = "\u0192",
    quoteleft = "\u2018", quoteright = "\u2019",
    quotedblleft = "\u201c", quotedblright = "\u201d",
    quotedblbase = "\u201e", quotesinglbase = "\u201a",
    guillemotleft = "\u00ab", guillemotright = "\u00bb",
    guilsinglleft = "\u2039", guilsinglright = "\u203a",
    dagger = "\u2020", daggerdbl = "\u2021", section = "\u00a7",
    paragraph = "\u00b6", periodcentered = "\u00b7", bullet = "\u2022",
    ellipsis = "\u2026", emdash = "\u2014", endash = "\u2013",
    exclamdown = "\u00a1", questiondown = "\u00bf", cent = "\u00a2",
    currency = "\u00a4", sterling = "\u00a3", yen = "\u00a5",
    copyright = "\u00a9", registered = "\u00ae",
    trademark = "\u2122", ordfeminine = "\u00aa",
    ordmasculine = "\u00ba", perthousand = "\u2030", lozenge = "\u25ca",
    germandbls = "\u00df", nbspace = "\u00a0",
    AE = "\u00c6", ae = "\u00e6", OE = "\u0152", oe = "\u0153",
    Oslash = "\u00d8", oslash = "\u00f8", Aring = "\u00c5",
    aring = "\u00e5", Ntilde = "\u00d1", ntilde = "\u00f1",
    Ccedilla = "\u00c7", ccedilla = "\u00e7", Atilde = "\u00c3",
    atilde = "\u00e3", Otilde = "\u00d5", otilde = "\u00f5",
    Agrave = "\u00c0", agrave = "\u00e0", Egrave = "\u00c8",
    egrave = "\u00e8", Igrave = "\u00cc", igrave = "\u00ec",
    Ograve = "\u00d2", ograve = "\u00f2", Ugrave = "\u00d9",
    ugrave = "\u00f9", Aacute = "\u00c1", aacute = "\u00e1",
    Eacute = "\u00c9", eacute = "\u00e9", Iacute = "\u00cd",
    iacute = "\u00ed", Oacute = "\u00d3", oacute = "\u00f3",
    Uacute = "\u00da", uacute = "\u00fa", Acircumflex = "\u00c2",
    acircumflex = "\u00e2", Ecircumflex = "\u00ca",
    ecircumflex = "\u00ea", Icircumflex = "\u00ce",
    icircumflex = "\u00ee", Ocircumflex = "\u00d4",
    ocircumflex = "\u00f4", Ucircumflex = "\u00db",
    ucircumflex = "\u00fb", Adieresis = "\u00c4", adieresis = "\u00e4",
    Edieresis = "\u00cb", edieresis = "\u00eb", Idieresis = "\u00cf",
    idieresis = "\u00ef", Odieresis = "\u00d6", odieresis = "\u00f6",
    Udieresis = "\u00dc", udieresis = "\u00fc", Ydieresis = "\u0178",
    ydieresis = "\u00ff", dotlessi = "\u0131", circumflex = "\u02c6",
    tilde = "\u02dc", macron = "\u00af", breve = "\u02d8",
    dotaccent = "\u02d9", ring = "\u02da", cedilla = "\u00b8",
    hungarumlaut = "\u02dd", ogonek = "\u02db", caron = "\u02c7",
    dieresis = "\u00a8", acute = "\u00b4", grave = "`",
    Delta = "\u0394", Gamma = "\u0393", Omega = "\u03a9", fi = "\ufb01",
    fl = "\ufb02", apple = "\uf8ff")

# spanvec is present in every TDA-generated file's prolog, but not
# assumed identical across TDA versions, so it is extracted fresh from
# the actual file's text rather than hardcoded. Only the glyph
# *names* (.PS_AGL) are assumed fixed and standard.
.ps_spanvec <- function(txt) {
    i <- grep("^/spanvec", txt)
    if (!length(i))
        return(character(0))
    j <- i[1L]
    win <- j:min(j + 400L, length(txt))
    end_rel <- which(grepl("\\]\\s*def", txt[win]))[1L]
    if (is.na(end_rel))
        return(character(0))
    block <- txt[j:(j + end_rel - 1L)]
    m <- regmatches(block, regexec("^8#([0-7]+)\\s+/([A-Za-z0-9]+)", block))
    m <- Filter(function(x) length(x) == 3L, m)
    if (!length(m))
        return(character(0))
    codes <- vapply(m, function(x) x[2L], character(1))
    names_ <- vapply(m, function(x) x[3L], character(1))
    uni <- unname(.PS_AGL[names_])
    uni[is.na(uni)] <- ""
    v <- uni
    names(v) <- codes
    v
}

# The %#text annotation records where a label goes, but for a y axis label it
# is not the baseline PostScript will use: `adjust` lowers the string by a
# third of lpt, while the annotation records only a fifth of it.  The rest has
# to come off here or the number sits above its tick instead of astride it.
.ps_adjust_y <- function(y, s, st) {
    if (!grepl("adjust", s))
        return(y)
    lpt <- st$lpt %||% st$fontsize %||% 8.5
    y - (lpt / 3 - lpt / 5)
}

# setdash takes a PostScript array of on/off lengths in points. grid's
# own lty only offers a few named keywords or a hex-string pattern
# whose digits are scaled by the current line width, not an absolute
# length -- for TDA's thin ~0.567pt lines that scaling crushes
# every one of its 9 line types down to barely-visible micro-dashes,
# indistinguishable from each other even though their real point
# lengths differ by a factor of 9 (checked: rendered side
# by side, all 9 looked like near-identical fine stippling, not the
# clearly separable dash lengths the manual's reference figure
# shows). Returning the raw on/off lengths here instead of mapping
# them through grid's relative encoding lets the actual line-drawing
# code break each path into real dash segments at its true length in
# points, independent of line width -- see .dash_segments().
.ps_dash <- function(spec) {
    v <- .ps_num(strsplit(gsub("[][]", " ", spec), "\\s+")[[1]])
    v <- v[!is.na(v) & v > 0]
    if (!length(v))
        "solid"
    else {
        if (length(v) == 1L)
            v <- c(v, v)
        v
    }
}

# Break a polyline into its actual dash segments at true, absolute
# lengths (points), alternating on/off along the path's cumulative
# arc length -- not through grid's line-width-relative hex encoding,
# which cannot represent TDA's real dash proportions for a thin line
# (see .ps_dash()'s comment). Returns x/y with NA separating each
# "on" run, which grid.lines() treats as a break between segments.
.dash_segments <- function(x, y, dash) {
    seg_len <- sqrt(diff(x)^2 + diff(y)^2)
    total <- sum(seg_len)
    if (total <= 0)
        return(list(x = x, y = y))
    cyc <- sum(dash)
    n_cyc <- ceiling(total / cyc) + 1L
    marks <- cumsum(rep(dash, n_cyc))
    marks <- marks[marks <= total]
    marks <- c(0, marks, total)
    on <- rep(c(TRUE, FALSE), length.out = length(marks) - 1L)
    cum <- c(0, cumsum(seg_len))
    ox <- numeric(0); oy <- numeric(0)
    for (i in which(on)) {
        a <- marks[i]; b <- marks[i + 1L]
        if (b <= a) next
        j0 <- max(1L, findInterval(a, cum))
        j1 <- min(length(x), findInterval(b, cum) + 1L)
        idx <- j0:j1
        px <- x[idx]; py <- y[idx]
        # trim the first/last point to the exact dash boundary
        if (length(idx) >= 2L) {
            f0 <- (a - cum[j0]) / seg_len[j0]
            f0 <- min(max(f0, 0), 1)
            px[1L] <- x[j0] + f0 * (x[j0 + 1L] - x[j0])
            py[1L] <- y[j0] + f0 * (y[j0 + 1L] - y[j0])
            jl <- idx[length(idx) - 1L]
            f1 <- (b - cum[jl]) / seg_len[jl]
            f1 <- min(max(f1, 0), 1)
            px[length(px)] <- x[jl] + f1 * (x[jl + 1L] - x[jl])
            py[length(py)] <- y[jl] + f1 * (y[jl + 1L] - y[jl])
        }
        if (length(ox)) { ox <- c(ox, NA); oy <- c(oy, NA) }
        ox <- c(ox, px); oy <- c(oy, py)
    }
    list(x = ox, y = oy)
}

#' Read TDA's PostScript output
#'
#' Parses the PostScript a plot command writes into a table of drawing
#' operations: line segments, polygons, symbols and text, with the graphics
#' state that applied to each.
#'
#' @param file a \code{.ps} file written by TDA, or a run whose directory
#'   contains one.
#' @param which when \code{file} is a run, the name of the PostScript file.
#' @return A list with \code{ops}, the drawing operations, and \code{bbox},
#'   the bounding box TDA declared. Each operation carries the plot command
#'   that produced it, taken from the \code{\%#} annotations TDA writes into
#'   the stream.
#' @family plotting
#' @examples
#' d <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
#' p <- tda_ps(d, xlim = c(0, 6), ylim = c(0, 10))
#' p <- tda_pl_axes(p)
#' p <- tda_pl_lines(p, "x", "y")
#' pdf(NULL)  # plot() runs the session and draws it; capture the run
#' p <- plot(p)
#' dev.off()
#' ps <- tda_read_ps(p$run, which = p$file)
#' ps$bbox
#' ps$ops[[3]]$op
#' @export
tda_read_ps <- function(file, which = NULL) {
    if (!is.character(file)) {
        d <- file$dir %||% file$run$dir
        if (is.null(d))
            stop("`file` must be a path or a TDA run")
        cand <- if (!is.null(which)) file.path(d, which)
                else list.files(d, "\\.ps$", full.names = TRUE)
        if (!length(cand))
            stop("no PostScript file in ", d)
        file <- cand[1L]
    }
    txt <- readLines(file, warn = FALSE)

    bb <- grep("^%%BoundingBox:", txt, value = TRUE)
    bbox <- if (length(bb))
        .ps_num(strsplit(trimws(sub("^%%BoundingBox:", "", bb[1L])),
                         "\\s+")[[1]]) else rep(NA_real_, 4L)

    st <- .ps_state()
    clipbox <- NULL
    source_cmd <- NA_character_
    pending_text <- NULL
    pending_num <- NULL
    pending_string <- NULL
    pending_val <- numeric(0)
    # dplot() (combining several already-created files into a grid,
    # examples/exam/plot7.cf) writes a global scale factor once, at
    # the very top of the file -- "% changed: global scaling / /sfx
    # 0.43 def / /sfy 0.45 def" -- shrinking each combined panel down
    # to fit its cell, applied later via "sfx sfy scale" (a
    # variable-name scale, resolved by the general ps_vars tracker
    # below, same as xorg/yorg's "xorg yorg translate"). That
    # definition sits before .ps_body()'s %%EndProlog cutoff, so
    # it is not a per-panel coordinate-system variable the main loop
    # ever sees itself -- without seeding it here first, "sfx sfy
    # scale" resolved to nothing at all (ps_vars[["sfx"]] was never
    # set), so no panel was ever actually shrunk to fit, and every
    # one kept its full, original size: checked, a
    # real four-panel file's two panels sharing a row had raw
    # x-extents 283.5pt wide each (each panel's un-scaled 90mm)
    # that overlapped one another by more than half, not the roughly
    # 122pt each a 0.43x shrink to fit a 283.5pt-wide row would give.
    # Scoped to before %%EndProlog specifically (not the whole file)
    # so a coincidentally-named variable defined naturally later, as
    # part of a real panel's coordinate setup, is never
    # overwritten by a stale, page-level value from above it.
    # spanvec (TDA's Times-Roman re-encoding for special characters,
    # see .ps_spanvec()'s comment) is extracted once, up front, from
    # the file's full text -- it lives in the prolog, before the
    # %%EndProlog cutoff .ps_body() applies below, so this has to happen
    # before that, not from the stripped body.
    spanvec_table <- .ps_spanvec(txt)
    ps_vars <- local({
        end <- grep("%%EndProlog", txt)
        head_txt <- if (length(end)) txt[seq_len(end[1L])] else character()
        defs <- grep("^/[A-Za-z_][A-Za-z0-9_]*\\s+[-0-9.eE+]+\\s+def\\b",
                    head_txt, value = TRUE)
        v <- list()
        for (d in defs) {
            nm <- sub("^/([A-Za-z_][A-Za-z0-9_]*).*", "\\1", d)
            val <- .ps_num(sub("^/[A-Za-z_][A-Za-z0-9_]*\\s+([-0-9.eE+]+).*",
                               "\\1", d))
            if (!is.na(val))
                v[[nm]] <- val
        }
        v
    })
    body <- .ps_body(txt)
    # TDA annotates its output with "%#N: command", which says which plot
    # command produced what follows.  Those are kept: they are how a drawing
    # operation is attributed back to its source, which is what a per-command
    # grid conversion needs.
    body <- body[nzchar(trimws(body))]
    body <- body[!grepl("^%", body) | grepl("^%#", body)]

    pending_sym <- NULL
    pending_sym_type <- NULL
    stack <- list()
    path_stack <- list()
    x <- y <- 0
    path <- list()
    path_filled <- FALSE
    ops <- list()
    # "%#Parameter: ..." is written once per coordinate system TDA sets
    # up -- once for an ordinary, single-panel plot, but once *per
    # panel* for a multi-panel layout (psetup() called more than once,
    # or dplot() combining several already-created files into a grid).
    # Tagging every op with which one was most recently seen is what
    # lets the render step lay panels out in their own, separate
    # sub-regions later (see tda_plot_ps()'s comment) instead of
    # treating a whole multi-panel file as one shared coordinate
    # space, which is wrong the moment more than one is on the page:
    # checked against Ghostscript's rendering of a real
    # 2x2 combined file, an evenly-proportioned grid of four panels,
    # not the wildly uneven layout untagged rendering produced. Not
    # "% define current coordinate system" itself, which carries the
    # actual xorg/yorg/pxlen/pylen values in a comment but does not
    # survive the general "%..." (unless "%#...") filter above, and
    # would be stale after dplot() repositions a panel besides (that
    # real value is only in "/xorg N def % changed", already handled
    # correctly by the general ps_vars tracker elsewhere).
    panel <- 0L
    add <- function(...) ops[[length(ops) + 1L]] <<- c(list(...), list(panel = panel))

    flush <- function(kind, keep_path = FALSE) {
        # A lone moveto draws nothing; emitting it would leave zero-length
        # segments in the output.
        if (length(path) < 2L) {
            if (!keep_path)
                path <<- list()
            return()
        }
        p <- do.call(rbind, path)
        add(op = kind, x = p[, 1L], y = p[, 2L], clip = st$clip,
            lwd = st$lwd, lty = st$lty,
            col = st$col, fontsize = st$fontsize, command = source_cmd)
        if (!keep_path)
            path <<- list()
        path_filled <<- keep_path
    }

    # 94% of a spatial plot's PostScript is "<x> <y> l" -- 49003 of 52305
    # lines for a 20000-polygon coastline -- and every one of them walks
    # the whole grepl() cascade below, matching none of it, before
    # reaching the coordinate handler at the bottom.  That cascade was
    # 58% of tda_read_ps()'s total time (Rprof), and it is the reason a
    # large map took minutes rather than seconds.
    #
    # Classified ONCE here, vectorised, instead of ~40 regex calls per
    # line.  The fast path below handles only the plain lineto -- no
    # comment suffix, nothing else on the line -- and everything else
    # falls through to the original code unchanged.
    # options(tdaR.ps_fastpath = FALSE) turns the short-cut off and runs
    # every line through the original cascade.  It exists so the two can
    # be compared: a fast path is only correct if it is a pure
    # short-cut, and the test suite checks that both produce identical
    # ops for every reference .ps in the tree.
    fast_l <- if (isFALSE(getOption("tdaR.ps_fastpath", TRUE)))
        rep(FALSE, length(body))
    else
        grepl("^\\s*[-0-9.eE+]+\\s+[-0-9.eE+]+\\s+l\\s*$", body)
    # The classification above is one vectorised call, but PARSING was
    # still per line: strsplit() plus two as.numeric() for each of the
    # 34820 lineto lines of a world map, which Rprof put at about a
    # fifth of the total.  scan() reads every fast line's two numbers in
    # a single pass instead, and the loop just indexes the result.
    # Kept in step with fast_l by construction: fast_xy has one row per
    # TRUE, in order, and fast_pos maps a body index to that row.
    # A plain "<x> <y> m|rm|rl" matches NONE of the forty-odd patterns in
    # the cascade below -- checked one by one: the only pattern that
    # could, "\\bm$", sits inside the "^stroke" branch and so is never
    # reached by a line starting with a digit.  Those lines can therefore
    # skip straight to the tokeniser, which is where they were always
    # handled, without running any of it.  The handler itself is
    # untouched: its flush/pending-text/fill bookkeeping is real state
    # and duplicating it in a short-cut would be a second copy to keep
    # in step.
    skip_cascade <- grepl("^\\s*[-0-9.eE+]+\\s+[-0-9.eE+]+\\s+(m|rm|rl)\\s*$",
                          body)
    fast_xy <- NULL
    fast_pos <- integer(0)
    if (any(fast_l)) {
        fi <- which(fast_l)
        fast_pos <- integer(length(body))
        fast_pos[fi] <- seq_along(fi)
        # Reading the trailing "l" as a third, character field lets
        # scan() do the whole line, instead of a sub() pass over all
        # 34820 strings first just to remove it.
        v <- suppressWarnings(tryCatch(
            scan(text = body[fi],
                 what = list(double(), double(), character()),
                 quiet = TRUE, na.strings = c("NA", "")),
            error = function(e) NULL))
        # scan() gives up at the first unreadable token, so fall back to
        # the per-line parse for the whole file if it did not return
        # exactly one value per line rather than silently misalign.
        if (!is.null(v) && length(v[[1L]]) == length(fi) &&
            length(v[[2L]]) == length(fi))
            fast_xy <- cbind(v[[1L]], v[[2L]])
        else
            fast_l <- rep(FALSE, length(body))
    }
    body_i <- 0L
    for (ln in body) {
        body_i <- body_i + 1L
        if (fast_l[body_i]) {
            a2 <- fast_xy[fast_pos[body_i], 1L]
            b2 <- fast_xy[fast_pos[body_i], 2L]
            if (anyNA(c(a2, b2)))
                next
            # .ps_map() inlined: it is two multiply-adds, and calling it
            # once per lineto was a quarter of the total time (Rprof) in
            # function-call overhead alone.
            cm <- st$ctm
            p2 <- c(cm[1L] * a2 + cm[3L] * b2 + cm[5L],
                    cm[2L] * a2 + cm[4L] * b2 + cm[6L])
            if (!length(path))
                path <- list(c(x, y))
            path[[length(path) + 1L]] <- p2
            x <- p2[1L]; y <- p2[2L]
            next
        }
        s <- trimws(ln)
        if (!skip_cascade[body_i]) {
        # "stroke" and "closepath" are about a sixth of the lines a map
        # emits (1000 and 993 of 38106 for a world coastline) and they
        # used to sit thirty-odd branches down the cascade, so each ran
        # thirty regexes to reach its own.  Hoisted to the front: no
        # branch that was above them can match a line beginning
        # "stroke" or "closepath" -- checked pattern by pattern -- and
        # the blocks are unchanged, moved verbatim.
        #
        # Unlike the lineto and moveto short-cuts this reorders live
        # control flow, so the fast/slow test is NOT enough on its own:
        # a reorder changes both paths and they would agree while both
        # being wrong.  It was checked by reading every reference .ps in
        # the tree with the pre-reorder build and again after, and
        # requiring identical ops.
        if (grepl("^closepath", s)) {
            if (length(path))
                path[[length(path) + 1L]] <- path[[1L]]
            next
        }
        if (grepl("^stroke", s)) {
            flush("lines")
            # "stroke m" is TDA's arrowhead idiom: "currentpoint /
            # stroke m" strokes the line just drawn, then immediately
            # starts a new path at the same point that line ended on
            # (PostScript's currentpoint, pushed before the stroke,
            # popped by this m) -- which is already exactly what x/y
            # hold here, tracked from the last m/l/rm/rl regardless, so
            # no separate operand-stack simulation is needed for this
            # specific idiom, just recognising the trailing m at all
            # (bare "m" with no coordinates of its own on the same line
            # never matched the general "n-2/n-1/n" tokeniser below,
            # so the moveto -- and the whole arrowhead path that
            # depends on starting from it -- was silently dropped).
            if (grepl("\\bm$", s))
                path <- list(c(x, y))
            next
        }

        if (grepl("^%#", s)) {
            # "%#Parameter: ..." marks a new coordinate system -- see the
            # panel/add() comment above for why this, not "%#N: command"
            # (source_cmd, below) or "% define current coordinate
            # system" (stripped before the main loop even sees it,
            # unless "%#..."), is what a multi-panel layout's
            # panel boundary is tracked from.
            if (grepl("^%#Parameter:", s))
                panel <- panel + 1L
            # Only a "%#N: command" annotation ends a path.  "%#stroke" and
            # "%#text:" are notes about what follows, and flushing on them
            # threw away the tick marks before their real stroke arrived.
            if (!grepl("^%#(stroke|text:)", s))
                flush("lines")
            # A "%#text:" annotation carries the position the string is to be
            # drawn at.  It has to: TDA emits "0 0 m" before the show, so the
            # current point is the origin and the coordinates exist nowhere
            # else in the stream.  Without this every label lands in the
            # bottom left corner.
            # Only "%#N: command" names a command.  "%#stroke" is a note
            # about the path that follows, and treating it as a command name
            # overwrote the real one, so most operations were attributed to
            # "stroke" instead of to the axis or series that drew them.
            # "%#symbol: n x y size" carries the position the symbol is
            # drawn at, the way %#text: does.  The symbol line itself is
            # "x y circle invers x y m cross" -- the coordinates are inline,
            # not the current point, so a parser that uses the current point
            # puts every marker of a scatterplot on top of the last polyline
            # vertex.
            if (grepl("^%#symbol:", s)) {
                v <- .ps_num(strsplit(trimws(sub("^%#symbol:", "", s)),
                                      "\\s+")[[1]])
                v <- v[!is.na(v)]
                if (length(v) >= 3L) {
                    # The x,y TDA writes into this comment are the
                    # coordinates it called ps_sym() with, in whatever
                    # coordinate system was current there -- for a plain
                    # symbol plot that is already the final device
                    # position (no transform in effect beyond the
                    # baseline one), but a sunflower plot's single-
                    # point glyph (scplot2a(), n==1) wraps its own
                    # "X Y translate" around a symbol drawn at the
                    # local (0,0), so treating this comment's value
                    # as an absolute, final position is only correct
                    # when nothing has been translated first. Mapping
                    # it through the current CTM handles both: a no-op
                    # for the plain case (checked: existing
                    # symbol-plot tests are unaffected), and the actual
                    # translated offset for the sunflower case, where
                    # every single-point cell previously collapsed onto
                    # the same, wrong (0,0)-mapped position instead of
                    # its grid cell, rendering as missing dots.
                    mapped <- .ps_map(st$ctm, v[2L], v[3L])
                    pending_sym <- mapped
                    pending_sym_type <- v[1L]
                }
                next
            }
            if (grepl("^%#text:", s)) {
                v <- .ps_num(strsplit(trimws(sub("^%#text:", "", s)),
                                      "\\s+")[[1]])
                v <- v[!is.na(v)]
                if (length(v) >= 2L) {
                    # x, y, and the rotation (5th field: "x y fontsize just
                    # rotation label", confirmed by comparing this line's
                    # own numbers against known rotate= calls) -- when
                    # present, grid.text() gets told to actually rotate the
                    # label, which it silently never did before this.
                    #
                    # The x,y here are whatever TDA called ps_text() with,
                    # in whatever coordinate system was current there --
                    # the same "already final" assumption %#symbol: made
                    # (see that comment), and wrong for the identical
                    # reason: a multi-panel layout's per-panel text
                    # ("gsave / x y translate / 0 0 m / show / grestore")
                    # writes this comment's x,y as the *local*
                    # translate offset, not an absolute page position.
                    # Checked against a real two-panel
                    # replication (examples/exam/plot6.cf): both panels'
                    # own labels came back at the identical, un-offset x,
                    # which is only possible if each panel's physical
                    # origin was never actually applied to it. Mapping
                    # through the current CTM is a no-op for the ordinary,
                    # single-panel case (confirmed: existing text-position
                    # tests are unaffected) and gives each panel's
                    # label its own, distinct position otherwise.
                    mapped <- .ps_map(st$ctm, v[1L], v[2L])
                    pending_text <- c(mapped, if (length(v) >= 5L) v[5L]
                                      else 0)
                }
            } else if (grepl("^%#\\d+:", s)) {
                source_cmd <- trimws(sub("^%#\\d+:\\s*", "", s))
            }
            next
        }
        # TDA defines a name for a numeric value it is about to use --
        # "/xorg 100 def" then later "xorg yorg translate" -- for
        # several different purposes (fsiz/ssiz/lpt below have their
        # own specific handling already; this is everything else, most
        # importantly the physical page origin of each panel in a
        # multi-panel layout, xorg/yorg, and the rotation/scale
        # factors rf/sfx/sfy that go with them). Recorded here as a
        # simple name-to-value lookup so a later "translate"/"rotate"/
        # "scale" whose own arguments are these names rather than
        # literal numbers can resolve them -- without this, such a
        # line's .ps_num() calls all returned NA, so the transform
        # silently did nothing at all. For a single-panel plot (the
        # overwhelming majority) that is invisible: the whole drawing
        # is uniformly missing the same one offset, and this renderer
        # scales/centres on the content's extent regardless of
        # its absolute page position. For a multi-panel layout
        # (psetup() called more than once at a different psorg=,
        # examples/exam/plot6.cf) every panel's content ends up at
        # the same, un-offset position instead of side by side --
        # checked, not by inspection: two panels' own
        # labels came back at the identical x, which is only possible
        # if the offset separating them was never applied.
        if (grepl("^/[A-Za-z_][A-Za-z0-9_]*\\s+[-0-9.eE+]+\\s+def\\b", s)) {
            vname <- sub("^/([A-Za-z_][A-Za-z0-9_]*).*", "\\1", s)
            vval <- .ps_num(sub("^/[A-Za-z_][A-Za-z0-9_]*\\s+([-0-9.eE+]+).*",
                                "\\1", s))
            if (!is.na(vval))
                ps_vars[[vname]] <- vval
        }
        # Several operators can share a line -- "/fsiz 8.50 def FT (x) show"
        # is one line -- so the state settings are stripped off and the rest
        # of the line carries on being examined rather than being skipped.
        if (grepl("^/fsiz", s)) {
            # Text size is drawn through the same transform everything
            # else here is: PostScript's "scale" affects how large
            # a shown glyph actually appears, the same as it affects
            # line coordinates, but the raw "/fsiz N def" value itself
            # never changes -- it is the font *before* that transform
            # is applied, not the true, on-page size. TDA's
            # dplot() (combining several already-created files into a
            # grid, examples/exam/plot7.cf) shrinks each panel with
            # exactly this kind of scale (sfx/sfy, often well under 1),
            # and reading fontsize as the raw, un-scaled "/fsiz" value
            # left every panel's text at its original, full size
            # relative to geometry that actually got smaller -- text
            # many times too large for its own panel, confirmed
            # directly: a real combined file's titles overlapped
            # across the whole image rather than sitting inside their
            # own, much smaller panels. sqrt(a^2+b^2) is the CTM's
            # linear part's magnitude (its rotation/scale, not
            # translation) -- the standard way a uniform-ish transform's
            # own effective scale factor is read off, and a no-op for
            # the ordinary, unscaled case (ctm's a=1,b=0 there).
            st$fontsize <- .ps_num(sub("^/fsiz\\s+([0-9.]+).*", "\\1", s)) *
                sqrt(st$ctm[1L]^2 + st$ctm[2L]^2)
            # "FT"/"FS" here is which font the show that follows draws
            # through -- Times-Roman (spanvec-encoded) or the standard
            # Symbol font -- and which of .PS_AGL/.PS_SYMBOL_FONT an
            # octal escape in that show's string decodes through
            # (see .ps_unescape()'s comment for why this matters
            # at all). Defaults to "FT", TDA's default text font,
            # for any /fsiz line that does not name one explicitly.
            fm <- regmatches(s, regexpr("(?<=def\\s)(FT|FS)(?=[\\s(])", s,
                                        perl = TRUE))
            if (length(fm) && nzchar(fm))
                st$font <- fm
            s <- trimws(sub("^/fsiz\\s+[0-9.]+\\s*def\\s*(FT|FS)?", "", s))
            if (!nzchar(s))
                next
        }
        if (grepl("^/ssiz", s)) {
            st$ssiz <- .ps_num(sub("^/ssiz\\s+([0-9.]+).*", "\\1", s))
            next
        }
        # lpt is the size the axis routines use for their labels, which is
        # what `adjust` divides to lower a y axis label onto its tick.
        if (grepl("^/lpt", s)) {
            st$lpt <- .ps_num(sub("^/lpt\\s+([0-9.]+).*", "\\1", s))
            next
        }
        if (grepl("setdash", s)) {
            st$lty <- .ps_dash(sub("\\s*0?\\s*setdash.*", "", s))
            next
        }
        if (grepl("setlinewidth", s)) {
            st$lwd <- max(0.1, .ps_num(sub("\\s*setlinewidth.*", "", s)))
            next
        }
        if (grepl("setgray", s)) {
            g <- .ps_num(sub("\\s*setgray.*", "", s))
            if (!is.na(g))
                st$col <- grDevices::gray(min(1, max(0, g)))
            # "0.8000 setgray fill grestore" is one line: consuming it here
            # and moving on would lose the fill, which is how shaded nodes
            # came back unshaded.
            s <- trimws(sub("^[-0-9.eE+]*\\s*setgray", "", s))
            if (!nzchar(s))
                next
        }
        if (grepl("\\b(translate|rotate|scale)\\b", s)) {
            for (op in c("translate", "rotate", "scale")) {
                if (!grepl(paste0("\\b", op, "\\b"), s))
                    next
                toks <- strsplit(trimws(sub(paste0("\\s*", op, ".*"),
                                            "", s)), "\\s+")[[1]]
                # A token here is either a literal number ("100 600
                # translate") or a name defined earlier with "/name
                # value def" ("xorg yorg translate", TDA's usual
                # form -- see the general variable tracker above for
                # why resolving it matters).
                v <- vapply(toks, function(tk) {
                    n <- .ps_num(tk)
                    if (!is.na(n)) n
                    else if (!is.null(ps_vars[[tk]])) ps_vars[[tk]]
                    else NA_real_
                }, numeric(1))
                v <- v[!is.na(v)]
                # "rotate" alone, with no angle on its own line, is TDA's
                # own arrowhead idiom -- the angle came from an earlier
                # "dx dy atan" instead (see the pending_val handling
                # above), left on the same one-value stack "rotate"
                # itself would pop from in real PostScript. Without this,
                # every arrowhead used the same, unrotated default
                # orientation regardless of which way its line
                # actually pointed -- checked against a real
                # eight-arrow test file (examples/exam/plot13.ps) and an
                # independent Ghostscript rendering of it.
                if (op == "rotate" && !length(v) && length(pending_val)) {
                    v <- pending_val[length(pending_val)]
                    pending_val <- pending_val[-length(pending_val)]
                }
                if (op == "translate" && length(v) >= 2L)
                    st$ctm <- .ps_concat(st$ctm,
                                         c(1, 0, 0, 1, v[1L], v[2L]))
                else if (op == "rotate" && length(v) >= 1L) {
                    th <- v[length(v)] * pi / 180
                    st$ctm <- .ps_concat(st$ctm,
                                         c(cos(th), sin(th),
                                           -sin(th), cos(th), 0, 0))
                } else if (op == "scale" && length(v) >= 2L)
                    st$ctm <- .ps_concat(st$ctm,
                                         c(v[1L], 0, 0, v[2L], 0, 0))
            }
            # The current point is in user space and does not move with the
            # transform; only its mapping changes.
            next
        }
        if (grepl("^gsave", s)) {
            stack[[length(stack) + 1L]] <- st
            # PostScript's gsave/grestore save and restore the
            # current path too, not just colour/CTM/etc -- confirmed
            # directly against a real file's output (examples/exam/
            # plot19.cf's plotsp): "[draw curve] / gsave / [extend
            # path down to the x axis, closing it] / fill grestore /
            # stroke" fills the full, closed polygon, then grestore
            # undoes the path extension, so the stroke that follows
            # only outlines the curve itself, not the polygon's
            # bottom and sides. path was never included in what gsave/
            # grestore save and restore here, so that stroke used
            # whatever path had accumulated by that point regardless --
            # the full, closed polygon, drawing a border TDA's
            # output never has.
            path_stack[[length(path_stack) + 1L]] <- path
            # "gsave 0.95 setgray" is one line -- TDA's plot-frame
            # background fill (plframe's gs=) does exactly this, and
            # the general setgray handling above only ever matches a
            # line that *starts* with a number or with "gsave" alone;
            # this trailing setgray, right after "gsave", never reached
            # it, so the fill that followed used whatever colour was
            # active before the gsave (black, checked: a
            # real background rectangle came back solid black instead
            # of the light grey plframe(gs=0.95) actually asked for --
            # about as visible a wrong colour as this parser can
            # produce) instead of the one this line just set.
            s <- trimws(sub("^gsave", "", s))
            if (grepl("setgray", s)) {
                g <- .ps_num(sub("\\s*setgray.*", "", s))
                if (!is.na(g))
                    st$col <- grDevices::gray(min(1, max(0, g)))
            }
            next
        }
        if (grepl("^grestore", s)) {
            if (length(stack)) {
                st <- stack[[length(stack)]]
                stack[[length(stack)]] <- NULL
            }
            if (length(path_stack)) {
                path <- path_stack[[length(path_stack)]]
                path_stack[[length(path_stack)]] <- NULL
            }
            next
        }
        # clip consumes the current path to set the clipping region and
        # newpath discards it -- neither draws anything.  Skipping them while
        # leaving the path in place meant TDA's clip rectangle was stroked at
        # the next flush, putting a frame round every plot that nobody asked
        # for.
        if (grepl("^clip|^newpath", s)) {
            # A clip path bounds everything drawn after it.  Remembering it
            # is what lets the viewport be the plot area TDA set up rather
            # than the extent of whatever was drawn: a series running past
            # the axis range is clipped by TDA in the PostScript, and without
            # this the replayed picture rescales to include it.
            if (grepl("^clip", s) && length(path) > 2L) {
                p <- do.call(rbind, path)
                st$clip <- c(range(p[, 1L]), range(p[, 2L]))
                # A multi-panel layout (psetup() called more than once, or
                # dplot() combining several already-created files into a
                # grid) sets a *different* clip region per panel -- keeping
                # only the first one seen (an earlier version of this)
                # forced every later panel's content to be clamped down
                # to the first panel's small area, discarding the rest
                # of the grid visually even though the data for it was
                # parsed and present: checked against a real
                # 2x2 combined file (examples/exam/plot7.cf) where only one
                # of four panels ever became visible. The union (min/max
                # across every clip region seen) is the correct viewport
                # for the whole file; a single-panel file has only one
                # clip region regardless, so this is unchanged there.
                clipbox <- if (is.null(clipbox)) st$clip
                    else c(min(clipbox[1L], st$clip[1L]),
                           max(clipbox[2L], st$clip[2L]),
                           min(clipbox[3L], st$clip[3L]),
                           max(clipbox[4L], st$clip[4L]))
            }
            path <- list()
            next
        }
        if (grepl("^showpage|^sclear", s)) next
        if (grepl("^fill", s)) {
            # A fill with no path belongs to the circle just emitted: TDA
            # draws a shaded node as "x y r 0 360 arc" then "gs setgray
            # fill".  Without this the node keeps the outline colour.
            if (length(path) < 2L && length(ops) &&
                identical(ops[[length(ops)]]$op, "circle")) {
                ops[[length(ops)]]$fill <- st$col
                path <- list()
            } else {
                # Real PostScript does not clear the current path on fill
                # (only newpath, or implicitly starting a fresh moveto,
                # does) -- a fill immediately followed by a stroke of the
                # same path, with no newpath between them, is a standard
                # "shade the region, then outline it" idiom, and TDA's
                # contour output (plotc) relies on exactly this: each band
                # is "gsave / setgray fill grestore / stroke / grestore".
                # Clearing the path here left that stroke with nothing to
                # draw at all -- checked (contour lines were
                # silently missing, not merely a wrong colour) -- so the
                # path survives a fill and is only cleared once something
                # else (stroke, a new moveto, clip/newpath) actually
                # consumes it.
                flush("polygon", keep_path = TRUE)
            }
            # "fill grestore" often shares one line with the setgray that
            # set the fill colour (the whole line is "0.9000 setgray fill
            # grestore"), and the general gsave/grestore handling above
            # only ever matches a line that *starts* with "grestore" --
            # this trailing one, embedded after "fill", never reached it,
            # so the colour push/pop around the fill never actually
            # popped: the stroke that followed kept the fill's grey
            # instead of the colour active before the fill, confirmed
            # directly (a real run's stroke came back the identical
            # colour as its fill, not the black TDA's line-drawing
            # colour otherwise always is). Popped here explicitly instead.
            s <- trimws(sub("^fill", "", s))
            if (grepl("^grestore", s) && length(stack)) {
                st <- stack[[length(stack)]]
                stack[[length(stack)]] <- NULL
                # This "grestore" is embedded after "fill" on the same
                # line ("0.9000 setgray fill grestore") rather than
                # starting its line, so it never reaches the main,
                # standalone "^grestore" handler above -- which is
                # also where path is popped back (see that handler's
                # own comment on why). Missing that here left
                # path_stack permanently one entry too deep after
                # every "fill grestore" line, so the *next* standalone
                # grestore restored the wrong saved path -- confirmed
                # directly: examples/exam/plot21.cf's contour
                # bands, which use exactly this combined-line form,
                # lost their outline stroke entirely once path
                # was made gsave/grestore-aware at all.
                if (length(path_stack)) {
                    path <- path_stack[[length(path_stack)]]
                    path_stack[[length(path_stack)]] <- NULL
                }
            }
            next
        }

        # An arrowhead ("currentpoint / stroke m / dx dy / atan / rotate")
        # pushes the direction vector as two bare numbers on their own
        # line, for atan to consume -- not to be confused with the
        # single-bare-number case just below (a bare axis-label value).
        if (grepl("^[-0-9.eE+]+\\s+[-0-9.eE+]+\\s*$", s)) {
            pending_val <- c(pending_val, .ps_num(strsplit(trimws(s),
                                                           "\\s+")[[1]]))
            next
        }
        if (grepl("^atan\\b", s)) {
            if (length(pending_val) >= 2L) {
                n <- length(pending_val)
                # PostScript's "num1 num2 atan angle" is atan2(num1,
                # num2) -- num1 pushed first/deeper on the stack, num2
                # second/on top, matching left-to-right reading order on
                # TDA's "dx dy atan" line. Swapped here at first
                # (atan2(num2, num1) instead), which for a line at
                # exactly 45 degrees gives back the same angle -- a fixed
                # point of swapping the two arguments -- so those
                # arrowheads looked right while every other angle came
                # out rotated 90 degrees off, checked against
                # the real eight-arrow test file once actually looked at
                # closely rather than glanced over as "now correct."
                ang <- atan2(pending_val[n - 1L], pending_val[n]) * 180 / pi
                if (ang < 0) ang <- ang + 360
                pending_val <- c(pending_val[seq_len(n - 2L)], ang)
            }
            next
        }
        # Axis labels are bare numbers, not parenthesised strings: pl_axis
        # prints the value with the current format and then calls center and
        # show, so the number sits on its own line.
        if (grepl("^[-0-9.eE+]+\\s*$", s)) {
            pending_num <- trimws(s)
            next
        }
        if (!is.null(pending_num) &&
            grepl("^(center|adjust|aright|show)\\b", s)) {
            just <- if (grepl("center", s)) "centre"
                    else if (grepl("adjust|aright", s)) "right" else "left"
            # The %#text annotation carries the position the label is finally
            # drawn at, after the rm that centres or right-adjusts it.  The
            # current point is where the tick is, which for a y axis label
            # sits a third of the font size too high.
            tx <- if (!is.null(pending_text)) pending_text[1L] else x
            ty <- if (!is.null(pending_text)) pending_text[2L] else y
            add(op = "text", x = tx, y = .ps_adjust_y(ty, s, st),
                label = pending_num, just = just,
                rot = if (!is.null(pending_text)) pending_text[3L] else 0,
                lwd = st$lwd, lty = st$lty, col = st$col,
                fontsize = st$fontsize, command = source_cmd)
            if (grepl("show", s)) {
                pending_num <- NULL
    pending_sym <- NULL
                pending_text <- NULL
    pending_num <- NULL
    pending_sym <- NULL
            }
            next
        }

        # "(string)" with no show/center/adjust/aright keyword on the same
        # line, and the standalone keyword-only line that follows it: TDA
        # sometimes splits these across two lines rather than one --
        # checked, examples/exam/plots.cf's font-switching
        # text ("/fsiz 8.50 def FS (\141)" on its own line, "show" on the
        # next) -- and the combined "(...)  show" pattern just below this
        # one only ever matches when both are on the same line, so this
        # whole segment (and, for a label built from several font
        # switches in a row, every earlier segment before the last one)
        # was silently dropped rather than merely mis-decoded.
        if (grepl("^\\([^()]*\\)\\s*$", s)) {
            pending_string <- list(
                text = sub("^\\((.*)\\)\\s*$", "\\1", s), font = st$font)
            next
        }
        if (!is.null(pending_string) &&
            grepl("^(show|center|adjust|aright)\\b", s)) {
            lab <- .ps_unescape(pending_string$text, font = pending_string$font,
                                spanvec_table = spanvec_table)
            just <- if (grepl("center", s)) "centre"
                    else if (grepl("adjust|aright", s)) "right" else "left"
            lx <- if (!is.null(pending_text)) pending_text[1L] else x
            ly <- if (!is.null(pending_text)) pending_text[2L] else y
            add(op = "text", x = lx, y = .ps_adjust_y(ly, s, st),
                label = lab, just = just,
                rot = if (!is.null(pending_text)) pending_text[3L] else 0,
                lwd = st$lwd, lty = st$lty, col = st$col,
                fontsize = st$fontsize, command = source_cmd)
            pending_string <- NULL
            if (grepl("show", s)) {
                pending_text <- NULL
                pending_num <- NULL
                pending_sym <- NULL
            }
            next
        }

        # text: (string) show, optionally through center or adjust -- or,
        # with TDA's wf=1 (a white box behind the text), (string)
        # sclear show: an extra token between the string and the keyword
        # that a plain "\\)\\s*keyword" match does not allow for, which
        # silently dropped every wf=1 label until this was loosened to
        # allow for one extra token in between. That extra-token group has
        # to be optional and non-greedy, or it swallows the keyword itself
        # when there is no extra token (plain "(string) center"), which a
        # first attempt at this fix did -- caught by testing that case
        # again after the fix still working.
        if (grepl("\\)\\s*(?:\\S+\\s+)?(show|center|adjust|aright)", s,
                  perl = TRUE)) {
            lab <- .ps_unescape(
                sub(paste0("^.*?\\((.*)\\)\\s*(?:\\S+\\s+)?",
                          "(show|center|adjust|aright).*$"),
                    "\\1", s, perl = TRUE),
                font = st$font, spanvec_table = spanvec_table)
            just <- if (grepl("center", s)) "centre"
                    else if (grepl("adjust|aright", s)) "right" else "left"
            # The %#text annotation gives where the label is finally drawn,
            # after the rm that centres or right-adjusts it.  The current
            # point is the tick, which for a y axis label is a third of the
            # font size too high.
            lx <- if (!is.null(pending_text)) pending_text[1L] else x
            ly <- if (!is.null(pending_text)) pending_text[2L] else y
            add(op = "text", x = lx, y = .ps_adjust_y(ly, s, st),
                label = lab, just = just,
                rot = if (!is.null(pending_text)) pending_text[3L] else 0,
                lwd = st$lwd, lty = st$lty, col = st$col,
                fontsize = st$fontsize, command = source_cmd)
            pending_text <- NULL
    pending_num <- NULL
    pending_sym <- NULL
            next
        }

        # ploto draws a circle as a raw arc rather than through the circle
        # procedure: "x y r a1 a2 arc".  Without this the commonest symbol in
        # TDA's output is silently dropped.
        if (grepl("\\barc\\b", s)) {
            v <- .ps_num(strsplit(trimws(sub("\\s*arc.*", "", s)),
                                  "\\s+")[[1]])
            v <- v[!is.na(v)]
            if (length(v) >= 3L) {
                # "x y r a1 a2 arc": a node is 0 to 360, a full circle, but a
                # curved edge is a segment -- 0 to 90 for a quarter.  Drawing
                # every arc as a circle puts a stray ring on the plot wherever
                # an edge is bent.
                a1 <- if (length(v) >= 4L) v[4L] else 0
                a2 <- if (length(v) >= 5L) v[5L] else 360
                cc <- .ps_map(st$ctm, v[1L], v[2L])
                # The radius alone is only correct for a uniformly-scaled
                # circle (a plain marker, always drawn this way). plote's
                # own technique for an ellipse is a unit circle under a
                # non-uniform scale (and possibly a rotation) -- collapsing
                # that to a single radius via the determinant, as this used
                # to, silently rendered every non-circular ellipse as a
                # circle of some averaged size. The linear part of the CTM
                # is kept too, so a true (possibly rotated) ellipse can be
                # drawn when it is not just a uniform scale.
                add(op = "circle", x = cc[1L], y = cc[2L],
                    r = v[3L] * sqrt(abs(st$ctm[1L] * st$ctm[4L] -
                                         st$ctm[2L] * st$ctm[3L])),
                    r0 = v[3L], m = st$ctm[1L:4L],
                    a1 = a1, a2 = a2,
                    lwd = st$lwd, lty = st$lty, col = st$col,
                    fontsize = st$fontsize, command = source_cmd,
                    clip = st$clip)
                x <- v[1L]; y <- v[2L]
                path <- list()
                next
            }
        }

        if (!is.null(pending_sym_type)) {
            sx <- if (!is.null(pending_sym)) pending_sym[1L] else x
            sy <- if (!is.null(pending_sym)) pending_sym[2L] else y
            # sx and sy come from the %#symbol annotation, which is already
            # in device coordinates, as are x and y once a moveto has mapped
            # them.  Mapping either again applies the transform twice.
            # The symbol-drawing line itself (e.g. "circle invers ... m
            # cross") is TDA's PostScript procedure calls, not
            # something meant to be read back token by token -- it is
            # identified and skipped here by having just seen the
            # %#symbol: comment that always immediately precedes it, using
            # the type number TDA itself put there (see .SYMBOL_TABLE)
            # rather than guessing from whichever procedure names happen to
            # appear in the line, which previously picked "circle" from
            # "circle invers ... cross" every time and silently dropped
            # both that the circle is filled white (invers) and that a
            # cross is drawn over it -- on a white background the result
            # was a correctly-positioned but entirely invisible marker.
            add(op = "symbol", x = sx, y = sy,
                type = pending_sym_type,
                symbol = (.SYMBOL_TABLE[[as.character(pending_sym_type)]] %||%
                          list(shape = "circle"))$shape,
                size = st$ssiz,
                lwd = st$lwd, lty = st$lty, col = st$col,
                fontsize = st$fontsize, command = source_cmd,
                clip = st$clip)
            pending_sym <- NULL
            pending_sym_type <- NULL
            next
        }

        }

        tok <- strsplit(s, "\\s+")[[1]]
        n <- length(tok)
        # TDA marks a point it adds purely to close a path to the plot's
        # own edge or corner for filling purposes with a trailing
        # "% ccc"/"% bbb"-style comment. Stripping it is required just to
        # parse the line at all -- without it, "283.50 0 l % ccc" has
        # "ccc" as its final token, not "l", so the line silently
        # matched no known path command and was dropped outright: the
        # path was missing its last point, and grid's implicit
        # "close back to the start" then drew a straight line across the
        # gap instead of the real point TDA actually placed there. The
        # point itself is real and is kept, the same as any other --
        # checked against an independent PDF rendering of a
        # real plotc output (examples/exam/plot21.ps) that every band,
        # including the ones this comment marks, is both filled and
        # outlined in full.
        if (n >= 2L && tok[n - 1L] == "%") {
            tok <- tok[seq_len(n - 2L)]
            n <- length(tok)
        }
        if (n >= 3L && tok[n] %in% c("m", "l", "rm", "rl")) {
            a <- .ps_num(tok[n - 2L])
            b <- .ps_num(tok[n - 1L])
            if (anyNA(c(a, b)))
                next
            # A path is built in device space point by point.  The matrix in
            # force when a point is emitted is the one that maps it: an
            # arrowhead does "x y m" before its rotate, so the start is in the
            # untransformed system and the rl offsets that follow are in the
            # rotated one.  Mapping the finished path with the final matrix
            # sends the start somewhere else entirely.
            if (tok[n] %in% c("m", "l")) {
                p2 <- .ps_map(st$ctm, a, b)
            } else {
                # rm and rl are offsets, so only the linear part applies.
                p2 <- c(x + st$ctm[1L] * a + st$ctm[3L] * b,
                        y + st$ctm[2L] * a + st$ctm[4L] * b)
            }
            if (tok[n] %in% c("m", "rm")) {
                # flush() ends the path; the pending %#text position belongs
                # to the label about to be shown and has to survive the rm
                # that centres or adjusts it.
                keep_text <- pending_text
                # A path already emitted once by fill (kept around only in
                # case a stroke immediately followed it, TDA's "shade
                # then outline" idiom -- see fill's handling above) and
                # never actually stroked is discarded here, not re-emitted:
                # flushing it again would draw the same shape a second time,
                # as a spurious, uncoloured "lines" op nothing in the
                # original PostScript ever asked for.
                if (path_filled) path <- list() else flush("lines")
                pending_text <- keep_text
                path <- list(p2)
                path_filled <- FALSE
            } else {
                if (!length(path))
                    path <- list(c(x, y))
                path[[length(path) + 1L]] <- p2
            }
            x <- p2[1L]; y <- p2[2L]
        }
    }
    flush("lines")

    # Consecutive text ops at the same position are TDA's multi-font
    # text, split into one "show" per font switch (see .ps_unescape()'s
    # own comment for why): "alpha", " = 2.5 (bullet) ", "beta" are three
    # separate text ops here, each already correctly decoded on its own,
    # but meant to be read as the single label "alpha = 2.5 (bullet)
    # beta", not three overlapping ones at the same spot. Merged here,
    # after parsing, rather than tracked through the main loop's
    # state machine, since by this point every segment's label is
    # already final and only needs concatenating.
    ops <- local({
        out <- list()
        i <- 1L
        n <- length(ops)
        while (i <= n) {
            o <- ops[[i]]
            if (identical(o$op, "text")) {
                j <- i
                lab <- o$label
                while (j < n && identical(ops[[j + 1L]]$op, "text") &&
                       isTRUE(all.equal(unname(ops[[j + 1L]]$x), unname(o$x))) &&
                       isTRUE(all.equal(unname(ops[[j + 1L]]$y), unname(o$y))) &&
                       identical(ops[[j + 1L]]$just, o$just)) {
                    j <- j + 1L
                    lab <- paste0(lab, ops[[j]]$label)
                }
                o$label <- lab
                out[[length(out) + 1L]] <- o
                i <- j + 1L
            } else {
                out[[length(out) + 1L]] <- o
                i <- i + 1L
            }
        }
        out
    })

    list(ops = ops, bbox = bbox, clip = clipbox, file = file)
}

#' Draw TDA's PostScript output with grid
#'
#' Replays what \code{\link{tda_read_ps}} parsed onto the current grid device,
#' so a plot TDA produced appears in an R graphics device rather than only in
#' a \code{.ps} file.
#'
#' @param x a path, a TDA run, or the result of \code{\link{tda_read_ps}}.
#' @param newpage start a new page before drawing.
#' @param cex character expansion for the text, relative to the size TDA
#'   asked for. The sizes in the file are PostScript points against a plot of
#'   perhaps 80mm; drawn into a device several times that, they come out
#'   small, so the default scales them with the drawing, times
#'   \code{getOption("tdaR.ps.cex", 1)}.
#' @param ... passed to \code{\link{tda_read_ps}}.
#' @return The parsed operations, invisibly.
#' @family plotting
#' @examples
#' # a PostScript file with something in it: five points on a quadratic.
#' # psetup needs its ranges -- axes alone draw nothing, and tda_plot_ps
#' # then has nothing to render.
#' r <- tda_run(c("nvar(noc = 5, X = case, Y = case * case);",
#'                "psfile = p.ps;",
#'                "psetup(pxlen = 90, pylen = 50, pxa = 0,6, pya = 0,30);",
#'                "plxa(sc = 1); plya(sc = 10); plframe;",
#'                "plot = X,Y;", "psclose;"))
#' p <- tda_plot_ps(file.path(r$dir, "p.ps"))
#' @export
tda_plot_ps <- function(x, newpage = TRUE, cex = NULL, ...) {
    p <- if (is.list(x) && !is.null(x$ops)) x else tda_read_ps(x, ...)
    if (!length(p$ops))
        stop("nothing to draw: the PostScript file has no drawing operations")

    # A circle extends a radius beyond its centre, so the range has to allow
    # for that or one drawn near an edge is clipped in half -- but only for
    # a full circle. A partial arc (a1/a2 not spanning 360 degrees,
    # TDA's technique for a shallow curve: a tiny slice of a huge-radius
    # circle) has a real extent nowhere near that big, and using the full
    # radius for it was inflating the whole plot's bounding box to fit an
    # almost entirely off-screen circle, shrinking everything else into a
    # corner. For those, the extent comes from points actually sampled
    # along the visible arc instead.
    #
    # o$r (R's partial-matching $) silently matched a text op's rot
    # field when no op$r existed, treating a rotation angle as if it were
    # a radius and corrupting the bounding box for any plot with a
    # rotated label in it -- caught by running with
    # options(warnPartialMatchDollar = TRUE) rather than found by
    # inspection. Every access below is guarded by op == "circle" first,
    # and uses [[ instead of $, which never partially matches.
    arc_extent <- function(o) {
        ang <- seq(o[["a1"]], o[["a2"]], length.out = 12L) * pi / 180
        list(x = o[["x"]] + o[["r"]] * cos(ang),
            y = o[["y"]] + o[["r"]] * sin(ang))
    }
    is_circle <- function(o) identical(o[["op"]], "circle")
    is_partial <- function(o) is_circle(o) &&
        !is.null(o[["a1"]]) && !is.null(o[["a2"]]) &&
        abs((o[["a2"]] - o[["a1"]]) %% 360) > 1e-6 &&
        abs((o[["a2"]] - o[["a1"]]) %% 360 - 360) > 1e-6
    rad <- function(o) if (!is_circle(o) || is_partial(o)) 0 else o[["r"]]
    # A text op's o$y is not where it actually ends up drawn: the
    # render step below shifts it down by 0.32 * fontsize to undo
    # PostScript's baseline-vs-centre convention (see that comment for
    # why). Left out of the extent here, that shift is invisible to
    # the range/pad computation that follows -- harmless whenever
    # there is enough other content to give the range real spread, but
    # for a plot whose only content is a single text op (confirmed
    # directly: examples/tests/p3.ps, one label and nothing else) xr/yr
    # collapse to that one point, pad is a few hundredths of a point,
    # and the fixed, points-sized text offset then lands the label
    # entirely outside the viewport's native range -- rendering
    # nothing, not merely a wrong position.
    text_y_off <- function(o)
        if (identical(o[["op"]], "text")) 0.32 * (o[["fontsize"]] %||% 8.5) else 0
    # A label occupies width, and o$x is only its anchor: a y axis tick
    # label is right-adjusted, so it runs leftwards from the tick and the
    # whole of its footprint lies outside the extent o$x alone reports.
    # The viewport is then sized to content that stops at the axis, and
    # the labels -- the one part of a plot guaranteed to sit outside it --
    # are cut off at the device edge. TDA sets text in the same units as
    # everything else it draws, so the footprint is estimable from the
    # label's length and font size; 0.6 is about the average advance
    # of the Helvetica named in TDA's prologue, near enough for a
    # bounding box. A rotated label spends that footprint on the other
    # axis instead.
    text_span <- function(o) {
        if (!identical(o[["op"]], "text")) return(c(0, 0))
        w <- nchar(o[["label"]] %||% "") * 0.6 * (o[["fontsize"]] %||% 8.5)
        switch(o[["just"]] %||% "left",
               right  = c(-w, 0),
               centre = c(-w / 2, w / 2),
               c(0, w))
    }
    upright <- function(o) {
        rot <- o[["rot"]] %||% 0
        abs(((rot + 45) %% 180) - 45) <= 45
    }
    text_x_span <- function(o) if (upright(o)) text_span(o) else c(0, 0)
    text_y_span <- function(o) if (upright(o)) c(0, 0) else text_span(o)
    extent <- function(ops) {
        # A circle with its own recorded clip region (see .draw_arc()'s
        # own comment on clip-drawing for why a circle needs this at
        # all) contributes only the part of itself that clip allows,
        # not its full, unclipped radius -- otherwise a circle far
        # larger than the plot itself (TDA's technique for a
        # near-straight edge or a full-bleed fill, checked:
        # examples/exam/gd18.cf's rd=-3 edges use a similarly
        # oversized radius) inflates the whole viewport to fit content
        # that was never actually visible past the plot's
        # boundary, shrinking everything else into a sliver.
        clip_extent <- function(o) {
            cl <- o[["clip"]]
            xr <- c(o[["x"]] - rad(o), o[["x"]] + rad(o))
            yr <- c(o[["y"]] - rad(o), o[["y"]] + rad(o))
            if (!is.null(cl) && all(is.finite(cl))) {
                xr <- c(max(xr[1L], cl[1L]), min(xr[2L], cl[2L]))
                yr <- c(max(yr[1L], cl[3L]), min(yr[2L], cl[4L]))
            }
            list(x = xr, y = yr)
        }
        xs <- unlist(lapply(ops, function(o) {
            if (is_partial(o)) arc_extent(o)$x
            else if (is_circle(o) && !is.null(o[["clip"]])) clip_extent(o)$x
            else if (identical(o[["op"]], "text")) o[["x"]] + text_x_span(o)
            else c(o[["x"]] - rad(o), o[["x"]] + rad(o))
        }))
        ys <- unlist(lapply(ops, function(o) {
            if (is_partial(o)) arc_extent(o)$y
            else if (is_circle(o) && !is.null(o[["clip"]])) clip_extent(o)$y
            else if (identical(o[["op"]], "text"))
                o[["y"]] + c(0, text_y_off(o)) + text_y_span(o)
            else c(o[["y"]] - rad(o), o[["y"]] + rad(o))
        }))
        list(xr = range(xs, na.rm = TRUE), yr = range(ys, na.rm = TRUE))
    }

    e <- extent(p$ops)
    xr <- e$xr
    yr <- e$yr
    # Where TDA set a clipping region, that is the plot area it intended, so
    # the view is taken from it rather than from the drawn extent -- a series
    # running past the axis range would otherwise shrink everything else.
    if (!is.null(p$clip) && all(is.finite(p$clip))) {
        xr <- range(xr[1L], p$clip[1:2])
        yr <- range(yr[1L], p$clip[3:4])
        xr[2L] <- max(p$clip[2L], min(xr[2L], p$clip[2L]))
        yr[2L] <- max(p$clip[4L], min(yr[2L], p$clip[4L]))
    } else if (!is.null(p$bbox) && length(p$bbox) == 4L &&
               all(is.finite(p$bbox))) {
        # No clip: a graph drawing sets none.  The file's BoundingBox
        # is then TDA's statement of the area it drew into, whitespace
        # included -- gd31.cf declares 0-8 by -1-5 and puts its four nodes
        # in the middle third, so fitting the view to the nodes alone
        # blows them up to fill a frame TDA leaves mostly empty.
        xr <- range(xr, p$bbox[c(1L, 3L)])
        yr <- range(yr, p$bbox[c(2L, 4L)])
    } else if (!is.null(p$bbox) && length(p$bbox) == 4L &&
               all(is.finite(p$bbox)) && p$bbox[3L] > p$bbox[1L] &&
               p$bbox[4L] > p$bbox[2L]) {
        # No clip: the file's BoundingBox is what TDA declared the
        # plot area to be, and it is not the drawn extent -- a graph whose
        # nodes sit in the middle of its coordinate range leaves margin on
        # every side.  Shrink-wrapping to the content instead blows the
        # nodes and their labels up to fill the frame (gd31.cf: four nodes
        # spanning x 3 to 5 of a 0 to 8 range).
        xr <- range(xr, p$bbox[c(1L, 3L)])
        yr <- range(yr, p$bbox[c(2L, 4L)])
    }
    pad <- 0.04 * max(diff(xr), diff(yr), 1)

    if (newpage)
        grid::grid.newpage()
    # Keep TDA's proportions: a plot set up as 110 by 75 must not come back
    # square.  A respected layout gives the drawing area the aspect of the
    # coordinate range whatever the device shape is.
    #
    # diff(xr)/diff(yr) themselves, not the padded range, used to go
    # straight into the layout's widths/heights -- fine for an ordinary
    # plot, but a session with only a single point of content (one label,
    # nothing else to give the range any spread) has diff(xr) == 0, and a
    # zero-width/height grid.layout() is what grid calls a "non-finite
    # location and/or size for viewport", not a helpful message pointing
    # at the real cause. The pad already computed for xscale/yscale below
    # is applied here too, so the layout's proportions never collapse
    # to zero even in that degenerate case.
    grid::pushViewport(grid::viewport(
        width = 0.92, height = 0.92,
        layout = grid::grid.layout(1, 1,
            widths = grid::unit(diff(xr) + 2 * pad, "null"),
            heights = grid::unit(diff(yr) + 2 * pad, "null"),
            respect = TRUE)))
    grid::pushViewport(grid::viewport(
        layout.pos.row = 1, layout.pos.col = 1,
        xscale = c(xr[1L] - pad, xr[2L] + pad),
        yscale = c(yr[1L] - pad, yr[2L] + pad)))

    # TDA sizes text in the same units as everything else it draws: a node is
    # 11.34 points across and its number is set at 8.5, so the number fits
    # inside it.  The drawing is scaled to fill the device, so the text has to
    # be scaled by exactly the same factor or that relationship is lost --
    # anything else makes the labels overflow their nodes.
    #
    # This global span/avail pair is the right one for the single-panel
    # case (the overwhelming majority of files), used directly below.
    # For a multi-panel layout it is only a fallback default passed to
    # a caller's explicit cex=; the per-panel branch recomputes its
    # own cex fresh for each panel's own, much smaller sub-viewport
    # (see that branch's comment for why): using this same, global
    # cex inside a sub-viewport whose own available space is a fraction
    # of the whole page made every label enormous relative to its
    # panel, checked -- overlapping text across the entire
    # combined image, not merely oversized within its panel.
    # span/avail need to be paired by direction, not mixed: an
    # anisotropic plot (TDA's psetup() lets width/height and the
    # user coordinate ranges scale completely independently -- a
    # plot 150mm wide but only 15mm tall, checked against
    # an adversarially-constructed example built specifically to
    # trigger this) has a wide x-span paired with a wide device
    # width, and a short y-span paired with a short device height;
    # comparing the largest span (max(diff(xr), diff(yr)), here the
    # wide x-direction) against the smallest available dimension
    # (min of the device's width/height, here the short
    # y-direction) mixes two different directions' own numbers
    # together. That gave a cex far smaller than the plot's true
    # scale factor in either direction alone -- checked,
    # text that should have been a normal, readable 8.5pt rendered at
    # a small fraction of a point, illegible. The right comparison
    # pairs each axis with its own available space and takes
    # whichever axis is the tighter constraint.
    avail_x <- grid::convertWidth(grid::unit(0.92, "npc"), "points",
                                  valueOnly = TRUE)
    avail_y <- grid::convertHeight(grid::unit(0.92, "npc"), "points",
                                   valueOnly = TRUE)
    span_x <- diff(xr)
    span_y <- diff(yr)
    global_cex <- if (is.finite(span_x) && span_x > 0 && is.finite(avail_x) &&
                      is.finite(span_y) && span_y > 0 && is.finite(avail_y))
        max(0.5, min(6, min(avail_x / span_x, avail_y / span_y))) else 1
    # cex = NULL: the scale that keeps TDA's proportions, times the
    # option tdaR.ps.cex (default 1).  A document that finds TDA's 8.5 pt
    # labels heavy beside its prose sets the option once rather than
    # passing cex= to every plot() call.
    user_cex <- cex
    if (is.null(user_cex))
        cex <- global_cex * getOption("tdaR.ps.cex", 1)

    # A multi-panel layout (psetup() called more than once, or dplot()
    # combining several already-created files into a grid) has more
    # than one of TDA's coordinate systems sharing a single page --
    # tagged during parsing by op$panel (see tda_read_ps()'s
    # comment on why "%#Parameter:" is what that is tracked from).
    # Treating the whole file as one shared coordinate space (what
    # happens if every panel's ops are rendered directly into the
    # single, global viewport above, keyed only by the union of every
    # panel's extent) is wrong the moment there is more than one:
    # checked against Ghostscript's rendering of a real
    # 2x2 combined file, an evenly-proportioned grid of four panels,
    # not the wildly uneven layout naive rendering produced (one panel
    # dominating the canvas, the other three compressed into a
    # corner) -- because nothing tied each panel's own, much smaller
    # extent to its fair share of the overall space; whichever
    # panel happened to have the largest extent ended up occupying
    # most of the global scale's range, and the others were left
    # with whatever sliver of native-unit space their own, much
    # smaller extent implied. Giving each panel its own, separate
    # sub-viewport -- sized and positioned from its extent, in
    # "native" units of the outer, global xscale/yscale already
    # pushed above, which is exactly the coordinate space each
    # panel's extent already lives in -- reproduces the even,
    # separated grid a correct multi-panel rendering needs, without
    # touching anything about how a single panel (the overwhelming
    # majority of files) renders at all.
    panel_ids <- vapply(p$ops, function(o) o[["panel"]] %||% 1L, integer(1))
    panels <- unique(panel_ids)
    if (length(panels) <= 1L)
        .render_ps_ops(p$ops, cex)
    else {
        for (pnl in panels) {
            sub <- p$ops[panel_ids == pnl]
            ei <- extent(sub)
            xri <- ei$xr
            yri <- ei$yr
            if (!all(is.finite(c(xri, yri))))
                next
            padi <- 0.04 * max(diff(xri), diff(yri), 1)
            grid::pushViewport(grid::viewport(
                x = grid::unit(xri[1L] - padi, "native"),
                y = grid::unit(yri[1L] - padi, "native"),
                width = grid::unit(diff(xri) + 2 * padi, "native"),
                height = grid::unit(diff(yri) + 2 * padi, "native"),
                just = c("left", "bottom"),
                xscale = c(xri[1L] - padi, xri[2L] + padi),
                yscale = c(yri[1L] - padi, yri[2L] + padi)))
            # cex is the same, global one for every panel: font size
            # itself is already scaled by each panel's CTM at parse
            # time (see the "/fsiz" handler's comment for why),
            # putting it in the same, consistent device-point space as
            # the geometry: geometry, no matter which panel it belongs
            # to, is always exactly filling its sub-viewport (that
            # is what the sub-viewport's size, computed from this
            # same panel's extent, guarantees), so one, page-wide
            # cex is the correct scale for text too. Recomputing cex
            # per panel from its sub-viewport's available space
            # (an earlier version of this) double-counted the scaling
            # the fontsize fix already applies -- checked,
            # it left one panel enormous while the other three, whose
            # own extents happened to make the arithmetic cancel out
            # more forgivingly, looked correct by coincidence.
            .render_ps_ops(sub, cex)
            grid::popViewport(1L)
        }
    }
    grid::popViewport(2L)
    invisible(p)
}

# The actual per-op drawing loop, factored out of tda_plot_ps() so a
# multi-panel layout can call it once per panel's sub-viewport
# instead of duplicating it -- behaviour for a single panel (pushed
# straight into tda_plot_ps()'s outer viewport) is unchanged.
.render_ps_ops <- function(ops, cex) {
    for (o in ops) {
        # TDA relies on PostScript's clip operator to cut a series at the
        # plot boundary.  Recording the region is not enough -- the segments
        # have to be cut, or the line runs past the axis exactly as it does
        # in the file.
        if (!is.null(o$clip) && identical(o$op, "lines"))
            o <- .clip_line(o)
        # A symbol is a single point, not a path to cut -- it is either
        # inside its recorded clip region or it is not. This was
        # never checked at all before, so a scatterplot with any point
        # outside the plot's declared range (from noise pushing a
        # value past it, say) drew that point wherever its raw
        # coordinate happened to fall, including well outside the frame
        # -- checked by rendering deliberately out-of-range
        # data and seeing points outside the frame's border, not by
        # inspection.
        if (!is.null(o$clip) && identical(o$op, "symbol") &&
            (o$x < o$clip[1L] - 1e-9 || o$x > o$clip[2L] + 1e-9 ||
             o$y < o$clip[3L] - 1e-9 || o$y > o$clip[4L] + 1e-9))
            next
        if (is.null(o) || (identical(o$op, "lines") && length(o$x) < 2L))
            next
        gp <- grid::gpar(lwd = o$lwd, lty = if (is.numeric(o$lty)) "solid" else o$lty,
                         col = o$col, fontsize = o$fontsize * cex)
        switch(o$op,
            lines = if (length(o$x) > 1L) {
                if (is.numeric(o$lty)) {
                    # A real dash pattern (see .ps_dash()'s
                    # comment): grid's lty can only encode it as a
                    # line-width-relative hex string, which crushes
                    # TDA's true, absolute point lengths into
                    # indistinguishable micro-dashes for a thin line.
                    # Drawing the actual on-segments as solid runs,
                    # computed at their real length, is what makes the
                    # 9 line types separable the way the manual's
                    # reference figure shows them.
                    seg <- .dash_segments(o$x, o$y, o$lty)
                    if (length(seg$x) > 1L)
                        grid::grid.lines(seg$x, seg$y, default.units = "native", gp = gp)
                } else
                    grid::grid.lines(o$x, o$y, default.units = "native", gp = gp)
            },
            # A filled path is filled and not stroked: TDA does
            # "gsave / setgray fill / grestore / newpath", discarding the
            # path without a stroke.  Giving the polygon a border in the fill
            # colour draws an outline TDA never asked for.
            polygon = grid::grid.polygon(o$x, o$y, default.units = "native",
                                         gp = grid::gpar(col = NA,
                                                         fill = o$col)),
            # PostScript's show draws from the baseline, so TDA lowers a
            # centred label by about a third of the font size to put it in the
            # middle of a node.  grid.text centres vertically instead, so that
            # offset has to be undone or every label sits low.
            text = grid::grid.text(o$label, o$x,
                                   o$y + 0.32 * (o$fontsize %||% 8.5),
                                   default.units = "native",
                                   just = o$just, rot = o$rot %||% 0, gp = gp),
            # grid.circle() defaults to npc units, not native: without
            # default.units the centre is interpreted as a fraction of the
            # viewport and the circle lands off-screen.
            circle = .draw_arc(o),
            symbol = .draw_symbol(o))
    }
    invisible(NULL)
}

# Cut a polyline at a rectangular clip region, interpolating where a segment
# crosses the boundary so the line stops on it rather than at the last point
# inside.  Segments wholly outside are dropped.
.clip_line <- function(o) {
    b <- o$clip
    inside <- function(x, y) x >= b[1L] - 1e-9 & x <= b[2L] + 1e-9 &
                             y >= b[3L] - 1e-9 & y <= b[4L] + 1e-9
    x <- o$x; y <- o$y
    n <- length(x)
    if (n < 2L)
        return(if (inside(x, y)) o else NULL)
    keep <- inside(x, y)
    if (all(keep))
        return(o)
    cross <- function(x1, y1, x2, y2) {
        # walk the parameter range and take the part inside
        t0 <- 0; t1 <- 1
        for (k in 1:4) {
            p <- switch(k, x1 - x2, x2 - x1, y1 - y2, y2 - y1)
            q <- switch(k, x1 - b[1L], b[2L] - x1, y1 - b[3L], b[4L] - y1)
            if (abs(p) < 1e-12) {
                if (q < 0) return(NULL)
            } else {
                r <- q / p
                if (p < 0) { if (r > t1) return(NULL); if (r > t0) t0 <- r }
                else       { if (r < t0) return(NULL); if (r < t1) t1 <- r }
            }
        }
        c(x1 + t0 * (x2 - x1), y1 + t0 * (y2 - y1),
          x1 + t1 * (x2 - x1), y1 + t1 * (y2 - y1))
    }
    px <- py <- numeric(0)
    for (i in seq_len(n - 1L)) {
        seg <- cross(x[i], y[i], x[i + 1L], y[i + 1L])
        if (is.null(seg))
            next
        if (!length(px) || abs(px[length(px)] - seg[1L]) > 1e-9 ||
            abs(py[length(py)] - seg[2L]) > 1e-9) {
            px <- c(px, seg[1L]); py <- c(py, seg[2L])
        }
        px <- c(px, seg[3L]); py <- c(py, seg[4L])
    }
    if (length(px) < 2L)
        return(NULL)
    o$x <- px; o$y <- py
    o
}

# A full turn is a circle; anything less is an arc, which grid has no
# primitive for, so it is drawn as a polyline along the angle range.
.draw_arc <- function(o) {
    a1 <- o$a1 %||% 0
    a2 <- o$a2 %||% 360
    # A circle recorded with a clip region (see the extent() comment,
    # in tda_plot_ps() itself, on why this is tracked at all) needs
    # that clip applied to its drawing too, not only to the
    # viewport's sizing -- otherwise a circle far larger than the
    # plot itself (checked, examples/exam/gd18.cf's
    # rd=-3 edges) paints straight over the plot's frame and axis
    # instead of stopping at the boundary TDA itself clips it to.
    clipped <- !is.null(o$clip) && all(is.finite(o$clip))
    if (clipped)
        grid::pushViewport(grid::viewport(
            x = grid::unit(o$clip[1L], "native"),
            y = grid::unit(o$clip[3L], "native"),
            width = grid::unit(diff(o$clip[1:2]), "native"),
            height = grid::unit(diff(o$clip[3:4]), "native"),
            just = c("left", "bottom"),
            xscale = o$clip[1:2], yscale = o$clip[3:4], clip = TRUE))
    on.exit(if (clipped) grid::popViewport(1L))
    # A real dash pattern (see .ps_dash()'s comment) is a numeric
    # on/off array now, not a grid-compatible keyword. An earlier
    # version of this fell back to a solid outline whenever lty was
    # numeric, for every circle or arc regardless of shape -- noted in
    # a comment as a known gap but never actually finished. Confirmed
    # directly: examples/exam/plot14.cf's small dashed arc (the
    # angle marker between two lines) rendered solid instead. Fixed
    # the same way straight lines already are (.dash_segments()):
    # sample the curve as a dense polyline and dash-segment that,
    # rather than passing grid a dash value it cannot use on a curve
    # at all.
    numeric_lty <- is.numeric(o$lty)
    lty <- if (numeric_lty) "solid" else o$lty
    gp <- grid::gpar(lwd = o$lwd, lty = lty, col = o$col,
                     fill = o$fill %||% NA)
    # A uniform scale (m proportional to a rotation matrix) is a true
    # circle and grid.circle() draws it correctly and more cheaply; a
    # non-uniform one (plote's technique for an ellipse) is not, and
    # needs the actual linear transform applied point by point rather than
    # a single averaged radius, or every non-circular ellipse renders as
    # a plain circle -- confirmed by rendering one with distinct axis
    # lengths and finding exactly that before this fix.
    is_ellipse <- !is.null(o$m) && !is.null(o$r0) &&
        abs(o$m[1L]^2 + o$m[2L]^2 - (o$m[3L]^2 + o$m[4L]^2)) > 1e-6 * o$r0^2
    if (abs(a2 - a1) >= 359.9) {
        if (is_ellipse) {
            th <- seq(0, 360, length.out = 73)[-73] * pi / 180
            pts <- vapply(th, function(t)
                .ps_map(c(o$m, 0, 0), o$r0 * cos(t), o$r0 * sin(t)),
                numeric(2))
            grid::grid.polygon(o$x + pts[1L, ], o$y + pts[2L, ],
                               default.units = "native", gp = gp)
        } else if (numeric_lty) {
            th <- seq(0, 360, length.out = 145) * pi / 180
            seg <- .dash_segments(o$x + o$r * cos(th), o$y + o$r * sin(th),
                                  o$lty)
            if (length(seg$x) > 1L)
                grid::grid.lines(seg$x, seg$y, default.units = "native",
                                 gp = grid::gpar(lwd = o$lwd, col = o$col))
        } else {
            grid::grid.circle(o$x, o$y, r = o$r, default.units = "native",
                              gp = gp)
        }
    } else {
        # "270 180 arc" runs anticlockwise from 270 to 180, the short way, not
        # the long way round through 0.  PostScript's arc always sweeps
        # increasing angles, so a decreasing pair means the sweep passes 360.
        if (a2 < a1)
            a2 <- a2 + 360
        th <- seq(a1, a2, length.out = 48) * pi / 180
        ax <- o$x + o$r * cos(th)
        ay <- o$y + o$r * sin(th)
        if (numeric_lty) {
            seg <- .dash_segments(ax, ay, o$lty)
            if (length(seg$x) > 1L)
                grid::grid.lines(seg$x, seg$y, default.units = "native",
                                 gp = grid::gpar(lwd = o$lwd, col = o$col))
        } else {
            grid::grid.lines(ax, ay, default.units = "native",
                             gp = grid::gpar(lwd = o$lwd, lty = o$lty,
                                             col = o$col))
        }
    }
}

# Draws one of TDA's 17 numbered symbols (see .SYMBOL_TABLE), in native
# (PostScript-point) units the same way .draw_arc does, using the same
# relative-offset geometry as the PostScript procedures themselves
# (t_psinit.c) so the shapes match, not just approximate, TDA's.
.draw_symbol <- function(o) {
    d <- .SYMBOL_TABLE[[as.character(o$type)]] %||%
        list(shape = "circle", outline = TRUE)
    s <- o$size
    fillcol <- switch(d$fill %||% "none", white = "white",
                      black = o$col, NA)
    outline_col <- if (isTRUE(d$outline) || is.null(d$fill)) o$col else NA
    # See .draw_arc()'s comment: a real dash pattern is a numeric
    # array here, not something grid's lty can take directly, and
    # dash-segmenting a symbol's curved or multi-part outline is
    # not implemented, so this falls back to solid rather than
    # passing grid a value it cannot use.
    lty <- if (is.numeric(o$lty)) "solid" else o$lty
    gp <- grid::gpar(lwd = o$lwd, lty = lty, col = outline_col,
                     fill = fillcol)
    switch(d$shape,
        circle = grid::grid.circle(o$x, o$y, r = s, default.units = "native",
                                   gp = gp),
        square = grid::grid.polygon(o$x + s * c(-1, 1, 1, -1),
                                    o$y + s * c(-1, -1, 1, 1),
                                    default.units = "native", gp = gp),
        rhomb = grid::grid.polygon(o$x + s * c(-1, 0, 1, 0),
                                   o$y + s * c(0, 1, 0, -1),
                                   default.units = "native", gp = gp),
        # trian1/trian2: an up- and a down-pointing triangle, matching the
        # relative offsets in t_psinit.c's /trian1, /trian2 procedures.
        trian1 = grid::grid.polygon(o$x + s * c(-1, 1, 0),
                                    o$y + s * c(-1, -1, 1.5),
                                    default.units = "native", gp = gp),
        trian2 = grid::grid.polygon(o$x + s * c(-1, 1, 0),
                                    o$y + s * c(1, 1, -1.5),
                                    default.units = "native", gp = gp))
    mgp <- grid::gpar(lwd = o$lwd, lty = "solid", col = o$col)
    if (identical(d$mark %||% "", "cross") || identical(d$mark, "both")) {
        grid::grid.lines(o$x + s * c(-1, 1), o$y + s * c(0, 0),
                         default.units = "native", gp = mgp)
        grid::grid.lines(o$x + s * c(0, 0), o$y + s * c(-1, 1),
                         default.units = "native", gp = mgp)
    }
    if (identical(d$mark %||% "", "xsym") || identical(d$mark, "both")) {
        k <- s * cos(pi / 4)
        grid::grid.lines(o$x + k * c(-1, 1), o$y + k * c(-1, 1),
                         default.units = "native", gp = mgp)
        grid::grid.lines(o$x + k * c(-1, 1), o$y + k * c(1, -1),
                         default.units = "native", gp = mgp)
    }
}


#' Audit a parse against the PostScript it came from
#'
#' Counts the drawing operators in the file and the operations the parser
#' produced, and reports anything unaccounted for. Reading a plot back is easy
#' to get half right -- an operator that is silently ignored looks exactly
#' like one that is not there -- so this compares the two directly rather than
#' leaving it to the eye.
#'
#' @param file a \code{.ps} file, a TDA run, or a parse from
#'   \code{\link{tda_read_ps}}.
#' @param ... passed to \code{\link{tda_read_ps}}.
#' @return A data frame with one row per operator: how often it appears in the
#'   file, whether the parser handles it, and how many operations resulted.
#' @family plotting
#' @examples
#' p <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
#' p <- tda_pl_circle(p, at = c(5, 5), r = 2)
#' tda_check_ps(tda_ps_file(p))
#' @export
tda_check_ps <- function(file, ...) {
    parsed <- if (is.list(file) && !is.null(file$ops)) file
              else tda_read_ps(file, ...)
    txt <- readLines(parsed$file, warn = FALSE)
    body <- .ps_body(txt)
    body <- body[!grepl("^%", body)]

    # Every operator TDA's writers can emit, from the fprintf calls in
    # t_psf.c, t_plot.c, t_plot3.c, t_cplot.c, t_map.c, t_tree.c and the rest.
    known <- c("m", "l", "rm", "rl", "stroke", "fill", "closepath", "newpath",
               "clip", "gsave", "grestore", "arc", "show", "center", "adjust",
               "aright", "setgray", "setlinewidth", "setdash", "translate",
               "rotate", "scale", "scalefont", "setfont", "showpage",
               "circle", "cross", "xsym", "square", "trian1", "trian2",
               "rhomb", "invers", "FT", "FS", "sclear")
    # What this parser acts on.  The rest are recognised but ignored, which is
    # only safe for the ones that do not move anything.
    handled <- c("m", "l", "rm", "rl", "stroke", "fill", "closepath",
                 "translate", "rotate", "scale",
                 "arc", "show", "center", "adjust", "aright", "setgray",
                 "setlinewidth", "setdash", "gsave", "grestore",
                 "circle", "cross", "xsym", "square", "trian1", "trian2",
                 "rhomb")
    ignorable <- c("newpath", "clip", "showpage", "FT", "FS", "sclear",
                   "scalefont", "setfont", "invers")

    n <- vapply(known, function(op)
        sum(grepl(sprintf("(^|[^a-zA-Z])%s([^a-zA-Z]|$)", op), body)),
        integer(1))
    out <- data.frame(operator = known, in_file = unname(n),
                      handled = known %in% handled,
                      ignorable = known %in% ignorable,
                      stringsAsFactors = FALSE)
    out$unaccounted <- out$in_file > 0L & !out$handled & !out$ignorable
    kinds <- vapply(parsed$ops, `[[`, character(1), "op")
    attr(out, "produced") <- table(kinds)
    out[out$in_file > 0L, ]
}
