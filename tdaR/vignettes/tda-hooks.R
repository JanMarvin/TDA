# Shared by the three vignettes; each setup chunk source()s it.

knitr::opts_chunk$set(collapse = FALSE, comment = "", prompt = FALSE,
                      fig.width = 6.4, fig.height = 3.2, fig.align = "center",
                      dev = "svg", dev.args = list(bg = "transparent"),
                      out.width = "600px")

local({
    n <- 0L
    base <- knitr::knit_hooks$get("chunk")
    knitr::knit_hooks$set(chunk = function(x, options) {
        x <- base(x, options)
        if (is.null(options$box)) return(x)
        n <<- n + 1L
        sprintf(paste0('<p class="boxcap"><span class="lab">Box %d</span>',
                       '%s</p>\n<div class="box">\n%s\n</div>'),
                n, options$box, x)
    })
})

local({
    # Plots go into the page as inline SVG rather than an image, so the
    # stylesheet can restyle them: black strokes and text become the
    # page's ink, the grey band fills its band colour, in either setting.
    # cairo names its glyph outlines glyph-0-0, glyph-0-1, ... in every
    # file, so each figure's ids get a prefix or a later plot would draw
    # with the first plot's letters.
    #
    # The figure's width is out.width, a px size taken from the plot's
    # own PostScript dimensions, written as rem at the stylesheet's base
    # font size of 15px.  Everything else on the page is sized in rem
    # too, so the A-/A+ buttons, Cmd/Ctrl +/- and a trackpad pinch, all
    # of which move the root font size, scale the figure with the text.
    # The viewBox stays, so the height follows the width; the px width
    # and height attributes cairo writes would pin it and go.
    n <- 0L
    knitr::knit_hooks$set(plot = function(x, options) {
        n <<- n + 1L
        svg <- readLines(x, warn = FALSE)
        svg <- svg[!grepl("^<\\?xml", svg)]
        svg <- gsub("glyph-", sprintf("g%d-", n), svg, fixed = TRUE)
        w <- options$out.width
        if (is.null(w)) w <- "600px"
        if (grepl("^[0-9.]+px$", w)) {
            w <- sprintf("%.4grem", as.numeric(sub("px$", "", w)) / 15)
        }
        svg[1L] <- sub(' width="[0-9.]+" height="[0-9.]+"', "", svg[1L])
        svg[1L] <- sub("<svg ", sprintf('<svg class="plot" style="width:%s" ',
                                        w), svg[1L])
        paste0("<p class=\"figure\">\n", paste(svg, collapse = "\n"), "\n</p>")
    })
})

# cairo's SVG surface places every glyph at a whole-point x position, so
# at 8.5 pt the letters of a label drift by up to half a point and an
# enlarged plot reads "Y ax is".  The device is opened four times the
# figure size and displayed at the figure size through its viewBox, so
# that rounding is a quarter point on the page.  Text scales with the
# device by itself (tda_plot_ps sizes it from the drawing area); line
# widths do not, so every page starts with lex = 4.
local({
    K <- 4
    knitr::opts_hooks$set(dev = function(options) {
        options$fig.width <- options$fig.width * K
        options$fig.height <- options$fig.height * K
        options
    })
    setHook("grid.newpage",
            function() grid::pushViewport(grid::viewport(gp = grid::gpar(lex = K))),
            "replace")
})

options(width = 100, tdaR.ps.cex = 0.8)
