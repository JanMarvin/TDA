# Build the pkgdown site.  The home page is the repository README; the
# vignettes are laid out by pkgdown/templates/content-article.html, which
# keeps their own stylesheet and scripts (tda-manual.css, theme-toggle.html)
# inside pkgdown's page frame.
#
#     source("pkgdown/build.R"); build_pkgdown_site()

build_pkgdown_site <- function(pkg = ".") {
  pkg <- normalizePath(pkg)
  readme <- file.path(pkg, "..", "README.md")
  index <- file.path(pkg, "pkgdown", "index.md")
  if (!file.exists(readme)) stop("README.md not found at ", readme)

  txt <- readLines(readme, warn = FALSE)
  repo <- "https://github.com/JanMarvin/TDA/blob/tdaR/"
  for (f in c("CONTRIBUTING.md", "ATTRIBUTIONS.md", "COPYING", "doc/changes-from-tda.md"))
    txt <- gsub(paste0("`", f, "`"), sprintf("[`%s`](%s%s)", f, repo, f),
                txt, fixed = TRUE)
  writeLines(txt, index)
  on.exit(unlink(index))

  pkgdown::build_site_github_pages(pkg, new_process = FALSE, install = FALSE)
  invisible(NULL)
}
