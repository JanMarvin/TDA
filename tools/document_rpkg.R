#!/usr/bin/env Rscript
# Regenerate man/ and NAMESPACE from the roxygen comments in R/.
#
# The committed man/ pages were written by hand, because roxygen2 was not
# available where the package was built.  The roxygen blocks in R/ are the
# source of record from here on; run this after changing them.
if (!requireNamespace("roxygen2", quietly = TRUE))
    stop("roxygen2 is needed: install.packages('roxygen2')")
roxygen2::roxygenise(".", clean = TRUE)
cat("man/ and NAMESPACE regenerated\n")
