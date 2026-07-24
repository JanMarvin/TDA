library(testit)

# One pass, on the export path: every number a wrapper returns comes
# through TDA's export channel, and test-audit-precision.R checks that
# it does.  The text-only path (options(tdaR.use_exports = FALSE)) is
# not run any more; it exists for diagnosing the exporters, not as a
# second implementation to keep green.
test_pkg("tdaR")
