# Generator for the interval_wages example data (reproducible; the
# function in R/data.R rebuilds it identically at run time, so nothing
# binary needs shipping).  Kept here as the reference recipe.
make_interval_wages <- function() {
    set.seed(1897)
    n <- 40L
    school_true <- runif(n, 8, 18)
    wage_true <- exp(0.8 + 0.11 * school_true + rnorm(n, 0, 0.25))
    # schooling reported in whole years: [floor, floor+1]
    school_lo <- round(school_true)   # schooling asked and answered in
    school_hi <- school_lo            # whole years: a point, realistically
    # income reported in brackets that widen with income, top bracket open
    br <- c(0, 5, 8, 12, 18, 27, 40)
    k <- findInterval(wage_true, br)
    wage_lo <- br[k]
    wage_hi <- c(br[-1], 80)[k]
    data.frame(wage_lo = wage_lo, wage_hi = wage_hi,
               school_lo = school_lo, school_hi = school_hi)
}
