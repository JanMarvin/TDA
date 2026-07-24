# Install every package in tdaR's Suggests, binaries only.
#   Linux (r2u): CRAN packages come from apt as r-cran-<name>;
#   Windows/macOS: CRAN binaries, plus readspss as a binary from r-universe.
pkgs <- strsplit(read.dcf("tdaR/DESCRIPTION", "Suggests")[1, 1], ",")[[1]]
pkgs <- trimws(sub("[[:space:]]*\\(.*\\)", "", pkgs))
pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (!length(pkgs)) quit(status = 0)
repos <- c("https://janmarvin.r-universe.dev", "https://cloud.r-project.org")
if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Linux") {
    cran <- setdiff(pkgs, "readspss")
    if (length(cran)) {
        r <- system2("sudo", c("apt-get", "install", "-y", "-qq",
                               paste0("r-cran-", tolower(cran))))
        if (r != 0) stop("apt-get failed")
    }
    if ("readspss" %in% pkgs)
        install.packages("readspss", repos = repos, type = "source")
} else {
    install.packages(pkgs, repos = repos, type = "binary")
}
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
