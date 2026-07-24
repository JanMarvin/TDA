#!/bin/sh
# Build tdaR instrumented and run EVERY harness into one profile.
#
# The profile lives outside the repository, at $COV_DIR (default
# /home/claude/covbuild), so that packaging the tree does not delete it.
# Earlier runs kept it in tdaR/src and lost it to the zip step's
# *.gcda exclusion -- measure, ship, and the numbers were gone.
#
#   sh tools/coverage.sh            build + run everything + report
#   sh tools/coverage.sh report     report from the existing profile
#
# TDA_EXT_INPUT matters: every gated test skips without it, which is why
# t_e00.c reads 0% on a run that does not set it.  Point it at a
# directory holding gshhs_l.b, test.e00 and friends.
set -e
COV_DIR=${COV_DIR:-/home/claude/covbuild}
LIB=$COV_DIR/lib
SRC=$COV_DIR/src

if [ "$1" != "report" ]; then
    rm -rf "$COV_DIR"; mkdir -p "$LIB" "$SRC"
    mkdir -p "$HOME/.R"
    cp "$HOME/.R/Makevars" "$COV_DIR/Makevars.bak" 2>/dev/null || true
    printf 'CFLAGS = -std=gnu99 -g -O0 --coverage -fno-strict-aliasing\nLDFLAGS = --coverage\n' > "$HOME/.R/Makevars"
    rm -f tdaR/src/*.o tdaR/src/*.so tdaR/src/*.gcda tdaR/src/*.gcno
    R CMD INSTALL --no-docs --no-html -l "$LIB" tdaR >/dev/null 2>&1
    cp "$COV_DIR/Makevars.bak" "$HOME/.R/Makevars" 2>/dev/null || rm -f "$HOME/.R/Makevars"

    ( cd tdaR/tests && R_LIBS=$LIB TDA_EXAMPLES=$(cd ../../examples && pwd) \
        Rscript test-all.R >/dev/null 2>&1 ) || true

    #  The manual-chapter and example scripts in examples/rscripts.
    #  They use paths like ../examples/exam/x.dat, so they run from a
    #  scratch directory whose PARENT holds an examples symlink.
    RS=$COV_DIR/run; rm -rf "$RS"; mkdir -p "$RS/x"
    ln -sfn "$(cd examples && pwd)" "$RS/examples"
    cp examples/rscripts/*.R "$RS/x"/ 2>/dev/null || true
    ( cd "$RS/x" && R_LIBS=$LIB Rscript -e '
        suppressMessages(library(tdaR)); ok <- 0; bad <- 0
        for (f in list.files(".", pattern="[.]R$")) {
          e <- try(suppressWarnings(suppressMessages(
                 source(f, echo=FALSE, local=new.env()))), silent=TRUE)
          if (inherits(e,"try-error")) bad <- bad + 1 else ok <- ok + 1 }
        cat("rscripts run:", ok, " failed:", bad, "\n")' 2>&1 | tail -1 )

    #  THE ROXYGEN EXAMPLES.  117 of 118 man pages carry one, and they
    #  are working, documented exercises of the API -- but a coverage run
    #  that only executes the test suite never touches them.  Every
    #  measurement in this project before now understated the number for
    #  exactly that reason.  R CMD check runs them; so must this.
    R_LIBS=$LIB Rscript -e '
      suppressMessages(library(tdaR))   # examples assume the package is loaded
      ok <- 0; bad <- 0
      for (f in list.files("tdaR/man", pattern="[.]Rd$", full.names=TRUE)) {
        r <- try(suppressWarnings(suppressMessages(
               tools::Rd2ex(f, out <- tempfile(fileext=".R")))), silent=TRUE)
        if (inherits(r,"try-error") || !file.exists(out)) next
        e <- try(suppressWarnings(suppressMessages(
               source(out, echo=FALSE, local=new.env()))), silent=TRUE)
        if (inherits(e,"try-error")) bad <- bad + 1 else ok <- ok + 1
      }
      cat("examples run:", ok, " failed:", bad, "\n")' 2>&1 | tail -1
    R_LIBS=$LIB Rscript -e '
      suppressMessages(library(tdaR)); ok <- 0
      for (s in c("exam","ehhnew","tests","coverage"))
        for (f in list.files(file.path("examples",s), pattern="[.]cf$")) {
          d <- file.path(tempdir(),"c"); unlink(d, recursive=TRUE); dir.create(d)
          file.copy(list.files(file.path("examples",s), full.names=TRUE), d)
          if (!inherits(try(suppressWarnings(tda_run_cf(file.path(d,f))),
                            silent=TRUE), "try-error")) ok <- ok + 1 }
      cat("cf cases run:", ok, "\n")' 2>&1 | tail -1
    # Keep the profile where the zip step cannot reach it, then CLEAN
    # tdaR/src.  Leaving .gcda/.gcno behind makes R CMD check report
    # "Subdirectory 'src' contains" as a WARNING, which has bitten this
    # tree twice; the profile lives in $SRC and is read from there.
    cp tdaR/src/*.gcda tdaR/src/*.gcno "$SRC"/ 2>/dev/null || true
    cp tdaR/src/*.c tdaR/src/*.h "$SRC"/ 2>/dev/null || true
    rm -f tdaR/src/*.gcda tdaR/src/*.gcno tdaR/src/*.o tdaR/src/*.so
fi

cd "$SRC" && gcov -n *.c 2>/dev/null | awk '
/^File/ {f=$2; gsub(/\x27/,"",f)}
/^Lines executed:/ {split($0,a,":"); split(a[2],b,"% of "); p=b[1]+0; t=b[2]+0
  c=int(p*t/100+0.5)
  k = (f ~ /t_spss|t_xls|t_sdx|t_zoo|t_rzoo|t_e00/) ? "file IO" : "other"
  T[k]+=t; C[k]+=c; TT+=t; CC+=c
  printf "%7d %7d %5.1f %s\n", t-c, t, p, f > "/tmp/percov.txt" }
END { for (k in T) printf "%-8s %7d of %7d = %5.1f%%\n", k, C[k], T[k], 100*C[k]/T[k]
      printf "TOTAL    %7d of %7d = %5.1f%%\n", CC, TT, 100*CC/TT }'
echo "per-file detail: sort -rn /tmp/percov.txt | head"
