#!/bin/sh
# Incremental builds of the TDA standalone in persistent /tmp trees, so
# a one-file change costs one compile plus a link (~1 s), not a full
# rebuild -- for the sanitizer and coverage variants too.
#
#   tools/quickbuild.sh asan                 sync + incremental ASAN build
#   tools/quickbuild.sh cov t_matc.c ...     instrument ONLY those files;
#                                            everything else built plain
#
# Binaries land in /tmp/tda-asan/tda and /tmp/tda-cov/tda.  The plain
# binary needs no script: `make -j8` in src/ is already incremental.
#
# cov mode: run the suite with /tmp/tda-cov/tda, then
#   cd /tmp/tda-cov && gcov -f -n t_matc.c
# Only the named files produce coverage data; the rest of the binary
# runs uninstrumented (and faster).  Re-running with different files
# re-instruments only the difference.  Delete stale counts with
#   rm -f /tmp/tda-cov/*.gcda
# before a fresh measurement.

set -e
mode=$1; shift || true
root=$(cd "$(dirname "$0")/.." && pwd)

BASE="-std=gnu99 -DS_UNIX=1 -DS_DOS=0 -DS_XWIN=0 -ffp-contract=off \
-fno-strict-aliasing -Wno-implicit-function-declaration -Wno-unused-value \
-Wno-parentheses -Wno-format"

case "$mode" in
asan)
    d=/tmp/tda-asan
    mkdir -p "$d"
    cp -u "$root"/src/*.c "$root"/src/*.h "$root"/src/makefile "$d"/
    cd "$d"
    make -j8 CFLAGS="$BASE -g -O1 -fsanitize=address -fno-omit-frame-pointer" \
         LDFLAGS="-fsanitize=address" > /dev/null
    echo "$d/tda ready (run with ASAN_OPTIONS=detect_leaks=0)"
    ;;
cov)
    [ $# -ge 1 ] || { echo "cov mode needs file names, e.g. t_matc.c" >&2; exit 1; }
    d=/tmp/tda-cov
    mkdir -p "$d"
    cp -u "$root"/src/*.c "$root"/src/*.h "$root"/src/makefile "$d"/
    cd "$d"
    # everything plain and incremental first ...
    make -j8 CFLAGS="$BASE -O1" > /dev/null 2>&1 || true
    # ... then force-instrument only the named files and relink.
    for f in "$@"; do
        rm -f "${f%.c}.gcda" "${f%.c}.gcno"
        cc $BASE -O0 --coverage -c "$f" -o "${f%.c}.o"
    done
    cc $BASE --coverage -o tda *.o -lm
    echo "$d/tda ready (instrumented: $*)"
    ;;
*)
    echo "usage: $0 asan | cov file.c [file.c ...]" >&2
    exit 1
    ;;
esac
