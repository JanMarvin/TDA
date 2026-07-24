#!/bin/sh
# The C must build without a single warning under the project's flag set,
# in the standalone and in the package configuration.  -O2 so that gcc's
# flow-based warnings (-Wmaybe-uninitialized) take part.
#
#     sh tools/check_warnings.sh            # rc 1 if anything warns
#     CC=clang sh tools/check_warnings.sh
set -e
cc=${CC:-cc}
w="-Wall -Wextra -Wpedantic -Wshadow -Wunused -Wconversion -Wsign-compare \
-Wcast-align -Wnull-dereference -Wdouble-promotion -Wformat=2 \
-Wmisleading-indentation -Wno-ignored-attributes"
d=$(mktemp -d)
trap 'rm -rf "$d"' EXIT
cp src/*.c src/*.h src/makefile "$d"/
(cd "$d" && make -j8 CC="$cc" CFLAGS="-std=c99 -O2 -ffp-contract=off \
    -fno-strict-aliasing $w") > "$d/build.log" 2>&1 || { cat "$d/build.log"; exit 1; }
rinc=$(R CMD config --cppflags 2>/dev/null || true)
for f in tdaR/src/*.c; do
    $cc -O2 -c -o "$d/r.o" $rinc -DNDEBUG -DTDA_R_PACKAGE $w "$f"
done >> "$d/build.log" 2>&1
if grep -q "warning:" "$d/build.log"; then
    grep -A3 "warning:" "$d/build.log"
    echo "C warnings under the project's flag set"
    exit 1
fi
echo "C builds clean under the project's flag set ($cc, standalone and package)"
