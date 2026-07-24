#!/bin/sh
# Run the whole C corpus under AddressSanitizer.
#
# This needs NO instrumented R.  The standalone tda binary is the same C
# the package compiles, and tests/check.py drives it directly, so a
# second build of src/ is the whole setup.  That matters because
# building R with -fsanitize=address is tedious and rhub needs a public
# repository.
#
# Memory bugs here are worse than in the standalone program: the package
# runs TDA IN-PROCESS, so an overflow takes the user's R session down
# with no condition to catch.  Two such bugs (rxls, sdgshhs) were found
# only when a user handed over a real file.
#
# NOTE: tests/check.py discards stderr, so a case that dies under ASAN
# shows up only as "got <0 lines>".  Re-run that case directly, from its
# own suite directory, to see the report.
#
#     sh tools/check_asan.sh
set -e
d=${ASAN_BUILD:-/tmp/tda-asan}
rm -rf "$d"
cp -r src "$d"
(cd "$d" && rm -f ./*.o tda && make -j8 \
    CFLAGS="-std=gnu99 -g -O1 -fsanitize=address -fno-omit-frame-pointer \
        -fno-strict-aliasing -Wno-implicit-function-declaration \
        -Wno-unused-value -Wno-parentheses -Wno-format" \
    LDFLAGS="-fsanitize=address" >/dev/null)
ASAN_OPTIONS=detect_leaks=0 python3 tests/check.py "$d/tda" \
    examples/exam examples/ehhnew examples/tests examples/coverage
