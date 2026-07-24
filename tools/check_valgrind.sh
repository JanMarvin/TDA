#!/bin/sh
# Run the C corpus under Valgrind's memcheck and write a log.
#
# Complements tools/check_asan.sh rather than repeating it.  ASAN is
# much faster and catches overflows of heap and stack objects; memcheck
# is slower but also reports READS OF UNINITIALISED MEMORY, which ASAN
# does not see at all.  Given TDA's age and its habit of computing with
# whatever a partly-filled array holds, that class is worth a look.
#
# Neither needs an instrumented R: the standalone tda binary is the same
# C the package compiles.
#
#     sh tools/check_valgrind.sh            # -> doc/valgrind.log
#     VG_SUITES="examples/exam" sh tools/check_valgrind.sh
#
# The log is a plain concatenation, one section per case, with only the
# cases that reported something kept.
#
# Cost, measured rather than guessed: about six minutes for all 491
# cases.  Two dominate (qr6 46s, rtd2m 34s); the rest average well under
# a second.  Narrow with VG_SUITES while working on one area.
#
# --track-origins=yes is on below.  It roughly triples the cost and only
# enriches the report for uninitialised values -- detection is the same
# without it -- so drop it for a quick sweep and re-run any case that
# reports.
#
# NOTE stdin is redirected from /dev/null.  Without it a case that reads
# stdin blocks until the per-case timeout -- 600s each -- which is what
# made the exam suite look like an overnight job when it is about five
# minutes.
set -e
command -v valgrind >/dev/null 2>&1 || {
    echo "valgrind is not installed (apt-get install valgrind)" >&2
    exit 1
}
d=${VG_BUILD:-/tmp/tda-valgrind}
log=${VG_LOG:-doc/valgrind.log}
suites=${VG_SUITES:-"examples/exam examples/ehhnew examples/tests examples/coverage"}

# -O0 with frame pointers so the traces name real functions; no
# sanitizer here, the two instrumentations do not mix.
rm -rf "$d"
cp -r src "$d"
(cd "$d" && rm -f ./*.o tda && make -j8 \
    CFLAGS="-std=gnu99 -g -O0 -fno-omit-frame-pointer -fno-strict-aliasing \
        -Wno-implicit-function-declaration -Wno-unused-value \
        -Wno-parentheses -Wno-format" >/dev/null)

bin=$(cd "$d" && pwd)/tda
: > "$log"
n=0
flagged=0
for s in $suites; do
    for cf in "$s"/*.cf; do
        [ -f "$cf" ] || continue
        n=$((n + 1))
        case_name=$(basename "$cf")
        out=$(cd "$(dirname "$cf")" && timeout 600 valgrind \
            --error-exitcode=0 --tool=memcheck --leak-check=no \
            --track-origins=yes --num-callers=20 \
            "$bin" "cf=$case_name" 2>&1 >/dev/null </dev/null || true)
        # Keep only cases that actually reported an error.  Matching
        # any "==pid== Word" line catches Valgrind's own banner and
        # heap summary, which every run prints -- that reported 138 of
        # 138 cases as findings.  The error count is the thing to test.
        if printf '%s' "$out" | grep -q "ERROR SUMMARY: [1-9]"; then
            flagged=$((flagged + 1))
            {
                echo "=============================================="
                echo "$cf"
                echo "=============================================="
                printf '%s\n\n' "$out"
            } >> "$log"
        fi
        printf '\r%d cases, %d with findings' "$n" "$flagged" >&2
    done
done
echo >&2
echo "$n cases run, $flagged with findings -> $log"
