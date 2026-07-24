#!/bin/sh
# Run the C corpus under UndefinedBehaviorSanitizer.
#
# The third of the set, and it catches what the other two do not: signed
# overflow, shifts past the width of a type, misaligned or null pointer
# use, and conversions whose result is not representable.  ASAN sees
# out-of-bounds and use-after-free; memcheck adds uninitialised reads;
# UBSan is about operations that are undefined even when every address
# is valid.  That matters for code of TDA's age, which predates the
# compilers that exploit those rules.
#
# No instrumented R needed: the standalone tda binary is the same C.
#
#     sh tools/check_ubsan.sh                 # -> doc/ubsan.log
#     UB_SUITES="examples/exam" sh tools/check_ubsan.sh
#
# -fsanitize=undefined does not abort by default, it prints and carries
# on, so a case can report many findings and still produce correct
# output.  halt_on_error=0 keeps that behaviour so one case reports
# everything it hits rather than only the first.
#
# NOTE stdin is redirected from /dev/null.  Without it a case that reads
# stdin blocks until the per-case timeout -- 600s each -- which is what
# made the exam suite look like an overnight job when it is about five
# minutes.
set -e
d=${UB_BUILD:-/tmp/tda-ubsan}
log=${UB_LOG:-doc/ubsan.log}
suites=${UB_SUITES:-"examples/exam examples/ehhnew examples/tests examples/coverage"}

rm -rf "$d"
cp -r src "$d"
(cd "$d" && rm -f ./*.o tda && make -j8 \
    CFLAGS="-std=gnu99 -g -O1 -fsanitize=undefined -fno-omit-frame-pointer \
        -fno-strict-aliasing -Wno-implicit-function-declaration \
        -Wno-unused-value -Wno-parentheses -Wno-format" \
    LDFLAGS="-fsanitize=undefined" >/dev/null)

bin=$(cd "$d" && pwd)/tda
: > "$log"
n=0
flagged=0
for s in $suites; do
    for cf in "$s"/*.cf; do
        [ -f "$cf" ] || continue
        n=$((n + 1))
        case_name=$(basename "$cf")
        out=$(cd "$(dirname "$cf")" && \
            UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0 \
            timeout 300 "$bin" "cf=$case_name" 2>&1 >/dev/null </dev/null || true)
        if printf '%s' "$out" | grep -q "runtime error:"; then
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
