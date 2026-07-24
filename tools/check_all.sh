#!/bin/sh
# Every gate, in the order that fails cheapest first.
#
# The slow ones are separate scripts and are NOT run here by default:
# ASAN takes a few minutes, memcheck outlasts a session on the exam
# suite.  Run those when the C changes, not on every edit.
#
#     sh tools/check_all.sh
#     sh tools/check_asan.sh        # after any C change
#     sh tools/check_ubsan.sh       # after any C change
#     sh tools/check_valgrind.sh    # slow; narrow with VG_SUITES
set -e
echo "== C build warnings ==";        sh tools/check_warnings.sh
echo "== run-time formats ==";       sh tools/check_rtfmt.sh
echo "== roxygen helpers ==";         python3 tools/check_roxygen_helpers.py
echo "== API argument names ==";      Rscript tools/check_api_names.R | tail -1
echo "== command dispatch ==";        python3 tools/check_dispatch.py
echo "== src/ vs tdaR/src/ ==";       python3 tools/check_src_sync.py
echo "== C reference suite ==";       python3 tests/check.py src/tda \
    examples/exam examples/ehhnew examples/tests examples/coverage | tail -1
echo "== every case exits cleanly =="
bin=$(pwd)/src/tda
bad=0
for s in examples/exam examples/ehhnew examples/tests examples/coverage; do
    for cf in "$s"/*.cf; do
        b=$(basename "$cf")
        (cd "$s" && timeout 120 "$bin" "cf=$b" >/dev/null 2>&1) || {
            echo "  $cf exited non-zero"; bad=$((bad + 1)); }
    done
done
echo "  $bad cases with a non-zero exit"
[ "$bad" -eq 0 ] || exit 1
echo "done"
