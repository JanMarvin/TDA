#!/bin/sh
# tda_rtfmt.c must reproduce the C library's output for every run-time
# format TDA builds; see tools/check_rtfmt.c.
set -e
d=$(mktemp -d)
trap 'rm -rf "$d"' EXIT
cd "$(dirname "$0")/.."
${CC:-cc} -std=c99 -O1 -I src -Wall -Wextra -Wno-format-nonliteral \
    -o "$d/check_rtfmt" tools/check_rtfmt.c -lm
"$d/check_rtfmt"
