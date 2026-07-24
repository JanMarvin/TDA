#!/bin/sh
# Build Rohwer's own unmodified sources, for settling "is this ours?"
#
#   sh tools/build_pristine.sh <dir-of-original-sources> [out]
#
# The shipped makefile does not link: its object list is missing files
# (g_chull, in t_gm.c, among them) and it defaults to an X11 build whose
# sources need headers that are not in the archive. Compiling every .c
# except the two X11 ones and linking them together works.
#
# The result answers one question only, and answers it completely: run a
# command file through both binaries and the difference, if any, is ours.
# Reasoning from "our binary says so" is not the same thing and was wrong
# at least once (see the gdf entry in changes-from-tda.md).
set -e
SRC=${1:?usage: build_pristine.sh <src-dir> [out]}
OUT=${2:-$SRC/tda_orig}
cd "$SRC"
sed -i 's/^#define S_XWIN      1/#define S_XWIN      0/' tda.h 2>/dev/null || true
gcc -O2 -fsigned-char -c $(ls *.c | grep -vE '^t_xwin\.c$|^t_xlib\.c$')
gcc -o "$OUT" $(ls *.o | grep -vE '^t_xwin\.o$|^t_xlib\.o$') -lm
echo "built $OUT"
