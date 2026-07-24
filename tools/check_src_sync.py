#!/usr/bin/env python3
"""Every file present in both src/ and tdaR/src/ must be byte-identical.

The C policy says any change to the C goes into BOTH copies.  Nothing
was checking, and three files (t_gen.c, tda_out.c, tda_rhooks.h) drifted
apart unnoticed for several sessions: the package build compiles the
guarded code, the standalone build compiles it away, so the 491-test
suite stayed green while the copies said different things.

    python3 tools/check_src_sync.py      # rc 1 if anything differs
"""
import os, sys, filecmp

A, B = "src", "tdaR/src"
bad = []
for fn in sorted(os.listdir(A)):
    if not fn.endswith((".c", ".h")):
        continue
    b = os.path.join(B, fn)
    if os.path.exists(b) and not filecmp.cmp(os.path.join(A, fn), b, shallow=False):
        bad.append(fn)
if bad:
    print("DIVERGED between src/ and tdaR/src/:")
    for fn in bad:
        print("   ", fn)
    sys.exit(1)
print("src/ and tdaR/src/ agree on every shared file")
