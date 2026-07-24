#!/usr/bin/env python3
"""Does our box print the same TABLE as the manual's, not merely the
same numbers somewhere in it?

    python3 tools/check_manual_shape.py <manual_pages.tsv> <html>

check_vignette_manual.py asks, for each number the manual prints,
whether it appears anywhere in the matching vignette box. That question
cannot see two faults, and both were found by hand rather than by it:

  a box showing half of what the manual shows -- every number present,
  each one once instead of twice, so membership is satisfied (6.5.5's
  epsdat box, which printed the single-episode table and not the
  multi-episode one)

  a box showing a DIFFERENT table whose numbers happen to overlap --
  small state codes and times drawn from the same small pool (3.4.3's
  seqpd box, which ran one sequence structure where the manual runs two)

So this compares counts. For every number in the manual's box, how many
times does it appear there, and how many times in ours? A box that
prints half the manual's rows shows up as a long list of tokens wanted
twice and found once. Exact equality is not the bar -- our tables carry
row numbers and column headers the manual has not got, and the manual
truncates long tables with "...." -- so what is reported is the tokens
we print FEWER times than the manual does, which is the direction that
means something is missing.
"""
import collections, re, sys, os

here = os.path.dirname(os.path.abspath(__file__))
ns = {"__name__": "x"}
sys.argv = sys.argv[:3]
exec(open(os.path.join(here, "check_vignette_manual.py")).read()
     .rsplit("\nmain()", 1)[0], ns)

pages, boxes_of_section = ns["pages"], ns["boxes_of_section"]
numbers, kind, section_boxes = ns["numbers"], ns["kind"], ns["section_boxes"]
vignette_boxes = ns["vignette_boxes"]

doc = open(sys.argv[2], errors="replace").read()
bysec, vb = section_boxes(doc), vignette_boxes(doc)

# Boxes investigated and written up in doc/changes-from-tda.md: the difference is
# explained there, so they are listed apart rather than kept at the top
# of a queue of work. Each one was checked against the program, and in
# five cases against a build of Rohwer's own sources.
SETTLED = {
    "6.14.2 Box 2": "sig is not identified on these data",
    "5.1.4.11 Box 1": "defective eigenvalues, machine-dependent",
    "5.1.4.11 Box 2": "defective eigenvalues, machine-dependent",
    "5.1.4.11 Box 3": "defective eigenvalues, machine-dependent",
    "6.12.5 Box 2": "manual out of date (checked against Rohwer's binary)",
    "8.4.1 Box 2": "manual out of date (checked against Rohwer's binary)",
    "7.2.1.1 Box 1": "manual out of date (checked against Rohwer's binary)",
    "6.15.2.1 Box 4": "TDA reads its data file as 4-byte floats",
    "6.2.2 Box 12": "gdf4.cf round-trips its data through a [10.6] file",
    "7.6.1.2 Box 6": "the right-hand table is reproduced exactly from "
                     "gd11.cf; the left-hand one (IT 1) comes from a run "
                     "that is not among the shipped command files, and "
                     "mxit = 1 gives 0.505 where it prints 0.5000",
    "6.7.2.3 Box 12": "the D matrix is identical; the manual's box also "
                      "prints its row and column indices 0..8 as headers, "
                      "where ours labels by the sequence values",
    # the same reason for all three dtda boxes
    "6.12.2 Box 6": "dtda widths: we ask for full precision ([24.16]) "
                    "where the manual's run used [10.4]",
    "6.12.3 Box 5": "dtda widths: we ask for full precision ([24.16]) "
                    "where the manual's run used [10.4]",
    "6.12.4 Box 5": "dtda widths: we ask for full precision ([24.16]) "
                    "where the manual's run used [10.4]",
}

rows, settled = [], []
for sec in pages:
    if not bysec.get(sec):
        continue
    # Counting has to be numeric, not textual: the manual writes 0.00
    # where we write 0, and keying a counter on the token itself reports
    # every such pair as missing. The first version of this did exactly
    # that and made 6.5.1 look 165 short when its values are identical.
    pool = []
    for n in sorted(bysec[sec]):
        pool += [float(t) for t, _ in numbers(vb.get(n, ""))]
    for num, cap, body in boxes_of_section(sec):
        if kind(cap) != "output":
            continue
        # A box's caption runs to the next caption, so the prose between
        # boxes lands inside it: citations, page numbers, "convergence
        # in 12 function calls". Those are not output and no vignette
        # box will ever print them. A line of eight or more words is
        # prose, not a table row.
        # TDA echoes the settings it was given before it estimates --
        # the algorithm, the iteration cap, the tolerances. The residual
        # classifier counts those as settings rather than results, and
        # so does this: a vignette box prints the fit, not the echo.
        SETTING = re.compile(
            r"^\s*(Algorithm|Number of model parameters|Type of covariance|"
            r"Maximum number of|Convergence criterion|Tolerance|Mue of |"
            r"Minimum of step|Scaling factor|Method:|Model:)", re.I)
        want = []
        for line in body.split("\n"):
            if len(re.findall(r"[A-Za-z]{2,}", line)) >= 8:
                continue
            if SETTING.match(line):
                continue
            parts = line.split()
            if len(parts) >= 5 and all(re.fullmatch(r"-?\d+", x) for x in parts):
                step = [int(parts[i + 1]) - int(parts[i])
                        for i in range(len(parts) - 1)]
                if len(set(step)) == 1 and step[0] > 0:
                    continue
            want += [(float(t), d) for t, d in numbers(line)]

        # Two passes, exact first. A coarse match must not consume a
        # value that an exact match needs: the manual writes 1.0308e-01
        # where a fit prints 0.1031, and matching that loosely first
        # could eat the 0.1031 another token wanted exactly.
        left = sorted(pool)
        unmatched = []
        for x, d in want:
            tol = 0.5 * 10 ** (-d) + 1e-12
            hit, best = None, None
            for k, y in enumerate(left):
                dy = abs(y - x)
                if dy <= tol and (best is None or dy < best):
                    hit, best = k, dy
            if hit is None:
                unmatched.append(x)
            else:
                del left[hit]
        short = coarse = 0
        for x in unmatched:
            # four significant figures: our printed precision against
            # the manual's is the usual reason a value does not match
            hit, best = None, None
            for k, y in enumerate(left):
                dy = abs(y - x)
                if dy <= max(abs(x), abs(y)) * 1e-3 + 1e-9 and \
                   (best is None or dy < best):
                    hit, best = k, dy
            if hit is None:
                short += 1
            else:
                coarse += 1
                del left[hit]

        if short:
            key = "%s Box %s" % (sec, num)
            (settled if key in SETTLED else rows).append(
                (sec, num, short, coarse, len(want), cap[:40]))

rows.sort(key=lambda r: -r[2])
print("boxes printing some token fewer times than the manual: %d\n" % len(rows))
print("%-10s %-7s %6s %6s %6s  %s"
      % ("section", "box", "short", "prec", "total", "caption"))
for sec, num, short, kinds, tot, cap in rows[:30]:
    print("%-10s Box %-3s %6d %6d %6d  %s" % (sec, num, short, kinds, tot, cap))
if settled:
    print("\nand %d box(es) already settled, with the reason:" % len(settled))
    for sec, num, short, _, _, _ in sorted(settled):
        print("  %-10s Box %-3s %4d short  %s"
              % (sec, num, short, SETTLED["%s Box %s" % (sec, num)]))
