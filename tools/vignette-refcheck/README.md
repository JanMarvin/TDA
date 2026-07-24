# Checking a vignette against TDA's own output

`.ref` files beside the command files record what TDA printed. These three
scripts compare the numbers in a knitted vignette's boxes against them.

    Rscript -e 'knitr::knit("vignettes/ehhnew.Rmd", output="/tmp/e.md")'

    # every number in a result table of the ref
    REFDIR=examples/ehhnew python3 cmp2.py /tmp/e.md 12:ehd2.ref 33:ehi1.ref

    # the Coeff column alone, which is the check that matters
    REFDIR=examples/ehhnew python3 coefchk.py /tmp/e.md 12:ehd2.ref 33:ehi1.ref

A spec is `boxes:reffile`; several boxes may feed one ref, as
`2,55,56:rt1.ref`. `boxes.py <file.md>` lists the box numbers and captions.

Matching is at the ref's own printed precision: a vignette value counts as
a hit if rounding it to the ref's decimal places reproduces the ref's.

Known non-findings, so they are not chased again:

- the minimiser's iteration protocol, memory figures and timestamps are
  filtered out of the ref; the nvar variable table is not entirely, so
  small integers (variable counts, column widths) show up as misses
- `Signif` in TDA against `Pr(>|t|)` in `summary()` is the same quantity
  reported the other way round
- a box using `head(coef(f), 8)` on a longer model will miss the rest
- `coefchk.py` takes the fourth-from-last number on a table row, so a ref
  line that is not a coefficient row can contribute a stray value

## Against the manual itself

`tman1.ps` (the manual, 2-up pages) becomes sequential text with

    gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=tman1.pdf tman1.ps
    pdftotext -layout -x 0   -y 0 -W 298 -H 842 tman1.pdf left.txt
    pdftotext -layout -x 298 -y 0 -W 297 -H 842 tman1.pdf right.txt
    # interleave left/right page by page (form feeds) into manual.txt

`cmp_manual.py` then takes every decimal number in a manual section and
looks for it in the knitted vignette's section, at the manual's printed
precision; `manbox.sh SECTION` prints the section's text for reading.
Both expect `manual.txt` at `/home/claude/man/manual.txt` and the
knitted vignette at `/tmp/vig/t.md` -- edit the two paths at the top.
Misses that are not errors: format codes (10.4, 13.10), section numbers
quoted in prose, and statistics TDA prints that the wrapper's print does
not (norm of residuals, log likelihood at the starting values).
