# Contributing

## Building

The standalone program:

    cd src && make            # gcc; make CC=clang also works

The R package:

    R CMD INSTALL tdaR

`tdaR/src` is a committed copy of `src/`. A change to the C goes into
both; `python3 tools/check_src_sync.py` fails on any difference.

The tests and the vignettes read the examples tree; set `TDA_EXAMPLES`
to `examples/` if it is not found from the working directory. Tests
that need third-party input files (SPSS, Excel, shapefiles, ArcInfo,
the `deha1.zoo` archive) run only when `TDA_EXT_INPUT` points at a
directory holding them; see `tests/fixtures/README.md`.

## Checks before a change is done

    sh tools/check_warnings.sh
    python3 tools/check_src_sync.py
    python3 tests/check.py src/tda examples/exam examples/ehhnew \
        examples/tests examples/coverage
    cd tdaR/tests/testit && Rscript -e 'testit::test_pkg("tdaR", dir = ".")'
    R CMD build tdaR && R CMD check --as-cran tdaR_*.tar.gz

`check_warnings.sh` builds with gcc and clang under `-Wall -Wextra
-Wpedantic -Wshadow -Wconversion` and the rest of the flag set it
lists; no `#pragma GCC diagnostic` anywhere. `tests/check.py` compares
the program's output with TDA's reference runs at printed precision;
`--regen` rewrites a reference only when its compared content changed.

When a vignette or anything feeding one changed, render the three
vignettes and run the vignette checks:

    for v in tdaR ehhnew beyond-the-manual; do
      Rscript -e "rmarkdown::render('tdaR/vignettes/$v.Rmd', output_dir = 'vout')"
    done
    python3 tools/vignette_fingerprint.py check \
        tests/vignette-qa/vignette-fingerprint.tsv vout/tdaR.html vout/ehhnew.html
    Rscript tools/check_gallery.R
    Rscript tools/check_vignette_ps.R

The fingerprint records the numbers every box prints; a change is
either a fault or an intended improvement, and is re-frozen with
`vignette_fingerprint.py write`. `check_gallery.R` and
`check_vignette_ps.R` compare the PostScript the gallery scripts and
the vignette plots produce with TDA's own `examples/exam/*.ps` byte
for byte.

When `t_zoo.c` or `tdaR/src/tda_zoo.c` changed, build the reference
zoo (`git clone https://github.com/troglobit/zoo`, `./autogen.sh &&
./configure && make`) and run

    Rscript tools/check_zoo_ref.R /path/to/zoo/src/zoo

which archives long names and directories with both writers and
compares extraction and every directory-entry name field.

## The C

The C is Rohwer and Pötter's. Reproduce, not improve: TDA's output
against its own reference runs is the contract, and an apparent oddity
is assumed intentional until a test shows otherwise. Three kinds of
change exist:

1. the context-threading that made the code re-entrant: every global is
   a member of `TDAContext` (`src/tda_context.h`), every function takes
   `TDAContext *ctx` first;
2. package-only additions behind `TDA_R_PACKAGE`: `rdataframe` (an
   in-memory data matrix from R), the export channel (`tda_export_row()`
   beside the printers, so results reach R at full precision), the
   interrupt check in the long searches, and `exit()` turned into a
   `longjmp` back to the `.Call` entry;
3. repairs of demonstrated faults, applied to both builds, each with a
   control file under `examples/tests/` and an entry in
   `doc/changes-from-tda.md`.

TDA has no error channel: errors are printed lines and the exit code
is always 0. `TDAContext.ErrCnt` counts them and the R side returns
them as `$errors`. What a TDA command does is answered from the C and
the manual, not by probing the binary.

## The R side

Use TDA's own command wherever TDA has one, never an R stand-in.
Results come through the export channel, never by parsing numbers from
printed text; the console text may locate a line, not supply a value.
Output files are requested only when the user asks for one (PostScript,
the SPSS/Stata/CSV exports). Every fix leaves a test or a control file.

Nothing that is not TDA's own ships in the package. The external test
fixtures and Blossfeld and Rohwer's book (the `ehhnew` vignette is
checked against it) stay outside the tree; `tools/check_vignette_book.py`
takes the book's text as an argument and embeds none of it.

## Vignette conventions

Box captions name the command file a box reproduces (`binary logit
(qr1.cf)`); `tools/check_vignette_refs.py` and the other vignette
checkers key on that, never on box numbers, which shift when a box is
added. A chunk that draws is folded, a chunk that fits and then draws is
split into an unfolded fit box and a folded plot box. The three
vignettes share `vignettes/tda-hooks.R` (box and figure hooks; figure
widths in rem so the A−/A+ buttons and the RStudio viewer zoom scale
plots with the text), `tda-manual.css` and `theme-toggle.html`.

## The site

`tdaR/pkgdown/build.R` builds the pkgdown site; the vignettes are laid
out by `tdaR/pkgdown/templates/content-article.html`, which keeps their
own stylesheet and scripts inside pkgdown's page frame.
