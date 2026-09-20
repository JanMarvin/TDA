# tdaR — TDA, Transition Data Analysis, as an R package

TDA is Götz Rohwer’s and Ulrich Pötter’s program for the statistical
analysis of event histories: transition rate models, parametric and
semi-parametric, with episode splitting, competing risks, multi-episode
data and time-dependent covariates, and around them a working
statistician’s toolkit — descriptive statistics, life tables and
product-limit estimation, regression of the ordinary kinds, sequence
analysis, loglinear models, graphs, spatial data, and its own PostScript
plotting. It is a command-file program: a text file of commands, a data
file, and a listing of what it computed.

tdaR is that program compiled into an R package. The C is TDA’s own,
made re-entrant so it can run inside R and be called more than once;
every result is computed by Rohwer’s code. Each R function writes the
command file TDA would have been given, runs it in-process, and returns
what TDA produced — the printed listing and the numbers behind it — as
an R object with the usual methods. This repository is maintained in
memory of Götz Rohwer.

``` r

library(tdaR)

rr  <- tda_rrdat()                        # the manual's example data (GLHS)
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr, model = "exponential")
fit$estimates
coef(fit); vcov(fit); logLik(fit)

tda_ple(Surv(TFP, DES) ~ factor(SEX), rr)   # product-limit estimates
tda_help("rate")                            # TDA's own manual entry
```

## Installation

``` r

# install.packages("remotes")
remotes::install_github("JanMarvin/TDA", ref = "tdaR", subdir = "tdaR")
```

The package compiles TDA’s C sources; a C compiler is needed (Rtools on
Windows, Xcode’s command-line tools on macOS). Nothing else is required.

## Documentation

Three vignettes:

- **Following the TDA User’s Manual** walks the manual chapter by
  chapter, with the manual’s own examples rewritten against the package,
  so that a page of the manual and the corresponding box of the vignette
  can be read side by side.
- **Techniques of Event History Modeling, in R** does the same for
  Blossfeld and Rohwer’s book, one command file at a time.
- **Beyond the Manual** covers what TDA grew after the manual stopped:
  spatial data, the readers for later file formats, and TDA’s online
  help.

Every exported function has a manual page with a runnable example;
[`?tdaR`](https://janmarvin.github.io/TDA/reference/tdaR-package.md)
gives the tour. `tda_help("<command>")` shows TDA’s help text for any
command under its TDA name.

## Layout of the repository

``` R
src/        TDA's C sources and the makefile for the standalone program
tdaR/       the R package; tdaR/src is a byte-identical copy of src/
examples/   TDA's shipped example command files with reference outputs
tests/      the C regression suite, its fixtures, and the vignette QA records
tools/      build and check scripts
doc/        changes-from-tda.md and the maintainer's notes
```

[`CONTRIBUTING.md`](https://github.com/JanMarvin/TDA/blob/tdaR/CONTRIBUTING.md)
describes building, the test gates, and how changes to the C are
handled.

## What has changed in the C

As little as possible. The rule is *reproduce, not improve*: TDA’s
output against its own reference runs is the contract, and every
apparent oddity is assumed intentional until a test shows otherwise.
Three kinds of change exist:

1.  the mechanical transformation that made the code a re-entrant
    library (every global became a member of one context struct;
    behaviour-neutral by construction and by test);
2.  package-only additions behind `TDA_R_PACKAGE`: an in-memory data
    matrix handed over from R, an export channel through which results
    reach R at full precision rather than through printed text, an
    interrupt check in the long searches, and error handling that
    returns to R instead of exiting the process;
3.  repairs of demonstrated faults, applied to both builds and listed
    one by one with their evidence in
    [`doc/changes-from-tda.md`](https://github.com/JanMarvin/TDA/blob/tdaR/doc/changes-from-tda.md).
    The version banner reads 6.4q so a patched build is distinguishable
    from Rohwer and Pötter’s last release.

## History

The account follows TDA’s manual (preface, 1998) and its own
`history.txt`, which ships with the package.

Rohwer began TDA in 1989, with early support from the Hamburger Institut
für Sozialforschung. The program grew inside Hans-Peter Blossfeld’s
research project *Household Dynamics and Social Inequality*, supported
by the European University Institute (Florence) and the European
Commission, and development then continued at the University of Bremen.
Version 5.7 (June 1994) was the first the authors called mature;
Blossfeld and Rohwer’s *Techniques of Event History Modeling* (1995; 2nd
ed. 2002) is written around it and shipped it on an accompanying disk.
The 6.x line was a redesign — a clean split of commands from parameters
— and from version 6.2 the program and manual have two authors, Rohwer
and Ulrich Pötter. The releases 6.1 to 6.4p (1997 to March 2009) fall in
Rohwer’s years as professor at Ruhr-Universität Bochum (1997–2013),
where TDA’s homepage still is. TDA has always been a non-commercial
project, released under the GNU General Public License. Its homepage is
<https://www.stat.rub.de/tda.html>.

## Licence

GPL-2, as TDA always was
([`COPYING`](https://github.com/JanMarvin/TDA/blob/tdaR/COPYING)). TDA
incorporates routines adapted from published algorithms and other
programs;
[`ATTRIBUTIONS.md`](https://github.com/JanMarvin/TDA/blob/tdaR/ATTRIBUTIONS.md)
(shipped in the package as `inst/COPYRIGHTS`) lists them by source file.
The example data, `rrdat.1`, are 600 job episodes from the German Life
History Study (Max-Planck-Institut für Bildungsforschung, Berlin; Mayer
and Brückner 1989), provided to TDA by Karl Ulrich Mayer and Hans-Peter
Blossfeld. The TDA User’s Manual is Rohwer and Pötter’s work and is
quoted only in brief, attributed excerpts.
