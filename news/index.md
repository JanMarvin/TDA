# Changelog

## tdaR 0.1.0

First release.

TDA 6.4 (Rohwer and Pötter, last release 6.4p of March 2009) compiled
into the package and run in-process; the build carries the banner 6.4q
because of the repairs listed in `doc/changes-from-tda.md` of the source
repository.

Wrapped: the episode-data and transition-rate commands (`edef`, `rate`
with every model, `ple`, `ltb`, `prate`, `epsdat`, `ejoin`, `frml`,
`freg`, `fml`); the regression commands (`lsreg`, `lsreg1`, `l1reg`,
`nlreg`, `npreg`, `glm`, `qreg` and its models, `ivreg`, `ireg`,
`loglin`); descriptive statistics, frequency tables, quantiles,
inequality and segregation measures; the sequence commands (`seqdef`,
`seqpd`, `seqsi`, `seqm`, `seqev`, `seqmd`, `seqgc`, `seqpm`); the
matrix language and the interval-arithmetic commands; graph analysis,
multidimensional scaling, clustering, assignment problems; the spatial
data commands; TDA’s PostScript plotting and a renderer for its output
in R; the readers for zoo archives, SPSS, Stata, Excel, DBF, ArcInfo E00
and shapefiles.

Results reach R through an export channel added to the C: every number a
function returns is TDA’s own value at double precision, never read back
from printed text.

Three vignettes: the TDA User’s Manual chapter by chapter, Blossfeld and
Rohwer’s *Techniques of Event History Modeling* command file by command
file, and the commands the manual does not cover.
