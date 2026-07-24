Copyrights and sources of code in tdaR
=======================================

TDA, the C program this package compiles, is
Copyright (C) 1989, 1991-2009 Goetz Rohwer (and, for the 6.x manual and
program, Ulrich Poetter), and is distributed under the GNU General Public
License, version 2 (file COPYING).  The R code, the vignettes, the build
and check scripts, and the changes to the C listed in
doc/changes-from-tda.md are Copyright (C) 2025-2026 Jan Marvin Garbuszus
and are distributed under the same licence.

TDA itself incorporates routines adapted from published algorithms and
from other programs.  This file records them, by source file, as the
source comments cite them.  Where a source carries a notice of its own it
is quoted.

Programs
--------

src/t_zoo.c, src/t_rzoo.c
    Reading zoo archives and LZD decompression, adapted from zoo, the
    archive package written by Rahul Dhesi (the archive header text is
    "ZOO 2.10 Archive").  The source comment reads: "Copyright is due to
    Rahul Dhesi.  See the copyright statement in the ZOO package."
    tdaR/src/tda_zoo.c (the package's own zoo()/unzoo()) follows the
    same directory-entry layout, including zoo 2.1's long-name
    extension and its rule for the 8.3 short name (zoo's dosname()).

src/t_tri.c
    Voronoi diagrams and Delaunay triangulations.  Several functions are
    "taken from Fortune's C program" (S.J. Fortune); the source comment
    reads: "The copyright is, of course, with Fortune."

src/t_gm.c  (polygon clipping)
    "Adapted from the MultClip program written and copyrighted (according
    to GPL) by Michael Leonov and Alexey Nikitin."

src/t_gm.c, src/t_top.c  (segment intersection and point-in-polygon)
    "Based on a proposal by Dan Sunday (April-B 2001 Algorithm)."

src/t_help.c  (s_match)
    Wildcard matching adapted from Mike Cornelison, "Two Wildcard
    Matching Utilities", C/C++ Users Journal 13(4), April 1995, 55-60.

src/t_lin.c
    Linear algebra routines marked "jack dongarra, linpack, 3/11/78".

src/t_gcmd.c  (machine constants)
    Adapted from the Fortran function DLAMCH (CMACH), "preliminary
    version March 26, 1990, Univ. of Tennessee, Oak Ridge National Lab,
    Argonne National Lab, Courant Institute, NAG Ltd., and Rice
    University" (the LAPACK routine; the source comment spells the
    package name "PAPACK").

src/t_gg.c  (gis, GRASP for maximum independent set)
    Adapted from a Fortran program by Mauricio G.C. Resende, Thomas A.
    Feo and Stuart H. Smith; see Feo, Resende, Smith, "A greedy
    randomized adaptive search procedure for maximum independent set",
    Operations Research 42 (1994), 860-878.  src/t_ass.c cites the same
    authors' code as ACM Algorithm 754 (M.G.C. Resende, AT&T Bell
    Laboratories).

src/t_nlreg.c  (orthogonal distance regression)
    Adapted from ACM Algorithm 676 (ODRPACK) by Paul T. Boggs, Richard H.
    Byrd, Janet R. Donaldson and Robert B. Schnabel (1987).

src/t_areg.c  (lowess)
    The LOWESS and LOWEST routines, robust locally weighted regression,
    carrying the routines' original header documentation; the source
    names no author.

src/t_loglin.c
    Log-linear fitting adapted from S.J. Haberman, "Log-Linear Fit for
    Contingency Tables", AS 51, Applied Statistics 21 (1972), 218-225;
    tests of marginal and partial association after E.D. Lustbader and
    R.K. Stodola, AS 160, Applied Statistics 30 (1981), 97-105, and M.B.
    Brown, "Screening Effects in Multidimensional Contingency Tables",
    Applied Statistics 25 (1976), 37-46; symbolic expansion of model
    formulas after M.J. Levine, CACM 13 (1970), 191-192.

src/t_int.c  (qtrap, qsimp)
    Numerical integration adapted from W.H. Press, B.P. Flannery, S.A.
    Teukolsky and W.T. Vetterling, Numerical Recipes in C, Cambridge
    University Press 1988, pp. 121-123.

Collected Algorithms of the ACM (CACM, TOMS)
--------------------------------------------

src/t_ass.c        548  G. Carpaneto, P. Toth, assignment problem
                   608  D.H. West, quadratic assignment (approximate)
                   754  GRASP, see t_gg.c above
src/t_cdf.c        304  normal distribution function
                   442  inverse normal distribution function
                   299  I.D. Hill, M.C. Pike, chi-squared integral
                   322  E. Dorrer, F-distribution
                   462  bivariate normal distribution
src/t_com.c        72, 371, 477  combinations, partitions in natural
                                 order, permutations
src/t_cplot.c      531  W.V. Snyder, contour plotting, TOMS 4(3)
src/t_gf.c         291  M.C. Pike, I.D. Hill, log gamma function
src/t_gg.c         457  C. Bron, J. Kerbosch, finding all cliques
src/t_gio.c        354  M.D. McIlroy, generator of spanning trees;
                        K. Paton, fundamental set of cycles, CACM 12 (1969)
src/t_gm.c         523  convex hull
src/t_graph.c           J.C. Tiernan, elementary circuits, CACM 13 (1970);
                        R.W. Floyd, shortest path
src/t_int.c        468  T.N.L. Patterson, automatic numerical integration
src/t_intp.c       624  R.J. Renka, triangulation and interpolation
                   526  H. Akima, bivariate interpolation
src/t_l1reg.c      478  I. Barrodale, F.D.K. Roberts, L1 norm regression
src/t_lp.c         333  R.C. Salazar, S.K. Sen, MINIT linear programming
                   449  linear programming in 0-1 variables
src/t_lsei.c       587  R.J. Hanson, K.H. Haskell, constrained least squares
src/t_mata.c,      575  I.S. Duff, row permutation for a zero-free diagonal
src/t_matc.c       529  I.S. Duff, J.K. Reid, block triangular form
src/t_mds.c        608, 754  see t_ass.c
src/t_min.c        178  A.F. Kaupe, direct search, with the remarks by
                        M. Bell and M.C. Pike (1966), R. de Vogelaere (1968),
                        F.K. Tomlin and L.B. Smith (1969)
src/t_rand.c       266  M.C. Pike, I.D. Hill, pseudo-random numbers (with
                        the 32-bit remark); 334 J.R. Bell, normal deviates;
                        369 Poisson deviates; 488 random numbers
src/t_spl.c        476  A.K. Cline, splines under tension; 433 H. Akima,
                        interpolation and smooth curve fitting
src/t_svd.c        343  singular value decomposition; 538 P.J. Nikolai,
                        eigenvalues
src/t_tmin.c       738  T. Chow, E. Eskow, R. Schnabel, tensor methods
                        for unconstrained optimization

Applied Statistics (Royal Statistical Society) algorithms
---------------------------------------------------------

src/t_cdf.c        AS 111  J.D. Beasley, S.G. Springer, percentage points
                           of the normal distribution
src/t_freq.c       AS 119  B.L. Leathers, sparse joint frequency
                           distributions
src/t_gf.c         AS 103  J.M. Bernardo, digamma; AS 121 B.E. Schneider,
                           trigamma; AS 147 Chi-Leung Lau, incomplete gamma;
                           R.J. Moore, derivatives of the incomplete gamma
                           integral, Applied Statistics 31 (1982)
src/t_loglin.c     AS 51, AS 160  see above (Haberman; Lustbader, Stodola)
src/t_svd.c        AS 60   D.N. Sparks, A.D. Todd, latent roots of a
                           symmetric matrix
src/t_tree.c       AS 13   minimum spanning tree

Data
----

inst/extdata/rrdat.1
    600 job episodes of 201 randomly selected respondents from the German
    Life History Study (GLHS), collected at the Max-Planck-Institut fuer
    Bildungsforschung, Berlin, in 1981-1983 (Mayer, K. U., Brueckner, E.
    (1989), Lebensverlaeufe und Wohlfahrtsentwicklung, Berlin: MPI fuer
    Bildungsforschung).  TDA's manual thanks Karl Ulrich Mayer and
    Hans-Peter Blossfeld for providing the data set, which has shipped
    with TDA since 1994 and is the example throughout the manual and
    Blossfeld and Rohwer, Techniques of Event History Modeling.

inst/extdata/exam/*.dat, deha1.zoo, tda.hlp
    TDA's own distribution files.
