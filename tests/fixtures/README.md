# Test fixtures

Everything here is generated, so it carries the same licence as the rest of the
sources:

    tda.zoo, tda.zad   a TDA data archive, built by tools/make_archive.py
    electric.dta       a minimal Stata file (format 104)
    electric.dbf       a minimal dBase III file

The remaining input formats need real third-party files, which are **not**
distributed here because they are not ours to redistribute. To exercise the
commands that read them, put them in a directory outside the repository and
point `TDA_EXT_INPUT` at it, or use the default location `../TDA_ext_input`:

    electric.sav       SPSS system file    (readspss, inst/extdata)
    electric.por       SPSS portable file  (readspss, inst/extdata)
    deaths.xls         Excel workbook      (readxl,   inst/extdata)
    nc.shp/.shx/.dbf   ESRI shapefile      (sf,       inst/shape)
    sample.e00         ArcInfo interchange export
    deha1.zoo          the example archive from the TDA teaching pages
                       at stat.rub.de (48 command files and their data);
                       tdaR's zoo tests, `inst/examples/deha1.R` and the
                       `examples/coverage` LZH archive test read it from
                       here

Without that directory, `tests/smoke2.py` reports the affected commands under
"need external input" rather than counting them as failures, and
`tools/make_tests.py` leaves them out of the generated suite.
