# External-fixture coverage tests

These `.cf` files exercise commands that need real third-party data
(`t_map.c`'s shapefile/GSHHS/DCW readers and drawing code) — not
synthesized fixtures. Per the same policy as `tests/fixtures/README.md`
and `tda-ext-fixtures/README.md`, none of that data is redistributed
here; that's why this directory is separate from `examples/coverage/`
and not part of the default `tools/check_all.sh` gate.

## Setup

Each `.cf` expects its fixture copied alongside it under the name it
opens, using the same files documented in `tda-ext-fixtures/README.md`:

- `mapshp.cf`, `mapshp_azi.cf`, `mapshp_azi2.cf` — `nc.shp`/`nc.shx`/
  `nc.dbf` from the R **sf** package (`inst/shape/nc.shp` in its own
  source tree), copied here as `shape.shp`/`shape.shx`/`shape.dbf`
  (TDA opens all three from one stem). 100 North Carolina counties,
  108 polygons. `mapshp.cf` uses the cylindrical projection (the same
  fixture and view/region/graticule already proven in tdaR's own
  `tda_map()` example, `R/spatial.R`); `mapshp_azi.cf`/`_azi2.cf` use
  the azimuthal orthographic projection and a geodesic-connected
  `sdpgeo` line, to reach the azimuthal-only code paths (`map_line_azi`,
  `gd_proj_inv`/`_p`, `mapgrid_*_azi*`, `mapgrid_bounds`/`_arc` via
  `cont=1`) that the cylindrical projection never touches.
- `mapgshhs.cf` — `gshhs_l.b` (GSHHS coastline database), copied here
  unchanged. Restricted to `level=1` and a small English Channel
  view/region so the polygon count stays manageable.
- `mapdcw.cf` — `zimbabwe.pnt` (DCW extract), copied here unchanged.
- `e00info.cf`, `e00pts.cf` — `test.e00` (the uploaded ArcInfo E00
  export, a water-well point coverage), copied here as `input.e00`.
  `e00info.cf` reads it in summary-only mode (no `df=`); `e00pts.cf`
  converts it to a TDA spatial file (`attr=2` also pulls in the
  attribute table, `dtda=` writes a TDA description file) and draws it
  with `sdpmap`.
- `e00poly.cf`, `e00single.cf`, `e00double.cf` — three more E00 sample
  files, from the same ESRI E00-format reference page the user's
  `test.e00` came from (its Appendix A/B worked examples), not
  invented: `poly.e00` (ARC+CNT+LAB+PAL+PRJ, a real polygon coverage —
  reaches `e00_arc`/`e00_cnt`/`e00_pal`/`e00_polygons`/`e00_prj`, none
  of which the point-only `test.e00` can reach), `single.e00`
  (ARC+LAB, no PAL — a line coverage, reaching `e00_lines`), and
  `double.e00` (the same shape as `single.e00` but using E00's
  double-precision ARC encoding, level 3 — one coordinate per line
  rather than level 2's four-per-line packed format, a genuinely
  different parse path inside `e00_arc`).
- `e00err.cf` — the error path (`e00_err`), reached by truncating a
  real E00 file mid-record (`truncated.e00`, the first 3 lines of
  `single.e00`) rather than fabricating a malformed one from scratch.
  `truncated.e00` ships alongside this README since it's a 3-line
  fragment, not the third-party file itself.

## Running

    cp $TDA_EXT_INPUT/nc_shape/nc.shp shape.shp   # (and .shx/.dbf)
    cp $TDA_EXT_INPUT/gshhs_l.b .
    cp $TDA_EXT_INPUT/zimbabwe.pnt .
    cp $TDA_EXT_INPUT/test.e00 input.e00
    cp $TDA_EXT_INPUT/e00-examples/poly.e00 .
    cp $TDA_EXT_INPUT/e00-examples/single.e00 .
    cp $TDA_EXT_INPUT/e00-examples/double.e00 .
    for f in *.cf; do /path/to/tda cf=$f; done

Coverage measured with these included is reported separately from the
`examples/coverage/` total in `doc/handover.md`, since it needs a setup
step `tests/check.py` doesn't do automatically.
