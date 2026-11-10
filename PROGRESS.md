# moongeodesy progress

Updated: 2026-11-11 (Asia/Shanghai). Version 0.2.1 in `moon.mod`
(0.2.0 is on mooncakes.io). Toolchain: moon 0.1.20261104; CI installs the
current MoonBit CLI.

## Status

- `moon test --target all`: every test passes on wasm, wasm-gc, js, and
  native (142 per backend).
- `bash tools/acceptance.sh` runs every automated item of
  [`docs/十一月验收清单.md`](docs/十一月验收清单.md); all pass.
- Accuracy and reference data: see the table in [README.md](README.md#accuracy-and-testing).
  This file does not repeat the numbers, so they cannot go stale here.

## Packages

- `geodesic`: Karney geodesics (series and `GeodesicExact`), lines,
  polygons, rhumb lines, intersections, Geocentric and LocalCartesian
  (GeographicLib port).
- `utm`: transverse Mercator, polar stereographic, UTM/UPS, MGRS.
- `conic` (November): Lambert conformal conic and Albers equal-area, EPSG
  grids 9801/9802/9822.
- `magnetic` (November): World Magnetic Model 2025.
- `cn`: WGS-84 / GCJ-02 / BD-09 with exact inverses; no offset on foreign
  land or within about 9 km of foreign coasts.

## History

- October 2026: library created for the October contest (see README,
  History).
- 2026-11-08: fixes for the October acceptance review (#32).
- 2026-11-10: `magnetic` (#34) and `conic` (#36); README audit and
  acceptance checklist (#37); 0.2.0 published.
- 2026-11-10, second audit: foreign coasts and islands no longer offset
  (GeoNames-tested boundary), Albers poles round-trip exactly, prolate
  inverse conformal latitude converges, third-party notices, remaining
  README and docs corrections; 0.2.1.
- 2026-11-11: references that do not come from GeographicLib for the
  GeographicLib ports: rhumb lines (textbook formulas), `GeodesicExact`
  (quadrature, f from -0.5 to 0.5), polygon area (closed forms,
  L'Huilier, Green's theorem), intersections (known answers); all in
  `mpmath` at 45 digits.

## Ideas not done

- Independent references for the reduced length `m12` and the scales
  `M12`, `M21` of `GeodesicExact` (still checked against `GeodSolve -E`
  only), and for intersections in general position (the independent cases
  have answers that follow from symmetry).
- Rhumb lines within 0.1 deg of a pole keep about 3e-13 relative accuracy
  in longitude (1.8 um after 20 000 km circling the pole), because the
  conformal latitude is carried as an angle; a co-latitude form would
  keep full relative accuracy there.
- `Intersect::All` (every intersection within a radius).
