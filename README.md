# moongeodesy

[![CI](https://github.com/ChrisTu117/moongeodesy/actions/workflows/ci.yml/badge.svg)](https://github.com/ChrisTu117/moongeodesy/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-mooncakes.io-green)](https://mooncakes.io/docs/ChrisTu117/moongeodesy)

Pure MoonBit geodesy: nanometre-accurate geodesics on the ellipsoid
(a port of Charles Karney's [GeographicLib](https://geographiclib.sourceforge.io/)
algorithms), UTM/UPS and MGRS grid references, Lambert conformal and Albers
equal-area conic projections, the World Magnetic Model 2025 for compass
declination, and the WGS-84 ↔ GCJ-02 ↔ BD-09 conversions needed for maps in
mainland China. No FFI; works on the `wasm`, `wasm-gc`, `js` and `native`
backends.

> 中文简介见[下方](#中文简介)。

| Package | What it does |
| --- | --- |
| `ChrisTu117/moongeodesy/geodesic` | Karney's geodesic algorithms (series and `GeodesicExact` elliptic integrals), rhumb lines, geodesic intersections, and local ENU (`LocalCartesian`): inverse & direct problems, geodesic lines, polygon area/perimeter, geocentric ECEF on any ellipsoid |
| `ChrisTu117/moongeodesy/utm` | Transverse Mercator (6th-order Krüger series), polar stereographic, UTM/UPS with standard zones, MGRS encode/decode |
| `ChrisTu117/moongeodesy/conic` | Lambert conformal conic and Albers equal-area conic (one or two standard parallels, any ellipsoid), EPSG-style grids with a false origin |
| `ChrisTu117/moongeodesy/magnetic` | World Magnetic Model 2025: declination, inclination, field components and their yearly change, grid variation, compass Blackout/Caution zones |
| `ChrisTu117/moongeodesy/cn` | WGS-84 → GCJ-02 (Gaode/Tencent), GCJ-02 → BD-09 (Baidu), with exact iterative inverses |

## Install

```bash
moon add ChrisTu117/moongeodesy
```

`conic` and `magnetic` need version 0.2.0 or later. Import the packages you
need in your `moon.pkg`:

```
import {
  "ChrisTu117/moongeodesy/geodesic",
  "ChrisTu117/moongeodesy/utm",
  "ChrisTu117/moongeodesy/conic",
  "ChrisTu117/moongeodesy/magnetic",
  "ChrisTu117/moongeodesy/cn",
}
```

## Usage

The comments state values; [`readme_claims_test.mbt`](src/examples/readme_check/readme_claims_test.mbt) asserts each of them.

```moonbit
let g = @geodesic.Geodesic::wgs84()

// Inverse problem: distance and bearings between two points
let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737) // Beijing → Shanghai
// r.s12 = 1066676.4 m, r.azi1 = 152.94°, r.azi2, r.m12, r.area12 ...

// Direct problem: go 1000 km due east from Shenzhen
let d = g.direct(22.5431, 114.0579, 90.0, 1000000.0) // d.lat2, d.lon2, d.azi2

// Geodesic line: waypoints at most 100 km apart, e.g. for drawing a route
let line = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
let pts = line.densify(100000.0)

// Polygon area (m², counter-clockwise positive) and perimeter (m)
let a = g.polygon_area([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])

// Rhumb line: constant-azimuth distance, London → Tokyo (longer than the geodesic)
let rh = @geodesic.Rhumb::wgs84()
let h = rh.inverse(51.4778, -0.0015, 35.6762, 139.6503) // h.s12, h.azi12, h.area12

// Intersection of two geodesics (signed metres from each start)
let ix = @geodesic.Intersect::wgs84()
let cross = ix.closest(0.0, 0.0, 45.0, 45.0, 10.0, 135.0) // cross.x, cross.y, cross.c

// Local East-North-Up relative to an origin (ellipsoidal height in metres)
let lc = @geodesic.LocalCartesian::wgs84(39.9087, 116.3975)
let enu = lc.forward(39.9097, 116.3985, 50.0) // enu.x east, enu.y north, enu.z up
let back = lc.reverse(enu.x, enu.y, enu.z) // back.lat, back.lon, back.h

// Other ellipsoids: CGCS2000 / GRS80 built in, or any (a, f)
let cgcs = @geodesic.Geodesic::cgcs2000()
let sphere = @geodesic.Geodesic::new(6371000.0, 0.0)

// High flattening (|f| above about 1/50): elliptic-integral geodesic
let flat = @geodesic.GeodesicExact::new(6378137.0, 0.1)
let hx = flat.inverse(40.6, -73.8, 1.5, 0.0) // hx.s12

// Chinese map coordinate systems
let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → Gaode
let bd = @cn.gcj02_to_bd09(gcj) // Gaode → Baidu
let wgs = @cn.bd09_to_wgs84(bd) // exact inverse (round trip < 1e-9°)

// UTM / UPS / MGRS (functions raise @utm.UtmError on bad input)
let u = @utm.from_latlon(39.9087, 116.3975) // zone 50N, easting/northing, convergence, scale
let m = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
let p = @utm.from_mgrs("50S MK 48502 17797") // centre of the 1 m square

// Conic grids: CGCS2000 Albers (equal area) with parallels 25°N/47°N, 105°E
let albers = @conic.ConicGrid::albers(
  6378137.0,
  1.0 / 298.257222101,
  lat_origin=0.0,
  lon_origin=105.0,
  stdlat1=25.0,
  stdlat2=47.0,
)
let gp = albers.forward(39.9087, 116.3975) // gp.easting = 956117.8 m, gp.northing, gp.gamma, gp.k

// Magnetic declination (World Magnetic Model 2025) for a compass
let wmm = @magnetic.MagneticModel::wmm2025()
let t = @magnetic.decimal_year(2026, 11, 10)
let mag = wmm.field(39.9087, 116.3975, 0.0, t) // mag.decl = -7.56° (west), mag.incl, mag.h ...
let compass = @magnetic.magnetic_bearing(r.azi1, mag.decl) // true bearing → compass bearing
```

A runnable tour lives in [`src/examples/demo`](src/examples/demo/main.mbt):

```bash
moon run src/examples/demo
```

A second example is a small command-line interface in [`src/examples/cli`](src/examples/cli/main.mbt). `geod` solves the WGS-84 geodesic inverse and direct problems. `rhumb` does the same for constant-azimuth rhumb lines. `intersect closest` finds the nearest intersection of two geodesics. `enu` converts between geodetic and local East-North-Up. `coord` converts one point among WGS-84, GCJ-02, and BD-09 (`gcj2wgs`, `bd2gcj`, and `bd2wgs` are the exact inverses). `mag` gives the WMM2025 declination and field for a place and date. `lcc` and `albers` project to and from conic grids.

```bash
moon run src/examples/cli -- geod inverse 39.9087 116.3975 31.2304 121.4737
moon run src/examples/cli -- geod direct 22.5431 114.0579 90 1000000
moon run src/examples/cli -- rhumb inverse 51.4778 -0.0015 35.6762 139.6503
moon run src/examples/cli -- intersect closest 0 0 45 45 10 135
moon run src/examples/cli -- enu forward 39.9087 116.3975 39.9097 116.3985 50
moon run src/examples/cli -- coord wgs2gcj 39.9087 116.3975
moon run src/examples/cli -- coord gcj2wgs 39.91010 116.40374
moon run src/examples/cli -- coord wgs2bd 39.9087 116.3975
moon run src/examples/cli -- coord bd2wgs 39.91645 116.41012
moon run src/examples/cli -- mag 39.9087 116.3975 2026-11-10
moon run src/examples/cli -- albers forward 25 47 0 105 39.9087 116.3975
moon run src/examples/cli -- lcc reverse 25 47 0 105 1000000 4000000
```

`coord` also accepts `gcj2bd` and `bd2gcj`. `rhumb direct`, `enu reverse`, and the rest are listed when you run `moon run src/examples/cli` with no arguments.

## Scenarios

Six end-to-end scenarios. CI runs each on all four backends, and each one aborts if a result differs from what this README says:

| Scenario | Run | What it shows |
| --- | --- | --- |
| [Map APIs → true distance](src/examples/s1_cn_distance/main.mbt) | `moon run src/examples/s1_cn_distance` | a Gaode (GCJ-02) shop and a Baidu (BD-09) customer converted to WGS-84, then the geodesic distance; mixing raw coordinates is off by ~764 m |
| [UTM ↔ MGRS field survey](src/examples/s2_utm_mgrs/main.mbt) | `moon run src/examples/s2_utm_mgrs` | UTM/UPS and MGRS for Beijing, Sydney, the Norway 32V and Svalbard 33X exceptions and the South Pole; UTM round trips below 1 µm, MGRS cell corners within √2 m |
| [Accuracy vs GeographicLib](src/examples/s3_accuracy/main.mbt) | `moon run src/examples/s3_accuracy` | 200 cases of the official GeodTest set, max inverse/direct errors per distance band (table below) |
| [Field and region area](src/examples/s4_polygon_area/main.mbt) | `moon run src/examples/s4_polygon_area` | a field traced in GCJ-02 converted to WGS-84, its geodesic area in m², ha and 亩, incremental `PolygonArea`, and a continental triangle whose area and perimeter agree with Python geographiclib to 1e-12 (relative) |
| [Compass bearing for a survey](src/examples/s5_compass/main.mbt) | `moon run src/examples/s5_compass` | GPS start, Gaode target: geodesic true bearing → compass bearing with the WMM2025 declination (263 m off over 2 km if ignored); the same bearing from a UTM map via grivation; declination from Mohe (−12.8°) to Kashgar (+4.1°) |
| [Areas on conic grids](src/examples/s6_conic_area/main.mbt) | `moon run src/examples/s6_conic_area` | 0.5° cells from Sanya to Mohe measured on CGCS2000 Albers (agrees with the ellipsoid to 1e-8; measured 1e-9), Lambert (off by up to 6%) and Web Mercator (off by up to 176%); Lambert scale and convergence by latitude |

## API reference

Signatures, units, and edge cases are in the package notes. `moon check` and `moon test` run the short example in each section.

- [`geodesic`](src/geodesic/README.mbt.md) — ellipsoids, inverse and direct, lines, polygons, rhumb lines, intersections, LocalCartesian ENU
- [`utm`](src/utm/README.mbt.md) — transverse Mercator, UTM/UPS, MGRS
- [`conic`](src/conic/README.mbt.md) — Lambert conformal conic, Albers equal-area, EPSG grids
- [`magnetic`](src/magnetic/README.mbt.md) — World Magnetic Model 2025, bearings, grivation
- [`cn`](src/cn/README.mbt.md) — WGS-84, GCJ-02, BD-09

Preview the generated symbol pages (`///` comments plus these notes):

```bash
moon doc --serve
```

## Accuracy and testing

Everything is checked against reference data that does not come from the
code we port, and `moon test --target all` runs every test on all four
backends. Each bound below is the one the test asserts; the measured maximum
is in brackets.

| Area | Reference | Asserted bound (measured) |
| --- | --- | --- |
| Geodesic inverse | 1,000 cases from Karney's [GeodTest-short.dat](https://doi.org/10.5281/zenodo.32156) (CC0) | distance < 8 nm (7.5 nm) |
| Geodesic direct | same data set | position < 10 nm (8.6–9.2 nm, by backend) |
| Geodesic inverse/direct | 2,020 cases vs Python `geographiclib` 2.1 | relative difference < 1e-11 |
| GeodesicExact | 111 cases vs GeographicLib 2.5 `GeodSolve -E`, oblate `f` up to 0.5, prolate down to −0.2 | distance < 8 nm (7.5 nm); position < 1e-13°, about 11 nm (6.3 nm) |
| Polygon area | 8 polygons (incl. pole-enclosing, antimeridian-crossing) vs `geographiclib` | < 1e-14 of the Earth's area |
| Rhumb inverse | 125 cases vs GeographicLib RhumbSolve 2.5 | distance < 4 nm (3.7 nm) |
| Rhumb direct | 296 cases vs the same RhumbSolve | position < 11 nm (10.5 nm); area < 5e-15 of the ellipsoid |
| Geodesic intersection | 139 cases (closest, offset, next, segment) vs GeographicLib IntersectTool 2.5; 2 ties compared by distance | \|Δx\|, \|Δy\| < 0.5 µm (0.45 µm); the curves meet to < 10 nm (7.9 nm) |
| Geocentric / LocalCartesian | 660 ECEF and 120 ENU cases (WGS-84, sphere, prolate, GRS80) vs 50-digit textbook formulas; 495 + 90 vs PROJ `cart` / `topocentric` | < 20 nm (7.5 nm) |
| Transverse Mercator / UTM / UPS | 600 points vs PROJ 9 (pyproj), up to ~930 km from the central meridian | < 10 nm (5.5 nm) |
| MGRS | 315 points vs NGA GeoTrans (`mgrs` package), precisions 0–5 | 315/315 identical strings |
| Lambert / Albers conic | 1592 points on 16 projections vs 50-digit evaluation of the EPSG GN 7-2 / Snyder formulas; PROJ 9.5; IOGP GIGS 5102/5103/5109; EPSG and Snyder worked examples | < 25 nm (19 nm); PROJ < 0.25 µm; GIGS within the files' rounding |
| WMM2025 | NOAA's 100 official double-precision test points; report Tables 3b, 3d, 6; 119 points vs an independent 40-digit evaluation | Y, Z < 2.5e-6 nT; X, H, F < 1e-3 nT (7.2e-4, see [notes](src/magnetic/README.mbt.md)); worked example to 10 decimals; independent < 1e-10 nT |
| GCJ-02 / BD-09 | 200 points vs `coord-convert`: GCJ-02 on the 141 inside China (the 59 on foreign land are left unchanged by design), BD-09 on all 200 | < 1e-11° |

Every reference data set has a generator in [`tools/`](tools/) (Python;
each script names its source, and the WMM, GeodTest-for-s3 and China-border
scripts check the SHA-256 of the downloaded files). The checklist that maps each README statement to its test is
[`docs/十一月验收清单.md`](docs/十一月验收清单.md).

Output of [`s3_accuracy`](src/examples/s3_accuracy/main.mbt) (200 cases of
GeodTest-short.dat, errors in nanometres; the example aborts if any is
above 8 nm):

| distance band | cases | max inverse \|Δs12\| (nm) | max direct position error (nm) |
|---|---|---|---|
| < 1 000 km | 21 | 0.96 | 1.77 |
| 1 000–10 000 km | 55 | 3.73 | 3.53 |
| > 10 000 km | 124 | 7.45 | 6.5 |
| all | 200 | 7.45 | 6.5 |

### Notes

- GCJ-02 is the published reverse-engineered algorithm used by most open-source
  libraries (Krasovsky ellipsoid, a = 6378245 m). The transform is the
  identity outside the usual China bounding box and on the land of other
  countries inside it (Seoul, Hanoi, ...), using Natural Earth 1:10m borders
  (China point of view). Open-source implementations are reported to agree
  with the official, unpublished cipher to about 1 m; this library matches
  them to 1e-11° but cannot check the official one.
- The series geodesic (`Geodesic`) assumes `|f|` up to about 1/50 and the
  rhumb series about 1/150 (GeographicLib's guidance). `GeodesicExact` is the
  elliptic-integral solution for larger flattening, including prolate
  ellipsoids.
- Nearly antipodal points may have more than one shortest geodesic; then the
  distance is unique but azimuths/area depend on which one is returned (same
  behaviour as GeographicLib).
- A rhumb line that crosses a pole has no single longitude; `lon2` and
  `area12` are NaN in that case.
- Two closest geodesic intersections can be the same L1 distance apart.
  `closest`, `next`, and `segment` then return one of them. `segment`
  needs one shortest geodesic between each pair of endpoints.
- Karney (2011) shows that the transverse Mercator series is accurate to a few
  nanometres within ~3900 km of the central meridian; the tests here reach
  ~930 km.
- WMM2025 is valid from 2025.0 to 2030.0 and from 1 km below to 850 km above
  the ellipsoid; `field` raises outside that range. The tests check that the
  code evaluates the model correctly; the model's own declination
  uncertainty is about 0.3°–0.4° away from the magnetic poles.

## Benchmarks

`moon bench -p ChrisTu117/moongeodesy/benchmarks --target native` (source:
[`src/benchmarks`](src/benchmarks/bench_test.mbt)). Mean time per call on one
core of a 4-core AMD EPYC cloud VM, moon 0.1.20261104, measured once on
2026-11-10. These numbers are indicative; CI only runs the benchmarks as a
smoke test.

| Operation | native | js |
| --- | ---: | ---: |
| geodesic inverse, Beijing → Shanghai | 1.37 µs | 1.86 µs |
| geodesic inverse, nearly antipodal | 1.44 µs | 2.05 µs |
| geodesic direct, 1 000 km | 0.45 µs | 1.01 µs |
| polygon area, 4 vertices | 2.99 µs | 5.16 µs |
| WGS-84 → GCJ-02 (with the border test) | 0.26 µs | 0.49 µs |
| GCJ-02 → WGS-84 (exact, iterative) | 0.63 µs | 0.91 µs |
| UTM from lat/lon | 0.64 µs | 0.52 µs |
| MGRS encode (1 m) | 0.78 µs | 0.91 µs |
| WMM2025 field, all elements and rates | 1.88 µs | 2.96 µs |
| Albers forward | 0.24 µs | 0.25 µs |
| Lambert forward | 0.27 µs | 0.24 µs |

## Port source and license

- `geodesic` (geodesics, `GeodesicExact`, polygons, rhumb lines,
  intersections, Geocentric and LocalCartesian), `utm`, and `conic` are
  ports of the algorithms in
  [GeographicLib](https://github.com/geographiclib/geographiclib) by Charles
  F. F. Karney (C++ and Python versions, v2.x), released under the
  **MIT/X11 License**, © Charles Karney. Rhumb lines follow GeographicLib's
  `Rhumb` class (order-6 series; Karney 2024,
  [doi:10.1007/s11200-024-0709-z](https://doi.org/10.1007/s11200-024-0709-z)).
  `GeodesicExact` follows GeographicLib's elliptic-integral geodesic
  (Karney 2012–2023; area via the DST in GeographicLib 2.5).
  Intersections follow GeographicLib's `Intersect` class (Karney 2024,
  [doi:10.1061/JSUED2.SUENG-1483](https://doi.org/10.1061/JSUED2.SUENG-1483);
  preprint [arXiv:2308.00495](https://arxiv.org/abs/2308.00495)).
  `conic` follows `LambertConformalConic` and `AlbersEqualArea`.
  Structure, series coefficients and edge-case handling follow GeographicLib;
  the MoonBit code is a rewrite, not a binding.
- `magnetic` is written from the equations in the WMM2025 technical report
  (Chulliat et al. 2025, [doi:10.25923/prbc-s316](https://doi.org/10.25923/prbc-s316)),
  not ported from a library. The coefficients and test values are from NOAA
  NCEI and the British Geological Survey
  ([doi:10.25921/aqfd-sd83](https://doi.org/10.25921/aqfd-sd83)); they are
  U.S. Government material, not subject to copyright protection
  (17 U.S.C. 403).
- GCJ-02/BD-09 use the widely published reverse-engineered formulas (as in
  `coord-convert`, MIT); the iterative exact inverses are our own.
- MGRS letter tables follow NGA's GeoTrans / GeographicLib `MGRS`.
- The China border data is Natural Earth 1:10m (public domain).
- Test data: GeodTest (Karney, CC0); IOGP GIGS test dataset 2.1 (© IOGP,
  used under the [GIGS Terms of Use](https://gigs.iogp.org/TermsOfUse), in
  `src/conic/gigs_data_test.mbt` only); the other vectors are generated by
  `tools/` with `geographiclib`, PROJ (`pyproj`), `mgrs`, `coord-convert`,
  GeographicLib command-line tools, and `mpmath`.
- This library is MIT licensed (see [LICENSE](LICENSE)).

## Code size

| Part | Lines of MoonBit |
| --- | --- |
| Library code (`geodesic`, `utm`, `conic`, `magnetic`, `cn`, `internal`) | ~9 200 |
| Generated data inside the library (China borders, WMM2025 coefficients) | ~1 600 |
| Hand-written tests | ~3 200 |
| Generated reference data (`*_data_test.mbt`, `vectors_test.mbt`) | ~38 000 |
| Package notes with checked examples (`README.mbt.md`) | ~860 |
| Examples (demo, CLI, six scenarios) | ~1 800 |

Most of the total is generated test vectors. Of the library code,
`GeodesicExact` (elliptic integrals, the line, and the sine transform for
area) is about 1.8k lines, `conic` about 1.4k and `magnetic` about 0.5k.

## History

### November 2026 contest period

- Fixes from the October acceptance review ([#32](https://github.com/ChrisTu117/moongeodesy/pull/32)):
  foreign land inside the China box is no longer offset, the exact GCJ-02
  inverse is exact next to the box edges and borders, UPS works at ±90°.
- World Magnetic Model 2025, package `magnetic`, with CLI `mag` and scenario
  s5 ([#33](https://github.com/ChrisTu117/moongeodesy/issues/33), [#34](https://github.com/ChrisTu117/moongeodesy/pull/34)).
- Lambert conformal conic and Albers equal-area projections, package
  `conic`, with CLI `lcc` / `albers` and scenario s6
  ([#35](https://github.com/ChrisTu117/moongeodesy/issues/35), [#36](https://github.com/ChrisTu117/moongeodesy/pull/36)).
- Every README statement mapped to a test or check
  ([`docs/十一月验收清单.md`](docs/十一月验收清单.md)): independent
  reference data for Geocentric/LocalCartesian, assertions tightened to the
  published bounds, scenarios that abort on wrong results and run on every
  backend in CI, and tests for API promises that had none.

### October 2026 contest period

The repository was created on 2026-10-09 for the MoonBit October hackathon.
In order: project scaffold and `cn` conversions → GeographicLib math helpers
→ geodesic inverse/direct/lines/polygons → transverse Mercator, UTM/UPS,
MGRS → GRS80/CGCS2000 ellipsoids → CLI example (#1) → API docs (#2) →
runnable scenarios (#3, #5) → this README (#4, #6) → benchmarks (#8) →
polygon area scenario (#10) → contest-list and PROGRESS refresh (#12) →
release 0.1.1 on [mooncakes.io](https://mooncakes.io/docs/ChrisTu117/moongeodesy)
(#14) → rhumb lines (#15) → geodesic intersections (#19) → LocalCartesian
ENU (#23) → CLI/demo coverage for rhumb, intersections, and ENU (#27) →
`GeodesicExact` for `|f| > 0.02` (#31). See the
[commit history](https://github.com/ChrisTu117/moongeodesy/commits/main)
and [pull requests](https://github.com/ChrisTu117/moongeodesy/pulls?q=is%3Apr).

## AI assistance

- **Goal.** A dependency-free MoonBit library for the geodesy that GIS, maps
  and logistics code needs in China: ellipsoidal distance/area, UTM/MGRS and
  conic grids, compass declination, and the GCJ-02/BD-09 offsets used by
  Gaode, Tencent and Baidu.
- **Approach.** An AI coding assistant translated GeographicLib's algorithms
  function by function, wrote the WMM2025 code from the NOAA report, and
  wrote tests; the author chose scope and package boundaries, reviewed the
  numerical code against the papers and sources, and decided on API shape
  and error handling. [`docs/十一月新增功能说明.md`](docs/十一月新增功能说明.md)
  explains the November code in Chinese.
- **How correctness is verified.** Not by reading the AI's output, but by
  reference data that does not share the ported code: Karney's official
  GeodTest data, Python `geographiclib`, PROJ, NGA GeoTrans,
  `coord-convert`, NOAA's WMM test values, IOGP GIGS, EPSG and Snyder worked
  examples, and 40–50-digit `mpmath` evaluations of textbook formulas
  (table above). Every README statement has a test or check, listed in the
  acceptance checklist; CI runs everything on wasm, wasm-gc, js and native on
  Linux and macOS. The geodesic code needs a correctly rounded `hypot` on very
  short lines; a named test checks two such lines against geographiclib.

## Related MoonBit libraries

moongeodesy covers the ellipsoid and grid systems; these libraries on
[mooncakes.io](https://mooncakes.io) work in the plane or on file formats and
pair well with it:

- [`totto2727/geo`](https://mooncakes.io/docs/totto2727/geo) — 2D planar geometry (port of georust/geo)
- [`cn-xjr/moongeokit`](https://mooncakes.io/docs/cn-xjr/moongeokit) — robust 2D geometry and spatial indexing
- [`LING71671/moon_strtree`](https://mooncakes.io/docs/LING71671/moon_strtree) — STR-packed spatial indexes
- [`NBB2006/routeweave`](https://mooncakes.io/docs/NBB2006/routeweave) — route encoding, GeoJSON, Web Mercator tiles
- [`2973181701/moongpx`](https://mooncakes.io/docs/2973181701/moongpx) — GPX parsing
- [`liyun6666/moon-shapefile`](https://mooncakes.io/docs/liyun6666/moon-shapefile) — Shapefile reader/writer

## Credits

- Geodesic and projection algorithms: C. F. F. Karney, *Algorithms for
  geodesics*, J. Geodesy 87, 43–55 (2013), *Transverse Mercator with an
  accuracy of a few nanometers*, J. Geodesy 85, 475–485 (2011), and
  *Geodesic intersections*, J. Surveying Eng. 150(3), 04024005 (2024);
  code structure follows GeographicLib (MIT/X11 License, © Charles Karney).
- Conic projection formulas: IOGP Guidance Note 7-2 (EPSG methods 9801,
  9802, 9822); J. P. Snyder, *Map Projections — A Working Manual*, USGS
  Professional Paper 1395 (1987).
- Geomagnetic model: A. Chulliat et al., *The US/UK World Magnetic Model for
  2025–2030: Technical Report*, NOAA NCEI (2025).
- Test data: GeodTest.dat, Karney (2010), doi:10.5281/zenodo.32156, CC0;
  WMM2025 test values, NOAA NCEI; GIGS test dataset 2.1, © IOGP.

## License

MIT

---

## 中文简介

**moongeodesy** 是一个纯 MoonBit 实现的大地测量库，无 FFI，支持 wasm / wasm-gc / js / native 全部后端：

- **`geodesic`**：移植 GeographicLib 的 Karney 测地线算法——椭球面上两点的最短距离与方位角（反算）、给定方位与距离求终点（正算）、测地线插点（用于地图画航线）、任意多边形的椭球面面积与周长、恒向线（rhumb line）、两条测地线的交点，以及相对原点的局部东北天直角坐标（`LocalCartesian`，经地心地固 ECEF）。级数解用于 `|f|` 约 1/50 以内（WGS-84 等）；`|f|` 更大时用 `GeodesicExact`（椭圆积分）。官方 GeodTest 测试集上反算距离误差 < 8 nm、正算位置误差 < 10 nm；恒向线与 RhumbSolve 2.5 对照 < 4 nm；交点与 IntersectTool 2.5 对照 < 0.5 µm；ECEF 与局部 ENU 与 50 位精度的教科书公式及 PROJ 对照 < 20 nm。
- **`utm`**：横轴墨卡托投影（Krüger 6 阶级数）、极球面投影、UTM/UPS 坐标（含挪威、斯瓦尔巴特特殊分带）以及 MGRS 军事格网参考系的编码与解码。与 PROJ 结果相差 < 10 nm，MGRS 字符串与 NGA GeoTrans 完全一致。
- **`conic`**（11 月新增）：兰勃特等角圆锥投影和阿尔伯斯等积圆锥投影，单/双标准纬线，任意椭球，带 EPSG 假原点的格网（如 CGCS2000、25°N/47°N、中央经线 105°E 的全国 Albers 等积投影，用于面积统计）。与 50 位精度的 EPSG/Snyder 公式对照 < 25 nm，与 PROJ 对照 < 0.25 µm，并通过 IOGP GIGS 官方测试。
- **`magnetic`**（11 月新增）：世界地磁模型 WMM2025，给出任意地点和日期的磁偏角、磁倾角、磁场分量及年变化率，罗盘方位与真方位、格网方位互换，以及磁极附近罗盘不可用的“禁用区/警告区”判断。与 NOAA 官方 100 个测试点、技术报告逐步算例对照，并与一个独立的 40 位精度实现对照 < 1e-10 nT。
- **`cn`**：国内地图坐标系转换——WGS-84（GPS）↔ GCJ-02（高德、腾讯，"火星坐标"）↔ BD-09（百度），提供迭代精确逆转换（往返误差 < 1e-9°）以及常用的单步近似公式。国外（含首尔、河内等框内国家）坐标不做偏移。

```moonbit
let g = @geodesic.Geodesic::wgs84()
let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737) // 北京 → 上海：r.s12 ≈ 1066.68 km
let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → 高德
let bd = @cn.gcj02_to_bd09(gcj) // 高德 → 百度
let mgrs = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
let rh = @geodesic.Rhumb::wgs84()
let h = rh.inverse(51.4778, -0.0015, 35.6762, 139.6503) // 伦敦 → 东京，恒向线
let ix = @geodesic.Intersect::wgs84()
let p = ix.closest(0.0, 0.0, 45.0, 45.0, 10.0, 135.0) // 两条测地线的最近交点
let lc = @geodesic.LocalCartesian::wgs84(39.9087, 116.3975) // 以天安门为原点的局部 ENU
let enu = lc.forward(39.9097, 116.3985, 50.0) // enu.x 东, enu.y 北, enu.z 天
let wmm = @magnetic.MagneticModel::wmm2025()
let mag = wmm.field(39.9087, 116.3975, 0.0, 2026.86) // 北京磁偏角 mag.decl ≈ -7.56°（偏西）
```

所有结果都与不共享移植代码的参考数据交叉验证（Karney 官方测试数据、geographiclib、PROJ、GeoTrans、coord-convert、NOAA 官方测试值、IOGP GIGS、EPSG 与 Snyder 算例、mpmath 高精度计算）；参考数据生成脚本见 `tools/`。README 里每条说法对应的测试见 [`docs/十一月验收清单.md`](docs/十一月验收清单.md)，11 月新增代码的中文说明见 [`docs/十一月新增功能说明.md`](docs/十一月新增功能说明.md)。

命令行示例在 `src/examples/cli`：

```bash
moon run src/examples/cli -- geod inverse 39.9087 116.3975 31.2304 121.4737
moon run src/examples/cli -- rhumb inverse 51.4778 -0.0015 35.6762 139.6503
moon run src/examples/cli -- intersect closest 0 0 45 45 10 135
moon run src/examples/cli -- enu forward 39.9087 116.3975 39.9097 116.3985 50
moon run src/examples/cli -- coord wgs2gcj 39.9087 116.3975
moon run src/examples/cli -- mag 39.9087 116.3975 2026-11-10
moon run src/examples/cli -- albers forward 25 47 0 105 39.9087 116.3975
```

`geod` 做测地线正反算，`rhumb` 做恒向线正反算，`intersect closest` 求两条测地线最近交点，`enu` 做局部东北天与经纬度互转。`coord` 支持 `wgs2gcj`、`gcj2wgs`、`gcj2bd`、`bd2gcj`、`wgs2bd`、`bd2wgs`。`mag` 计算某地某日的磁偏角，`lcc` 和 `albers` 做圆锥投影正反算。无参数时打印命令列表。

接口说明（单位、边界、短示例）在各包的 `README.mbt.md`。本地预览：`moon doc --serve`。
