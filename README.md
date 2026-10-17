# moongeodesy

[![CI](https://github.com/ChrisTu117/moongeodesy/actions/workflows/ci.yml/badge.svg)](https://github.com/ChrisTu117/moongeodesy/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-mooncakes.io-green)](https://mooncakes.io/docs/ChrisTu117/moongeodesy)

Pure MoonBit geodesy: nanometre-accurate geodesics on the ellipsoid
(a port of Charles Karney's [GeographicLib](https://geographiclib.sourceforge.io/)
algorithms), UTM/UPS and MGRS grid references, and the WGS-84 ↔ GCJ-02 ↔ BD-09
conversions needed for maps in mainland China. No FFI; works on the `wasm`,
`wasm-gc`, `js` and `native` backends.

> 中文简介见[下方](#中文简介)。

| Package | What it does |
| --- | --- |
| `ChrisTu117/moongeodesy/geodesic` | Karney's geodesic algorithms and rhumb lines: inverse & direct problems, geodesic lines, polygon area/perimeter on any ellipsoid |
| `ChrisTu117/moongeodesy/utm` | Transverse Mercator (6th-order Krüger series), polar stereographic, UTM/UPS with standard zones, MGRS encode/decode |
| `ChrisTu117/moongeodesy/cn` | WGS-84 → GCJ-02 (Gaode/Tencent), GCJ-02 → BD-09 (Baidu), with exact iterative inverses |

## Install

```bash
moon add ChrisTu117/moongeodesy
```

and import the packages you need in your `moon.pkg`:

```
import {
  "ChrisTu117/moongeodesy/geodesic",
  "ChrisTu117/moongeodesy/utm",
  "ChrisTu117/moongeodesy/cn",
}
```

## Usage

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

// Other ellipsoids: CGCS2000 / GRS80 built in, or any (a, f)
let cgcs = @geodesic.Geodesic::cgcs2000()
let sphere = @geodesic.Geodesic::new(6371000.0, 0.0)

// Chinese map coordinate systems
let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → Gaode
let bd = @cn.gcj02_to_bd09(gcj) // Gaode → Baidu
let wgs = @cn.bd09_to_wgs84(bd) // exact inverse (round trip < 1e-9°)

// UTM / UPS / MGRS (functions raise @utm.UtmError on bad input)
let u = @utm.from_latlon(39.9087, 116.3975) // zone 50N, easting/northing, convergence, scale
let m = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
let p = @utm.from_mgrs("50S MK 48502 17797") // centre of the 1 m square
```

A runnable tour lives in [`src/examples/demo`](src/examples/demo/main.mbt):

```bash
moon run src/examples/demo
```

A second example is a small command-line interface in [`src/examples/cli`](src/examples/cli/main.mbt). `geod` solves the WGS-84 inverse and direct problems. `coord` converts one point among WGS-84, GCJ-02, and BD-09 (`gcj2wgs`, `bd2gcj`, and `bd2wgs` are the exact inverses).

```bash
moon run src/examples/cli -- geod inverse 39.9087 116.3975 31.2304 121.4737
moon run src/examples/cli -- geod direct 22.5431 114.0579 90 1000000
moon run src/examples/cli -- coord wgs2gcj 39.9087 116.3975
moon run src/examples/cli -- coord gcj2wgs 39.91010 116.40374
moon run src/examples/cli -- coord wgs2bd 39.9087 116.3975
moon run src/examples/cli -- coord bd2wgs 39.91645 116.41012
```

`coord` also accepts `gcj2bd` and `bd2gcj`. Run `moon run src/examples/cli` to print the full command list.

## Scenarios

Four end-to-end scenarios, each runnable on every backend and run in CI:

| Scenario | Run | What it shows |
| --- | --- | --- |
| [Map APIs → true distance](src/examples/s1_cn_distance/main.mbt) | `moon run src/examples/s1_cn_distance` | a Gaode (GCJ-02) shop and a Baidu (BD-09) customer converted to WGS-84, then the geodesic distance; mixing raw coordinates is off by ~764 m |
| [UTM ↔ MGRS field survey](src/examples/s2_utm_mgrs/main.mbt) | `moon run src/examples/s2_utm_mgrs` | UTM/UPS and MGRS for Beijing, Sydney, the Norway 32V and Svalbard 33X exceptions and the South Pole, with round-trip residuals |
| [Accuracy vs GeographicLib](src/examples/s3_accuracy/main.mbt) | `moon run src/examples/s3_accuracy` | 200 cases of the official GeodTest set, max inverse/direct errors per distance band (table below) |
| [Field and region area](src/examples/s4_polygon_area/main.mbt) | `moon run src/examples/s4_polygon_area` | a field traced in GCJ-02 converted to WGS-84, its geodesic area in m², ha and 亩, incremental `PolygonArea`, and a continental triangle that matches Python geographiclib exactly (the example aborts if it does not) |

## API reference

Signatures, units, and edge cases are in the package notes. `moon check` and `moon test` run the short example in each section.

- [`geodesic`](src/geodesic/README.mbt.md) — ellipsoids, inverse and direct, lines, polygons, rhumb lines
- [`utm`](src/utm/README.mbt.md) — transverse Mercator, UTM/UPS, MGRS
- [`cn`](src/cn/README.mbt.md) — WGS-84, GCJ-02, BD-09

Preview the generated symbol pages (`///` comments plus these notes):

```bash
moon doc --serve
```

## Accuracy and testing

Everything is checked against independent reference implementations and
data sets (`moon test --target all`, 50 tests on 4 backends):

| Area | Reference | Result |
| --- | --- | --- |
| Geodesic inverse | 1,000 cases from Karney's [GeodTest-short.dat](https://doi.org/10.5281/zenodo.32156) (CC0) | max distance error 7.5 nm |
| Geodesic direct | same data set | max position error 8.6 nm |
| Geodesic inverse/direct | 2,020 cases vs Python `geographiclib` 2.1 | relative difference < 1e-11 |
| Polygon area | 8 polygons (incl. pole-enclosing, antimeridian-crossing) vs `geographiclib` | < 1e-14 of the Earth's area |
| Rhumb inverse | 125 cases vs GeographicLib RhumbSolve 2.5 | max distance error 3.73 nm |
| Rhumb direct | 296 cases vs the same RhumbSolve | max position error 10.5 nm |
| Rhumb area | those inverse and direct cases | < 4e-15 of the ellipsoid area |
| Transverse Mercator / UTM / UPS | 600 points vs PROJ 9 (pyproj) | max difference 5.5 nm |
| MGRS | 315 points vs NGA GeoTrans (`mgrs` package) | 315/315 identical strings |
| GCJ-02 / BD-09 | 200 points vs `coord-convert` | < 1e-11° |

The generators in [`tools/`](tools/) regenerate all reference vectors.

Output of [`s3_accuracy`](src/examples/s3_accuracy/main.mbt) (200 cases of
GeodTest-short.dat, errors in nanometres):

| distance band | cases | max inverse \|Δs12\| (nm) | max direct position error (nm) |
|---|---|---|---|
| < 1 000 km | 21 | 0.96 | 1.77 |
| 1 000–10 000 km | 55 | 3.73 | 3.53 |
| > 10 000 km | 124 | 7.45 | 6.5 |
| all | 200 | 7.45 | 6.5 |

### Notes

- GCJ-02 is the published reverse-engineered algorithm used by most open-source
  libraries (Krasovsky ellipsoid, a = 6378245 m); outside the rough mainland
  China bounding box the transform is the identity. Accuracy of the forward
  transform relative to the official one is ~1 m.
- Nearly antipodal points may have more than one shortest geodesic; then the
  distance is unique but azimuths/area depend on which one is returned (same
  behaviour as GeographicLib).
- A rhumb line that crosses a pole has no single longitude; `lon2` and
  `area12` are NaN in that case. The rhumb series assume `|f|` up to about
  1/150.
- The transverse Mercator series is accurate to a few nanometres within
  ~3900 km of the central meridian.

## Benchmarks

`moon bench -p ChrisTu117/moongeodesy/benchmarks --target native` (source:
[`src/benchmarks`](src/benchmarks/bench_test.mbt)). Mean time per call on one
x86-64 Linux core, moon 0.1.20260920; your numbers will differ.

| Operation | native | js |
| --- | ---: | ---: |
| geodesic inverse, Beijing → Shanghai | 3.0 µs | 3.6 µs |
| geodesic inverse, nearly antipodal | 3.0 µs | 3.0 µs |
| geodesic direct, 1 000 km | 0.66 µs | 2.3 µs |
| polygon area, 4 vertices | 4.6 µs | 7.8 µs |
| WGS-84 → GCJ-02 | 0.14 µs | 0.13 µs |
| GCJ-02 → WGS-84 (exact, iterative) | 0.80 µs | 0.62 µs |
| UTM from lat/lon | 0.86 µs | 0.83 µs |
| MGRS encode (1 m) | 1.5 µs | 2.3 µs |

So one core handles roughly 300 000 geodesic distances per second.

## Port source and license

- `geodesic` (geodesics, polygons, and rhumb lines) and `utm` are ports of
  the algorithms in
  [GeographicLib](https://github.com/geographiclib/geographiclib) by Charles
  F. F. Karney (C++ and Python versions, v2.x), released under the
  **MIT/X11 License**, © Charles Karney. Rhumb lines follow GeographicLib's
  `Rhumb` class (order-6 series; Karney 2024,
  [doi:10.1007/s11200-024-0709-z](https://doi.org/10.1007/s11200-024-0709-z)).
  Structure, series coefficients and edge-case handling follow GeographicLib;
  the MoonBit code is a rewrite, not a binding.
- GCJ-02/BD-09 use the widely published reverse-engineered formulas (as in
  `coord-convert`, MIT); the iterative exact inverses are our own.
- MGRS letter tables follow NGA's GeoTrans / GeographicLib `MGRS`.
- This library is MIT licensed (see [LICENSE](LICENSE)).

## Code size

| Part | Lines of MoonBit |
| --- | --- |
| Core library (`geodesic`, `utm`, `cn`, `internal`) | ~4 300 |
| Hand-written tests | ~730 |
| Generated reference data (`*_data_test.mbt`, `vectors_test.mbt`) | ~18 700 |
| Examples (demo, CLI, four scenarios) | ~1 400 |

The large total is almost entirely generated test vectors; the library itself
is about 4.3k lines.

## What is new in this contest period (October 2026)

The repository was created on 2026-10-09 for the MoonBit October hackathon;
all code is new in this period. In order: project scaffold and `cn`
conversions → GeographicLib math helpers → geodesic inverse/direct/lines/
polygons → transverse Mercator, UTM/UPS, MGRS → GRS80/CGCS2000 ellipsoids →
CLI example (#1) → API docs (#2) → runnable scenarios (#3, #5) → this README
(#4, #6) → benchmarks (#8) → polygon area scenario (#10) → contest-list and
PROGRESS refresh (#12) → release 0.1.1 on
[mooncakes.io](https://mooncakes.io/docs/ChrisTu117/moongeodesy) (#14) → rhumb
lines (#15). See the [commit history](https://github.com/ChrisTu117/moongeodesy/commits/main)
and [pull requests](https://github.com/ChrisTu117/moongeodesy/pulls?q=is%3Apr).

## AI assistance

- **Goal.** A dependency-free MoonBit library for the geodesy that GIS, maps
  and logistics code needs in China: ellipsoidal distance/area, UTM/MGRS grids
  and the GCJ-02/BD-09 offsets used by Gaode, Tencent and Baidu.
- **Approach.** An AI coding assistant translated GeographicLib's algorithms
  function by function and wrote tests; the author chose scope and package
  boundaries, reviewed the numerical code against Karney's papers and the
  GeographicLib sources, and decided on API shape and error handling.
- **How correctness is verified.** Not by reading the AI's output, but by
  independent references: Karney's official GeodTest data, Python
  `geographiclib`, PROJ 9, NGA GeoTrans and `coord-convert` (table above). The
  reference vectors are regenerated by scripts in `tools/`; CI runs every test
  on wasm, wasm-gc, js and native on Linux and macOS. Bugs found this way (e.g.
  the need for an accurate `hypot` on very short lines) were fixed and are
  covered by tests.

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
  geodesics*, J. Geodesy 87, 43–55 (2013) and *Transverse Mercator with an
  accuracy of a few nanometers*, J. Geodesy 85, 475–485 (2011); code structure
  follows GeographicLib (MIT/X11 License, © Charles Karney).
- Test data: GeodTest.dat, Karney (2010), doi:10.5281/zenodo.32156, CC0.

## License

MIT

---

## 中文简介

**moongeodesy** 是一个纯 MoonBit 实现的大地测量库，无 FFI，支持 wasm / wasm-gc / js / native 全部后端：

- **`geodesic`**：移植 GeographicLib 的 Karney 测地线算法——椭球面上两点的最短距离与方位角（反算）、给定方位与距离求终点（正算）、测地线插点（用于地图画航线）、任意多边形的椭球面面积与周长，以及恒向线（等方位角的 rhumb line：反算、正算、航线与赤道之间的面积）。测地线精度为纳米级（官方 GeodTest 测试集上最大误差 < 10 nm）；恒向线与 GeographicLib RhumbSolve 2.5 对照，距离最大误差 3.73 nm。
- **`utm`**：横轴墨卡托投影（Krüger 6 阶级数）、极球面投影、UTM/UPS 坐标（含挪威、斯瓦尔巴特特殊分带）以及 MGRS 军事格网参考系的编码与解码。与 PROJ 结果相差 < 10 nm，MGRS 字符串与 NGA GeoTrans 完全一致。
- **`cn`**：国内地图坐标系转换——WGS-84（GPS）↔ GCJ-02（高德、腾讯，"火星坐标"）↔ BD-09（百度），提供迭代精确逆转换（往返误差 < 1e-9°）以及常用的单步近似公式。

```moonbit
let g = @geodesic.Geodesic::wgs84()
let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737) // 北京 → 上海：r.s12 ≈ 1066.68 km
let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → 高德
let bd = @cn.gcj02_to_bd09(gcj) // 高德 → 百度
let mgrs = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
let rh = @geodesic.Rhumb::wgs84()
let h = rh.inverse(51.4778, -0.0015, 35.6762, 139.6503) // 伦敦 → 东京，恒向线
```

所有结果均与独立参考实现（GeographicLib Python、PROJ、GeoTrans、coord-convert）及 Karney 官方测试数据交叉验证；参考数据生成脚本见 `tools/`。

命令行示例在 `src/examples/cli`：

```bash
moon run src/examples/cli -- geod inverse 39.9087 116.3975 31.2304 121.4737
moon run src/examples/cli -- coord wgs2gcj 39.9087 116.3975
```

`geod inverse` 求 WGS-84 最短距离与方位角，`geod direct` 做正算。`coord` 支持 `wgs2gcj`、`gcj2wgs`、`gcj2bd`、`bd2gcj`、`wgs2bd`、`bd2wgs`。无参数时打印命令列表。

接口说明（单位、边界、短示例）在 `src/geodesic/README.mbt.md`、`src/utm/README.mbt.md`、`src/cn/README.mbt.md`。本地预览：`moon doc --serve`。
