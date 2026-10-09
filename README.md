# moongeodesy

Pure MoonBit geodesy: nanometre-accurate geodesics on the ellipsoid
(a port of Charles Karney's [GeographicLib](https://geographiclib.sourceforge.io/)
algorithms), UTM/UPS and MGRS grid references, and the WGS-84 ↔ GCJ-02 ↔ BD-09
conversions needed for maps in mainland China. No FFI; works on the `wasm`,
`wasm-gc`, `js` and `native` backends.

> 中文简介见[下方](#中文简介)。

| Package | What it does |
| --- | --- |
| `ChrisTu117/moongeodesy/geodesic` | Karney's geodesic algorithms: inverse & direct problems, geodesic lines, polygon area/perimeter on any ellipsoid |
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

// Other ellipsoids: GRS80, a sphere, Mars...
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

## Accuracy and testing

Everything is checked against independent reference implementations and
data sets (`moon test --target all`, 25 tests on 4 backends):

| Area | Reference | Result |
| --- | --- | --- |
| Geodesic inverse | 1,000 cases from Karney's [GeodTest-short.dat](https://doi.org/10.5281/zenodo.32156) (CC0) | max distance error 7.5 nm |
| Geodesic direct | same data set | max position error 8.6 nm |
| Geodesic inverse/direct | 2,020 cases vs Python `geographiclib` 2.1 | relative difference < 1e-11 |
| Polygon area | 8 polygons (incl. pole-enclosing, antimeridian-crossing) vs `geographiclib` | < 1e-14 of the Earth's area |
| Transverse Mercator / UTM / UPS | 600 points vs PROJ 9 (pyproj) | max difference 5.5 nm |
| MGRS | 315 points vs NGA GeoTrans (`mgrs` package) | 315/315 identical strings |
| GCJ-02 / BD-09 | 200 points vs `coord-convert` | < 1e-11° |

The generators in [`tools/`](tools/) regenerate all reference vectors.

### Notes

- GCJ-02 is the published reverse-engineered algorithm used by most open-source
  libraries (Krasovsky ellipsoid, a = 6378245 m); outside the rough mainland
  China bounding box the transform is the identity. Accuracy of the forward
  transform relative to the official one is ~1 m.
- Nearly antipodal points may have more than one shortest geodesic; then the
  distance is unique but azimuths/area depend on which one is returned (same
  behaviour as GeographicLib).
- The transverse Mercator series is accurate to a few nanometres within
  ~3900 km of the central meridian.

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

- **`geodesic`**：移植 GeographicLib 的 Karney 测地线算法——椭球面上两点的最短距离与方位角（反算）、给定方位与距离求终点（正算）、测地线插点（用于地图画航线）、任意多边形的椭球面面积与周长。精度为纳米级（官方 GeodTest 测试集上最大误差 < 10 nm）。
- **`utm`**：横轴墨卡托投影（Krüger 6 阶级数）、极球面投影、UTM/UPS 坐标（含挪威、斯瓦尔巴特特殊分带）以及 MGRS 军事格网参考系的编码与解码。与 PROJ 结果相差 < 10 nm，MGRS 字符串与 NGA GeoTrans 完全一致。
- **`cn`**：国内地图坐标系转换——WGS-84（GPS）↔ GCJ-02（高德、腾讯，"火星坐标"）↔ BD-09（百度），提供迭代精确逆转换（往返误差 < 1e-9°）以及常用的单步近似公式。

```moonbit
let g = @geodesic.Geodesic::wgs84()
let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737) // 北京 → 上海：r.s12 ≈ 1066.68 km
let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → 高德
let bd = @cn.gcj02_to_bd09(gcj) // 高德 → 百度
let mgrs = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
```

所有结果均与独立参考实现（GeographicLib Python、PROJ、GeoTrans、coord-convert）及 Karney 官方测试数据交叉验证；参考数据生成脚本见 `tools/`。
