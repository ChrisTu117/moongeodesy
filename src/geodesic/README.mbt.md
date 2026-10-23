# geodesic

椭球面测地线：正算、反算、航线插点、多边形面积、测地线交点、局部 ENU 直角坐标。

Package `ChrisTu117/moongeodesy/geodesic`. Angles are degrees. Lengths are metres. Areas are square metres. Azimuth is clockwise from north.

`internal/geomath` is an implementation detail. Import `geodesic` instead.

The command-line front end for inverse and direct is `src/examples/cli`.

## Ellipsoids

椭球。`wgs84`、`grs80`、`cgcs2000`，或任意 `(a, f)`。

```moonbit nocheck
pub let wgs84_a : Double
pub let wgs84_f : Double

pub fn Geodesic::wgs84() -> Geodesic
pub fn Geodesic::grs80() -> Geodesic
pub fn Geodesic::cgcs2000() -> Geodesic
pub fn Geodesic::new(a : Double, f : Double) -> Geodesic raise GeodesicError
pub fn Geodesic::equatorial_radius(Self) -> Double
pub fn Geodesic::flattening(Self) -> Double
pub fn Geodesic::ellipsoid_area(Self) -> Double
```

`a` is the equatorial radius in metres. `f` is the flattening. `new` raises `GeodesicError::InvalidEllipsoid` when `a` or the polar semi-axis `a * (1 - f)` is not finite and positive.

`cgcs2000()` uses the same `a` and `f` as `grs80()` (`a = 6378137` m, `f = 1/298.257222101`). CGCS2000 and GRS80 differ in GM and rotation rate, not in this ellipsoid.

The series assume `|f|` up to about 1/50. WGS-84, GRS80, and CGCS2000 are inside that range.

```mbt check
///|
test "ellipsoids" {
  let w = @geodesic.Geodesic::wgs84()
  assert_eq(w.equatorial_radius(), 6378137.0)
  let cgcs = @geodesic.Geodesic::cgcs2000()
  assert_eq(cgcs.flattening(), @geodesic.Geodesic::grs80().flattening())
  let sphere = @geodesic.Geodesic::new(6371000.0, 0.0)
  assert_eq(sphere.flattening(), 0.0)
  let bad = try {
    ignore(@geodesic.Geodesic::new(0.0, 0.0))
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
```

## Inverse and direct

反算与正算。纬度超出 ±90° 得到 NaN，不抛错。

```moonbit nocheck
pub fn Geodesic::inverse(
  Self, lat1 : Double, lon1 : Double, lat2 : Double, lon2 : Double,
  long_unroll? : Bool,
) -> InverseResult

pub fn Geodesic::distance(
  Self, lat1 : Double, lon1 : Double, lat2 : Double, lon2 : Double,
) -> Double

pub fn Geodesic::direct(
  Self, lat1 : Double, lon1 : Double, azi1 : Double, s12 : Double,
  long_unroll? : Bool,
) -> DirectResult

pub fn Geodesic::arc_direct(
  Self, lat1 : Double, lon1 : Double, azi1 : Double, a12 : Double,
  long_unroll? : Bool,
) -> DirectResult
```

`InverseResult` and `DirectResult` carry the endpoints, `s12` (m), `azi1` and `azi2` (degrees), auxiliary arc `a12` (degrees), reduced length `m12` (m), scales `scale12` and `scale21`, and equatorial area `area12` (m²).

`distance` returns only `s12`.

`arc_direct` takes the distance as an arc length in degrees on the auxiliary sphere.

Longitudes are reduced to [-180, 180]. Set `long_unroll=true` when the path crosses the antimeridian and `lon2` must stay continuous (`lon1` plus the signed difference).

Latitude ±90 is a valid pole. A latitude outside [-90, 90] is returned as NaN. No error is raised.

Nearly antipodal points can have more than one shortest geodesic. `s12` is unique. `azi1`, `azi2`, and `area12` follow one of those geodesics, as in GeographicLib.

```mbt check
///|
test "inverse and direct" {
  let g = @geodesic.Geodesic::wgs84()
  let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737)
  assert_true((r.s12 - 1066676.4).abs() < 1.0)
  let d = g.direct(22.5431, 114.0579, 90.0, 1000000.0)
  let back = g.inverse(22.5431, 114.0579, d.lat2, d.lon2)
  assert_true((back.s12 - 1000000.0).abs() < 1.0e-6)
  let pole = g.inverse(0.0, 0.0, 90.0, 0.0)
  assert_true((pole.s12 - 10001965.729).abs() < 1.0e-2)
  let bad = g.inverse(91.0, 0.0, 0.0, 0.0)
  assert_true(bad.lat1.is_nan())
}
```

## Geodesic lines

测地线上的点。先固定点 3，再取位置或加密。

```moonbit nocheck
pub fn Geodesic::line(Self, lat1 : Double, lon1 : Double, azi1 : Double) -> GeodesicLine
pub fn Geodesic::inverse_line(
  Self, lat1 : Double, lon1 : Double, lat2 : Double, lon2 : Double,
) -> GeodesicLine
pub fn Geodesic::direct_line(
  Self, lat1 : Double, lon1 : Double, azi1 : Double, s12 : Double,
) -> GeodesicLine

pub fn GeodesicLine::position(Self, s12 : Double, long_unroll? : Bool) -> DirectResult
pub fn GeodesicLine::arc_position(Self, a12 : Double, long_unroll? : Bool) -> DirectResult
pub fn GeodesicLine::distance(Self) -> Double
pub fn GeodesicLine::arc(Self) -> Double
pub fn GeodesicLine::set_distance(Self, s13 : Double) -> Unit
pub fn GeodesicLine::set_arc(Self, a13 : Double) -> Unit
pub fn GeodesicLine::densify(Self, max_step : Double) -> Array[(Double, Double)]
```

`inverse_line` and `direct_line` set point 3. `line` does not: `distance()` is NaN until `set_distance` or `set_arc`.

`position` takes metres from point 1. `arc_position` takes degrees of auxiliary arc. `densify` returns `(lat, lon)` points from point 1 to point 3, at most `max_step` metres apart. It needs a finite point 3.

```mbt check
///|
test "geodesic line" {
  let g = @geodesic.Geodesic::wgs84()
  let line = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
  let pts = line.densify(1000000.0)
  assert_true(pts.length() >= 2)
  let sep = g.distance(31.2304, 121.4737, 37.7749, -122.4194)
  assert_true((line.distance() - sep).abs() < 1.0e-6)
  let east = g.direct_line(0.0, 179.0, 90.0, 200000.0)
  let p = east.position(200000.0, long_unroll=true)
  assert_true(p.lon2 > 180.0)
}
```

## Polygons

多边形面积与周长。逆时针为正。可含极点，可跨日期变更线。

```moonbit nocheck
pub fn Geodesic::polygon(Self, polyline? : Bool) -> PolygonArea
pub fn Geodesic::polygon_area(Self, ArrayView[(Double, Double)]) -> PolygonResult

pub fn PolygonArea::add_point(Self, lat : Double, lon : Double) -> Unit
pub fn PolygonArea::add_edge(Self, azi : Double, s : Double) -> Unit
pub fn PolygonArea::clear(Self) -> Unit
pub fn PolygonArea::compute(Self, reverse? : Bool, sign? : Bool) -> PolygonResult
pub fn PolygonArea::test_point(
  Self, lat : Double, lon : Double, reverse? : Bool, sign? : Bool,
) -> PolygonResult
pub fn PolygonArea::test_edge(
  Self, azi : Double, s : Double, reverse? : Bool, sign? : Bool,
) -> PolygonResult
```

`PolygonResult` is `num` (vertex count), `perimeter` (m), and `area` (m²).

`polygon_area` closes the ring. Counter-clockwise order gives a positive area. `sign=true` (the default) keeps that sign. `reverse=true` makes clockwise positive. `sign=false` folds the area into one hemisphere of the ellipsoid.

`polygon(polyline=true)` sums length only. `area` is then NaN.

`add_edge` needs an existing vertex. A call on an empty polygon does nothing. `test_point` and `test_edge` do not modify the polygon. An empty vertex list gives area 0 and perimeter 0.

```mbt check
///|
test "polygon area" {
  let g = @geodesic.Geodesic::wgs84()
  let r = g.polygon_area([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])
  assert_true(r.area > 0.0)
  assert_true(r.perimeter > 0.0)
  let cw = g.polygon_area([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
  assert_true((cw.area + r.area).abs() < 1.0)
  let line = g.polygon(polyline=true)
  line.add_point(0.0, 179.0)
  line.add_point(0.0, -179.0)
  assert_true(line.compute().area.is_nan())
  assert_true(line.compute().perimeter > 0.0)
}
```

## Rhumb lines

恒向线（等方位角）。反算、正算、赤道与航线之间的面积。

```moonbit nocheck
pub fn Rhumb::wgs84() -> Rhumb
pub fn Rhumb::grs80() -> Rhumb
pub fn Rhumb::cgcs2000() -> Rhumb
pub fn Rhumb::new(a : Double, f : Double) -> Rhumb raise RhumbError

pub fn Rhumb::inverse(
  Self, lat1 : Double, lon1 : Double, lat2 : Double, lon2 : Double,
) -> RhumbInverseResult

pub fn Rhumb::direct(
  Self, lat1 : Double, lon1 : Double, azi12 : Double, s12 : Double,
  long_unroll? : Bool,
) -> RhumbDirectResult

pub fn Rhumb::distance(
  Self, lat1 : Double, lon1 : Double, lat2 : Double, lon2 : Double,
) -> Double

pub fn Rhumb::line(Self, lat1 : Double, lon1 : Double, azi12 : Double) -> RhumbLine
pub fn RhumbLine::position(Self, s12 : Double, long_unroll? : Bool) -> RhumbDirectResult
```

`RhumbInverseResult` carries the endpoints, `s12` (m), constant azimuth `azi12` (degrees), and `area12` (m²). `RhumbDirectResult` carries the same quantities for the point reached by travelling `s12` metres at that azimuth.

`area12` is the area between the rhumb line and the equator: the quadrilateral `(lat1, lon1)`, `(0, lon1)`, `(0, lon2)`, `(lat2, lon2)`, counter-clockwise positive.

The shortest rhumb line is returned. It wraps at most halfway around the earth. On opposite meridians the east-going line is chosen.

`distance` returns only `s12`. `line` fixes the azimuth; `position` then takes metres along it. `s12` may be negative.

Longitudes are reduced to [-180, 180]. Set `long_unroll=true` when the path crosses the antimeridian and `lon2` must stay continuous.

Latitude ±90 is a valid pole. A latitude outside [-90, 90] is returned as NaN. No error is raised. If the rhumb line crosses a pole, `lon2` and `area12` are NaN.

The series assume `|f|` up to about 1/150. WGS-84, GRS80, and CGCS2000 are inside that range.

```mbt check
///|
test "rhumb line" {
  let rh = @geodesic.Rhumb::wgs84()
  let r = rh.inverse(51.4778, -0.0015, 35.6762, 139.6503)
  let g = @geodesic.Geodesic::wgs84().distance(
    51.4778, -0.0015, 35.6762, 139.6503,
  )
  assert_true(r.s12 > g)
  assert_true(r.s12 > 1.1e7)
  let d = rh.direct(0.0, 179.0, 90.0, 500000.0, long_unroll=true)
  assert_true(d.lon2 > 180.0)
  let pole = rh.direct(80.0, 10.0, 0.0, 3000000.0)
  assert_true(pole.lon2.is_nan())
  assert_true(pole.area12.is_nan())
}
```

## Intersections

两条测地线的交点。最近点、下一点、线段。

```moonbit nocheck
pub fn Intersect::wgs84() -> Intersect
pub fn Intersect::new(Geodesic) -> Intersect raise IntersectError

pub fn Intersect::closest(
  Self,
  lat_x : Double, lon_x : Double, azi_x : Double,
  lat_y : Double, lon_y : Double, azi_y : Double,
  x0? : Double, y0? : Double,
) -> IntersectPoint raise IntersectError

pub fn Intersect::next(
  Self, lat : Double, lon : Double, azi_x : Double, azi_y : Double,
) -> IntersectPoint raise IntersectError

pub fn Intersect::segment(
  Self,
  lat_x1 : Double, lon_x1 : Double, lat_x2 : Double, lon_x2 : Double,
  lat_y1 : Double, lon_y1 : Double, lat_y2 : Double, lon_y2 : Double,
) -> SegmentIntersect raise IntersectError

pub fn Intersect::dist(Double, Double, x0? : Double, y0? : Double) -> Double
```

`IntersectPoint` is the signed distances `x` and `y` in metres, and `c`. `x` is measured along geodesic X from its start. `y` is measured along geodesic Y. `c` is 0 when the geodesics only cross. `c` is +1 when they lie on top of each other and point the same way. `c` is -1 when they lie on top of each other and point opposite ways.

`closest` takes a point and an azimuth for each geodesic. It minimises `|x - x0| + |y - y0|`. `x0` and `y0` default to 0.

`next` takes one point and two azimuths. The geodesics already meet there. The result is the next closest meeting, not `(0, 0)`.

`segment` takes the endpoints of two shortest segments. `segmode` is 0 when the intersection lies on both segments. Otherwise `segmode` is `3 * kx + ky`. `kx` is -1 when `x < 0`, 0 when `x` is on segment X, and 1 when `x` is past the end of X. `ky` is the same flag for Y.

`dist` is the L1 distance `|x - x0| + |y - y0|`.

`new` raises `IntersectError::TooEccentric` when the ellipsoid is too flat for the search. WGS-84, GRS80, and CGCS2000 are inside the accepted range. The series geodesic is the one used here. It is accurate for `|f|` up to about 1/50.

Two closest intersections can be the same distance apart. `closest`, `next`, and `segment` then return one of them. `segment` needs one shortest geodesic between each pair of endpoints.

```mbt check
///|
test "geodesic intersection" {
  let ix = @geodesic.Intersect::wgs84()
  let p = ix.closest(0.0, 0.0, 45.0, 45.0, 10.0, 135.0)
  let g = @geodesic.Geodesic::wgs84()
  let x = g.direct(0.0, 0.0, 45.0, p.x)
  let y = g.direct(45.0, 10.0, 135.0, p.y)
  assert_true(g.distance(x.lat2, x.lon2, y.lat2, y.lon2) < 1.0e-6)
  let s = ix.segment(0.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, 0.0)
  assert_eq(s.segmode, 0)
  let n = ix.next(0.0, 0.0, 45.0, 135.0)
  assert_true(n.x.abs() + n.y.abs() > 1.0e7)
}
```

## Local cartesian (ENU)

局部东北天直角坐标。相对原点的 `(x, y, z)`，经地心地固转换。

```moonbit nocheck
pub fn Geocentric::wgs84() -> Geocentric
pub fn Geocentric::grs80() -> Geocentric
pub fn Geocentric::cgcs2000() -> Geocentric
pub fn Geocentric::new(a : Double, f : Double) -> Geocentric raise GeocentricError
pub fn Geocentric::from_geodesic(Geodesic) -> Geocentric
pub fn Geocentric::forward(Self, lat : Double, lon : Double, h : Double) -> GeocentricPoint
pub fn Geocentric::reverse(Self, x : Double, y : Double, z : Double) -> GeodeticPosition

pub fn LocalCartesian::wgs84(lat0 : Double, lon0 : Double, h0? : Double) -> LocalCartesian
pub fn LocalCartesian::grs80(lat0 : Double, lon0 : Double, h0? : Double) -> LocalCartesian
pub fn LocalCartesian::new(
  a : Double, f : Double, lat0 : Double, lon0 : Double, h0? : Double,
) -> LocalCartesian raise GeocentricError
pub fn LocalCartesian::from_geodesic(
  Geodesic, lat0 : Double, lon0 : Double, h0? : Double,
) -> LocalCartesian
pub fn LocalCartesian::reset(Self, lat0 : Double, lon0 : Double, h0? : Double) -> Unit
pub fn LocalCartesian::forward(Self, lat : Double, lon : Double, h : Double) -> LocalCartesianPoint
pub fn LocalCartesian::reverse(Self, x : Double, y : Double, z : Double) -> GeodeticPosition
```

`h` and `h0` are ellipsoidal heights in metres. Local `x` is east, `y` is north, and `z` is up (normal to the ellipsoid). The origin is `(lat0, lon0, h0)`; `h0` defaults to 0. `Geocentric` is the ECEF conversion used internally; it is also useful on its own.

`new` raises `GeocentricError::InvalidEllipsoid` when `a` is not positive or `f >= 1`. Latitudes outside [-90°, 90°] yield NaN, as elsewhere in this package.

```mbt check
///|
test "local cartesian ENU" {
  let lc = @geodesic.LocalCartesian::wgs84(39.9087, 116.3975)
  let p = lc.forward(39.9097, 116.3985, 50.0)
  assert_true(p.x > 0.0 && p.y > 0.0 && p.z > 0.0)
  let g = lc.reverse(p.x, p.y, p.z)
  assert_true((g.lat - 39.9097).abs() < 1.0e-10)
  assert_true((g.lon - 116.3985).abs() < 1.0e-10)
  assert_true((g.h - 50.0).abs() < 1.0e-7)
  let ecef = @geodesic.Geocentric::wgs84().forward(39.9087, 116.3975, 0.0)
  assert_true(ecef.x.abs() > 1.0e6 && ecef.z.abs() > 1.0e6)
}
```
