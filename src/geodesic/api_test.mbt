///|
test "known distances" {
  let g = @geodesic.Geodesic::wgs84()
  // Beijing Tiananmen → Shanghai People's Square
  let d = g.distance(39.9087, 116.3975, 31.2304, 121.4737)
  assert_true((d - 1066676.4).abs() < 1.0)
  // quarter meridian (equator → pole) on WGS-84: 10001965.729 m
  let q = g.inverse(0.0, 0.0, 90.0, 0.0)
  assert_true((q.s12 - 10001965.729).abs() < 1.0e-3)
  assert_eq(q.azi1, 0.0)
  // half the equator
  let e = g.inverse(0.0, 0.0, 0.0, 179.0)
  assert_true(
    (e.s12 - 6378137.0 * 179.0 * 3.141592653589793 / 180.0).abs() < 1.0e-6,
  )
  // coincident points
  let z = g.inverse(12.5, 34.5, 12.5, 34.5)
  assert_eq(z.s12, 0.0)
}

///|
test "direct and inverse are consistent" {
  let g = @geodesic.Geodesic::wgs84()
  let d = g.direct(22.5431, 114.0579, 47.0, 8000000.0)
  let i = g.inverse(22.5431, 114.0579, d.lat2, d.lon2)
  assert_true((i.s12 - 8000000.0).abs() < 1.0e-8)
  assert_true((i.azi1 - 47.0).abs() < 1.0e-12)
  assert_true((i.azi2 - d.azi2).abs() < 1.0e-12)
}

///|
test "geodesic line and densify" {
  let g = @geodesic.Geodesic::wgs84()
  let l = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
  let pts = l.densify(100000.0)
  // 9899.5 km in steps of at most 100 km → 99 segments, 100 points
  assert_eq(pts.length(), 100)
  assert_true((pts[0].0 - 31.2304).abs() < 1.0e-12)
  assert_true((pts[99].0 - 37.7749).abs() < 1.0e-9)
  assert_true((pts[99].1 - -122.4194).abs() < 1.0e-9)
  let mid = l.position(l.distance() / 2.0)
  let a = g.distance(31.2304, 121.4737, mid.lat2, mid.lon2)
  assert_true((a - l.distance() / 2.0).abs() < 1.0e-6)
  let dl = g.direct_line(0.0, 0.0, 90.0, 1000.0)
  assert_true(
    (dl.arc() -
    1000.0 / 6378137.0 * 180.0 / 3.141592653589793 / (1.0 - 1.0 / 298.257223563)).abs() <
    1.0e-9,
  )
}

///|
test "custom ellipsoids" {
  // a sphere: geodesics are great circles
  let s = @geodesic.Geodesic::new(6371000.0, 0.0)
  let r = s.inverse(0.0, 0.0, 0.0, 90.0)
  assert_true((r.s12 - 6371000.0 * 3.141592653589793 / 2.0).abs() < 1.0e-6)
  // invalid parameters raise
  let bad = try {
    ignore(@geodesic.Geodesic::new(-1.0, 0.0))
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}

///|
test "predefined ellipsoids" {
  let w = @geodesic.Geodesic::wgs84()
  let c = @geodesic.Geodesic::cgcs2000()
  assert_eq(c.equatorial_radius(), 6378137.0)
  // WGS-84 and CGCS2000 differ by ~0.1 mm in the polar radius
  let dw = w.distance(39.9087, 116.3975, 31.2304, 121.4737)
  let dc = c.distance(39.9087, 116.3975, 31.2304, 121.4737)
  assert_true((dw - dc).abs() < 1.0e-3)
  assert_true(@geodesic.Geodesic::grs80().flattening() == c.flattening())
}
