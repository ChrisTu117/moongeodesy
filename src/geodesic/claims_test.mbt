// Small behaviours promised in README.mbt.md that had no test.

///|
test "Geodesic::line has no distance until set_distance or set_arc" {
  let g = @geodesic.Geodesic::wgs84()
  let line = g.line(30.0, 120.0, 45.0)
  assert_true(line.distance().is_nan())
  line.set_distance(1.0e6)
  assert_eq(line.distance(), 1.0e6)
  let d = g.direct(30.0, 120.0, 45.0, 1.0e6)
  assert_true((line.arc() - d.a12).abs() < 1.0e-12)
  line.set_arc(10.0)
  let p = line.arc_position(10.0)
  assert_true((line.distance() - p.s12).abs() < 1.0e-6)
  let q = g.direct(30.0, 120.0, 45.0, p.s12)
  assert_true(
    (q.lat2 - p.lat2).abs() < 1.0e-12 && (q.lon2 - p.lon2).abs() < 1.0e-12,
  )
}

///|
test "densify keeps consecutive points at most max_step apart" {
  let g = @geodesic.Geodesic::wgs84()
  let line = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
  let pts = line.densify(100000.0)
  for i in 0..<(pts.length() - 1) {
    let s = g.inverse(pts[i].0, pts[i].1, pts[i + 1].0, pts[i + 1].1).s12
    assert_true(s <= 100000.0 + 1.0e-6)
  }
}

///|
test "empty polygons and add_edge on an empty polygon" {
  let g = @geodesic.Geodesic::wgs84()
  let p = g.polygon()
  p.add_edge(0.0, 1000.0)
  assert_eq(p.compute().num, 0)
  let empty : Array[(Double, Double)] = []
  let r = g.polygon_area(empty[:])
  assert_eq(r.area, 0.0)
  assert_eq(r.perimeter, 0.0)
}

///|
test "bad ellipsoids raise; latitudes beyond 90 give NaN" {
  let raised = (a : Double, f : Double) => {
    try {
      ignore(@geodesic.Geodesic::new(a, f))
      false
    } catch {
      _ => true
    }
  }
  assert_true(raised(6378137.0, 1.0))
  assert_true(raised(@double.not_a_number, 0.0))
  assert_true(raised(@double.infinity, 0.0))
  let geoc = (a : Double, f : Double) => {
    try {
      ignore(@geodesic.Geocentric::new(a, f))
      false
    } catch {
      _ => true
    }
  }
  assert_true(geoc(6378137.0, 1.0))
  assert_true(geoc(-1.0, 0.0))
  let g = @geodesic.Geodesic::wgs84()
  let d = g.direct(91.0, 0.0, 10.0, 1000.0)
  assert_true(d.lat2.is_nan())
  assert_true(g.inverse(0.0, 0.0, 91.0, 0.0).s12.is_nan())
  let p = @geodesic.Geocentric::wgs84().forward(91.0, 0.0, 0.0)
  assert_true(p.x.is_nan() || p.z.is_nan())
}

///|
test "very short lines (accurate hypot) match geographiclib" {
  // geographiclib 2.1 (Python): Geodesic.WGS84.Inverse for points 1e-7 and
  // 1e-9 degrees apart. These need a correctly rounded hypot.
  let g = @geodesic.Geodesic::wgs84()
  let r1 = g.inverse(40.0, 116.0, 40.0000001, 116.0000001)
  assert_true((r1.s12 - 0.014007426953374454).abs() < 1.0e-15)
  let r2 = g.inverse(-33.0, 151.0, -33.000000001, 151.000000002)
  assert_true((r2.s12 - 2.1733383453160246e-4).abs() < 1.0e-17)
}
