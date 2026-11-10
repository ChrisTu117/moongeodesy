// Asserts the values stated in the comments of README.md's Usage tour.

///|
test "README usage comments" {
  let g = @geodesic.Geodesic::wgs84()
  let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737)
  // r.s12 = 1066676.4 m, r.azi1 = 152.94°
  assert_true((r.s12 - 1066676.4).abs() < 0.05)
  assert_true((r.azi1 - 152.94).abs() < 0.005)
  // waypoints at most 100 km apart
  let line = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
  let pts = line.densify(100000.0)
  for i in 0..<(pts.length() - 1) {
    assert_true(
      g.inverse(pts[i].0, pts[i].1, pts[i + 1].0, pts[i + 1].1).s12 <=
      100000.0 + 1.0e-6,
    )
  }
  // counter-clockwise positive
  let a = g.polygon_area([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])
  assert_true(a.area > 0.0) // east, north, west: counter-clockwise
  let cw = g.polygon_area([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
  assert_true(cw.area < 0.0)
  // London → Tokyo rhumb line is longer than the geodesic
  let h = @geodesic.Rhumb::wgs84().inverse(51.4778, -0.0015, 35.6762, 139.6503)
  assert_true(h.s12 > g.inverse(51.4778, -0.0015, 35.6762, 139.6503).s12)
  // round trip < 1e-9°
  let w = @cn.LatLon::new(39.9087, 116.3975)
  let back = @cn.bd09_to_wgs84(@cn.gcj02_to_bd09(@cn.wgs84_to_gcj02(w)))
  assert_true(
    (back.lat - w.lat).abs() < 1.0e-9 && (back.lon - w.lon).abs() < 1.0e-9,
  )
  // zone 50N, "50SMK4850217797", centre of the 1 m square
  let u = @utm.from_latlon(39.9087, 116.3975)
  assert_true(u.zone == 50 && u.north)
  assert_eq(@utm.to_mgrs(39.9087, 116.3975), "50SMK4850217797")
  let p = @utm.from_mgrs("50S MK 48502 17797")
  assert_eq(p.easting, 448502.5)
  // gp.easting = 956117.8 m
  let albers = @conic.ConicGrid::albers(
    6378137.0,
    1.0 / 298.257222101,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  )
  assert_true(
    (albers.forward(39.9087, 116.3975).easting - 956117.8).abs() < 0.05,
  )
  // mag.decl = -7.56° (west)
  let t = @magnetic.decimal_year(2026, 11, 10)
  let mag = @magnetic.MagneticModel::wmm2025().field(39.9087, 116.3975, 0.0, t)
  assert_true((mag.decl - -7.56).abs() < 0.005)
  // 中文简介: 北京磁偏角 ≈ -7.56° at 2026.86, r.s12 ≈ 1066.68 km
  let zh = @magnetic.MagneticModel::wmm2025().field(
    39.9087, 116.3975, 0.0, 2026.86,
  )
  assert_true((zh.decl - -7.56).abs() < 0.005)
  assert_true((r.s12 / 1000.0 - 1066.68).abs() < 0.005)
}
