// White-box test: round trips right next to every land border.

///|
/// Points 0.002° and 0.008° (about 200 m and 900 m) from every other border
/// vertex, in eight directions. Inside China the round trip must be exact;
/// everywhere the inverse must return a WGS-84 source of the GCJ-02 value.
test "round trip next to land borders" {
  let dirs : Array[(Double, Double)] = [
    (1.0, 0.0),
    (-1.0, 0.0),
    (0.0, 1.0),
    (0.0, -1.0),
    (0.7071, 0.7071),
    (0.7071, -0.7071),
    (-0.7071, 0.7071),
    (-0.7071, -0.7071),
  ]
  let n = foreign_xy.length() / 2
  let mut checked_inside = 0
  for v = 0; v < n; v = v + 2 {
    let lon0 = foreign_xy[2 * v]
    let lat0 = foreign_xy[2 * v + 1]
    for d in dirs {
      for r in [0.002, 0.008] {
        let lat = lat0 + r * d.1
        let lon = lon0 + r * d.0
        let w = LatLon::new(lat, lon)
        let g = wgs84_to_gcj02(w)
        let back = gcj02_to_wgs84(g)
        let again = wgs84_to_gcj02(back)
        assert_true(
          (again.lat - g.lat).abs() <= 1.0e-9 &&
          (again.lon - g.lon).abs() <= 1.0e-9,
          msg="source of g at \{lat},\{lon}",
        )
        if !out_of_china(lat, lon) {
          checked_inside += 1
          assert_true(
            (back.lat - lat).abs() <= 1.0e-9 && (back.lon - lon).abs() <= 1.0e-9,
            msg="round trip at \{lat},\{lon}",
          )
        }
      }
    }
  }
  assert_true(checked_inside > 1000, msg="inside points \{checked_inside}")
}
