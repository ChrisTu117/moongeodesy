///|
fn close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
test "wgs84_to_gcj02 matches coord-convert reference" {
  for v in cn_vectors {
    let (lat, lon, glat, glon, _, _) = v
    let g = @cn.wgs84_to_gcj02(@cn.LatLon::new(lat, lon))
    assert_true(close(g.lat, glat, 1.0e-11), msg="lat \{lat},\{lon}")
    assert_true(close(g.lon, glon, 1.0e-11), msg="lon \{lat},\{lon}")
  }
}

///|
test "gcj02_to_bd09 matches coord-convert reference" {
  for v in cn_vectors {
    let (_, _, glat, glon, blat, blon) = v
    let b = @cn.gcj02_to_bd09(@cn.LatLon::new(glat, glon))
    assert_true(close(b.lat, blat, 1.0e-11))
    assert_true(close(b.lon, blon, 1.0e-11))
  }
}

///|
test "exact inverses round-trip below 1e-9 degrees" {
  for v in cn_vectors {
    let (lat, lon, glat, glon, blat, blon) = v
    let w = @cn.gcj02_to_wgs84(@cn.LatLon::new(glat, glon))
    assert_true(close(w.lat, lat, 1.0e-9) && close(w.lon, lon, 1.0e-9))
    let g = @cn.bd09_to_gcj02(@cn.LatLon::new(blat, blon))
    assert_true(close(g.lat, glat, 1.0e-9) && close(g.lon, glon, 1.0e-9))
    let w2 = @cn.bd09_to_wgs84(@cn.LatLon::new(blat, blon))
    assert_true(close(w2.lat, lat, 1.0e-9) && close(w2.lon, lon, 1.0e-9))
  }
}

///|
test "approximate inverses are within a few metres" {
  for v in cn_vectors {
    let (lat, lon, glat, glon, blat, blon) = v
    let w = @cn.gcj02_to_wgs84_approx(@cn.LatLon::new(glat, glon))
    assert_true(close(w.lat, lat, 5.0e-5) && close(w.lon, lon, 5.0e-5))
    let g = @cn.bd09_to_gcj02_approx(@cn.LatLon::new(blat, blon))
    assert_true(close(g.lat, glat, 2.0e-5) && close(g.lon, glon, 2.0e-5))
  }
}

///|
test "identity outside China" {
  let p = @cn.LatLon::new(48.8584, 2.2945)
  assert_true(@cn.out_of_china(p.lat, p.lon))
  assert_eq(@cn.wgs84_to_gcj02(p), p)
  assert_eq(@cn.gcj02_to_wgs84(p), p)
  assert_false(@cn.out_of_china(39.9087, 116.3975))
}

///|
test "Tiananmen snapshot" {
  let g = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975))
  let b = @cn.gcj02_to_bd09(g)
  assert_true(close(g.lat, 39.91010, 1.0e-4) && close(g.lon, 116.40374, 1.0e-4))
  assert_true(close(b.lat, 39.91645, 1.0e-4) && close(b.lon, 116.41012, 1.0e-4))
}
