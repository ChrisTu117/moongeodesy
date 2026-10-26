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

///|
/// Just inside / just outside the usual mainland bounding box
/// (lon 72.004°..137.8347°, lat 0.8293°..55.8271°). Outside → GCJ-02 is
/// identity; inside → offset is applied. The box is not the national border.
test "out_of_china box edges" {
  // West edge at lon = 72.004
  assert_true(@cn.out_of_china(35.0, 72.003999))
  assert_false(@cn.out_of_china(35.0, 72.004))
  let west_out = @cn.LatLon::new(35.0, 72.003999)
  let west_in = @cn.LatLon::new(35.0, 72.004)
  assert_eq(@cn.wgs84_to_gcj02(west_out), west_out)
  assert_true((@cn.wgs84_to_gcj02(west_in).lon - west_in.lon).abs() > 1.0e-4)

  // East edge at lon = 137.8347
  assert_true(@cn.out_of_china(35.0, 137.8347001))
  assert_false(@cn.out_of_china(35.0, 137.8347))
  let east_out = @cn.LatLon::new(35.0, 137.8347001)
  let east_in = @cn.LatLon::new(35.0, 137.8347)
  assert_eq(@cn.wgs84_to_gcj02(east_out), east_out)
  assert_true((@cn.wgs84_to_gcj02(east_in).lon - east_in.lon).abs() > 1.0e-4)

  // South edge at lat = 0.8293
  assert_true(@cn.out_of_china(0.8292999, 105.0))
  assert_false(@cn.out_of_china(0.8293, 105.0))
  let south_out = @cn.LatLon::new(0.8292999, 105.0)
  let south_in = @cn.LatLon::new(0.8293, 105.0)
  assert_eq(@cn.wgs84_to_gcj02(south_out), south_out)
  assert_true((@cn.wgs84_to_gcj02(south_in).lat - south_in.lat).abs() > 1.0e-4)

  // North edge at lat = 55.8271
  assert_true(@cn.out_of_china(55.8271001, 105.0))
  assert_false(@cn.out_of_china(55.8271, 105.0))
  let north_out = @cn.LatLon::new(55.8271001, 105.0)
  let north_in = @cn.LatLon::new(55.8271, 105.0)
  assert_eq(@cn.wgs84_to_gcj02(north_out), north_out)
  assert_true((@cn.wgs84_to_gcj02(north_in).lat - north_in.lat).abs() > 1.0e-4)
}

///|
/// Hong Kong, Macau, and Taipei sit inside the usual box, so GCJ-02 is
/// applied even though they are often treated specially in map products.
/// Exact inverses must still round-trip.
test "HK Macau Taipei inside box and round-trip" {
  let cities : Array[(String, Double, Double)] = [
    ("Hong Kong", 22.3193, 114.1694),
    ("Macau", 22.1987, 113.5439),
    ("Taipei", 25.0330, 121.5654),
  ]
  for city in cities {
    let (name, lat, lon) = city
    assert_false(@cn.out_of_china(lat, lon), msg=name)
    let w = @cn.LatLon::new(lat, lon)
    let g = @cn.wgs84_to_gcj02(w)
    assert_true(
      (g.lat - w.lat).abs() > 1.0e-4 || (g.lon - w.lon).abs() > 1.0e-4,
      msg="\{name} GCJ offset",
    )
    let back = @cn.gcj02_to_wgs84(g)
    assert_true(close(back.lat, w.lat, 1.0e-9), msg="\{name} lat")
    assert_true(close(back.lon, w.lon, 1.0e-9), msg="\{name} lon")
    let bd = @cn.wgs84_to_bd09(w)
    let w2 = @cn.bd09_to_wgs84(bd)
    assert_true(close(w2.lat, w.lat, 1.0e-9), msg="\{name} bd lat")
    assert_true(close(w2.lon, w.lon, 1.0e-9), msg="\{name} bd lon")
  }
}

///|
/// A point south of lat 0.8293° (still in SE Asia longitudes): GCJ-02 is
/// identity; BD-09 is still applied everywhere.
test "south of box: GCJ identity, BD still offsets" {
  let p = @cn.LatLon::new(0.5, 105.0)
  assert_true(@cn.out_of_china(p.lat, p.lon))
  assert_eq(@cn.wgs84_to_gcj02(p), p)
  assert_eq(@cn.gcj02_to_wgs84(p), p)
  let bd = @cn.wgs84_to_bd09(p)
  assert_true(
    (bd.lat - p.lat).abs() > 1.0e-6 || (bd.lon - p.lon).abs() > 1.0e-6,
  )
  let back = @cn.bd09_to_wgs84(bd)
  assert_true(close(back.lat, p.lat, 1.0e-9))
  assert_true(close(back.lon, p.lon, 1.0e-9))
}
