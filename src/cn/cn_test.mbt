///|
fn close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
/// coord-convert only uses the bounding box, so it also offsets points on the
/// land of neighbouring countries. Those points must now be left unchanged;
/// every other point must still match the reference.
test "wgs84_to_gcj02 matches coord-convert reference" {
  let mut compared = 0
  for v in cn_vectors {
    let (lat, lon, glat, glon, _, _) = v
    let w = @cn.LatLon::new(lat, lon)
    let g = @cn.wgs84_to_gcj02(w)
    if @cn.out_of_china(lat, lon) {
      assert_eq(g, w)
      continue
    }
    compared += 1
    assert_true(close(g.lat, glat, 1.0e-11), msg="lat \{lat},\{lon}")
    assert_true(close(g.lon, glon, 1.0e-11), msg="lon \{lat},\{lon}")
  }
  assert_true(compared >= 120, msg="compared \{compared}")
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
    let g = @cn.bd09_to_gcj02(@cn.LatLon::new(blat, blon))
    assert_true(close(g.lat, glat, 1.0e-9) && close(g.lon, glon, 1.0e-9))
    // The reference also offsets foreign points; only Chinese ones apply.
    if !@cn.out_of_china(lat, lon) {
      let w = @cn.gcj02_to_wgs84(@cn.LatLon::new(glat, glon))
      assert_true(close(w.lat, lat, 1.0e-9) && close(w.lon, lon, 1.0e-9))
      let w2 = @cn.bd09_to_wgs84(@cn.LatLon::new(blat, blon))
      assert_true(close(w2.lat, lat, 1.0e-9) && close(w2.lon, lon, 1.0e-9))
    }
  }
}

///|
test "approximate inverses are within a few metres" {
  for v in cn_vectors {
    let (lat, lon, glat, glon, blat, blon) = v
    if !@cn.out_of_china(lat, lon) {
      let w = @cn.gcj02_to_wgs84_approx(@cn.LatLon::new(glat, glon))
      assert_true(close(w.lat, lat, 5.0e-5) && close(w.lon, lon, 5.0e-5))
    }
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
/// Just inside / just outside the usual bounding box
/// (lon 72.004°..137.8347°, lat 0.8293°..55.8271°). `out_of_china_bbox` is the
/// box alone. The "inside" points are over the sea, where GCJ-02 is applied;
/// on the land of other countries it is not (see boundary_test.mbt).
test "out_of_china box edges" {
  let edges : Array[(Double, Double, Double, Double)] = [
    // inside point, then the matching point just outside
    (19.0, 72.004, 19.0, 72.003999), // west edge, Arabian Sea
    (30.0, 137.8347, 30.0, 137.8347001), // east edge, Pacific
    (0.8293, 107.0, 0.8292999, 107.0), // south edge, South China Sea
    (55.8271, 137.8, 55.8271001, 137.8), // north edge, Sea of Okhotsk
  ]
  for e in edges {
    let (ilat, ilon, olat, olon) = e
    assert_false(@cn.out_of_china_bbox(ilat, ilon), msg="\{ilat},\{ilon}")
    assert_true(@cn.out_of_china_bbox(olat, olon), msg="\{olat},\{olon}")
    assert_false(@cn.out_of_china(ilat, ilon), msg="\{ilat},\{ilon}")
    assert_true(@cn.out_of_china(olat, olon), msg="\{olat},\{olon}")
    let inside = @cn.LatLon::new(ilat, ilon)
    let outside = @cn.LatLon::new(olat, olon)
    let g = @cn.wgs84_to_gcj02(inside)
    assert_true(
      (g.lat - ilat).abs() > 1.0e-4 || (g.lon - ilon).abs() > 1.0e-4,
      msg="offset at \{ilat},\{ilon}",
    )
    assert_eq(@cn.wgs84_to_gcj02(outside), outside)
    // The exact inverse returns the inside point even when its GCJ-02 value
    // has crossed the edge (acceptance feedback: ~0.0055° error before).
    let back = @cn.gcj02_to_wgs84(g)
    assert_true(close(back.lat, ilat, 1.0e-9), msg="lat \{ilat},\{ilon}")
    assert_true(close(back.lon, ilon, 1.0e-9), msg="lon \{ilat},\{ilon}")
  }
}

///|
/// Hong Kong, Macau, and Taipei are inside China in the boundary data, so
/// GCJ-02 is applied even though they are often treated specially in map
/// products.
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
