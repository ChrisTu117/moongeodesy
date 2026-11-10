// Behaviours promised in README.mbt.md that had no test or a looser one.

///|
test "LatLon shows as LatLon(lat, lon)" {
  inspect(@cn.LatLon::new(1.5, 2.0), content="LatLon(1.5, 2)")
}

///|
test "approximate inverses: measured maxima match the README" {
  let g = @geodesic.Geodesic::wgs84()
  let mut gcj_m = 0.0
  let mut bd_deg = 0.0
  for v in cn_vectors {
    let (lat, lon, glat, glon, blat, blon) = v
    if !@cn.out_of_china(lat, lon) {
      let w = @cn.gcj02_to_wgs84_approx(@cn.LatLon::new(glat, glon))
      let s = g.inverse(lat, lon, w.lat, w.lon).s12
      if s > gcj_m {
        gcj_m = s
      }
    }
    let b = @cn.bd09_to_gcj02_approx(@cn.LatLon::new(blat, blon))
    let d = (b.lat - glat).abs()
    let e = (b.lon - glon).abs()
    if d > bd_deg {
      bd_deg = d
    }
    if e > bd_deg {
      bd_deg = e
    }
  }
  println(
    "gcj02_to_wgs84_approx max error \{gcj_m} m; bd09_to_gcj02_approx max error \{bd_deg} deg",
  )
  assert_true(gcj_m < 3.5) // README: up to about 3 m
  assert_true(bd_deg < 1.5e-6) // README: up to about 1.5e-6 degrees
}
