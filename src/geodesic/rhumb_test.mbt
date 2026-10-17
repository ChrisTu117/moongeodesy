///|
fn rhumb_max(errs : Array[Double]) -> Double {
  let mut m = 0.0
  for e in errs {
    if e > m || e.is_nan() {
      m = e
    }
  }
  m
}

///|
fn rhumb_ang_err(a : Double, b : Double) -> Double {
  if a.is_nan() && b.is_nan() {
    return 0.0
  }
  let mut d = (a - b).abs()
  if d > 180.0 {
    d = 360.0 - d
  }
  d
}

///|
fn rhumb_abs_err(a : Double, b : Double) -> Double {
  if (a.is_nan() && b.is_nan()) || a == b {
    0.0
  } else {
    (a - b).abs()
  }
}

///|
test "rhumb inverse matches RhumbSolve" {
  let rh = @geodesic.Rhumb::wgs84()
  let ds = []
  let da = []
  let dS = []
  for row in rhumb_inverse {
    let r = rh.inverse(row[0], row[1], row[2], row[3])
    ds.push(rhumb_abs_err(r.s12, row[5]))
    dS.push(rhumb_abs_err(r.area12, row[6]) / rh.ellipsoid_area())
    if row[0].abs() < 89.9999 && row[2].abs() < 89.9999 && row[5] > 0.0 {
      da.push(
        rhumb_ang_err(r.azi12, row[4]) * row[5] * 3.141592653589793 / 180.0,
      )
    }
  }
  println(
    "rhumb inverse max |ds12| m = \{rhumb_max(ds)}, max azi err (m at distance) = \{rhumb_max(da)}, max |dS12|/area = \{rhumb_max(dS)}",
  )
  assert_true(rhumb_max(ds) < 1.0e-6)
  assert_true(rhumb_max(da) < 1.0e-6)
  assert_true(rhumb_max(dS) < 1.0e-15)
}

///|
fn rhumb_pos_err(
  lat : Double,
  lon : Double,
  rlat : Double,
  rlon : Double,
) -> Double {
  if lat.is_nan() || rlat.is_nan() {
    return if lat.is_nan() && rlat.is_nan() { 0.0 } else { 1.0e9 }
  }
  let dlat = (lat - rlat).abs() * 111319.5
  let dlon = if lon.is_nan() && rlon.is_nan() {
    0.0
  } else if lon.is_nan() || rlon.is_nan() {
    1.0e9
  } else {
    (lon - rlon).abs() * 111319.5 * @math.cos(rlat * 3.141592653589793 / 180.0)
  }
  dlat + dlon
}

///|
test "rhumb direct matches RhumbSolve" {
  let rh = @geodesic.Rhumb::wgs84()
  let ep = []
  let eS = []
  let eu = []
  for row in rhumb_direct {
    let r = rh.direct(row[0], row[1], row[2], row[3])
    ep.push(rhumb_pos_err(r.lat2, r.lon2, row[4], row[5]))
    eS.push(rhumb_abs_err(r.area12, row[6]) / rh.ellipsoid_area())
  }
  for row in rhumb_direct_unroll {
    let r = rh.direct(row[0], row[1], row[2], row[3], long_unroll=true)
    eu.push(rhumb_abs_err(r.lon2, row[5]))
    eS.push(rhumb_abs_err(r.area12, row[6]) / rh.ellipsoid_area())
  }
  println(
    "rhumb direct max pos err m = \{rhumb_max(ep)}, max |dlon| unrolled deg = \{rhumb_max(eu)}, max |dS12|/area = \{rhumb_max(eS)}",
  )
  assert_true(rhumb_max(ep) < 1.0e-6)
  assert_true(rhumb_max(eu) < 1.0e-9)
  assert_true(rhumb_max(eS) < 1.0e-14)
}

///|
test "rhumb direct and inverse agree" {
  let rh = @geodesic.Rhumb::wgs84()
  let d = rh.direct(22.5431, 114.0579, 47.0, 8000000.0)
  let i = rh.inverse(22.5431, 114.0579, d.lat2, d.lon2)
  assert_true((i.s12 - 8000000.0).abs() < 1.0e-6)
  assert_true(rhumb_ang_err(i.azi12, 47.0) < 1.0e-9)
  let back = rh.direct(d.lat2, d.lon2, i.azi12, -i.s12)
  assert_true((back.lat2 - 22.5431).abs() < 1.0e-9)
  assert_true((back.lon2 - 114.0579).abs() < 1.0e-9)
}

///|
test "rhumb edges: poles, antimeridian, east-west, sphere" {
  let rh = @geodesic.Rhumb::wgs84()
  let east = rh.inverse(0.0, 179.0, 0.0, -179.0)
  assert_eq(east.azi12, 90.0)
  assert_true(
    (east.s12 - 6378137.0 * 2.0 * 3.141592653589793 / 180.0).abs() < 1.0e-6,
  )
  assert_eq(east.area12, 0.0)
  let west = rh.inverse(45.0, 20.0, 45.0, 10.0)
  assert_eq(west.azi12, -90.0)
  let pole = rh.inverse(90.0, 0.0, 0.0, 25.0)
  assert_eq(pole.azi12, 180.0)
  assert_true((pole.s12 - 10001965.729).abs() < 1.0e-2)
  let same = rh.inverse(12.5, 400.0, 12.5, -320.0)
  assert_eq(same.s12, 0.0)
  assert_eq(same.lon1, 40.0)
  let cross = rh.direct(80.0, 10.0, 0.0, 3000000.0)
  assert_true(cross.lon2.is_nan())
  assert_true(cross.area12.is_nan())
  assert_true(cross.lat2 < 80.0)
  let un = rh.direct(0.0, 179.0, 90.0, 500000.0, long_unroll=true)
  assert_true(un.lon2 > 180.0)
  let bad = rh.inverse(91.0, 0.0, 0.0, 0.0)
  assert_true(bad.lat1.is_nan())
  assert_true(bad.s12.is_nan())
  let sphere = @geodesic.Rhumb::new(6371000.0, 0.0)
  let q = sphere.inverse(0.0, 0.0, 0.0, 90.0)
  assert_true((q.s12 - 6371000.0 * 3.141592653589793 / 2.0).abs() < 1.0e-6)
  let mer = sphere.inverse(0.0, 0.0, 90.0, 0.0)
  assert_true((mer.s12 - q.s12).abs() < 1.0e-6)
  assert_eq(mer.azi12, 0.0)
  let line = rh.line(10.0, 20.0, 450.0)
  assert_eq(line.azimuth(), 90.0)
  assert_eq(line.latitude(), 10.0)
  assert_eq(line.longitude(), 20.0)
  let cgcs = @geodesic.Rhumb::cgcs2000()
  let dw = rh.distance(39.9087, 116.3975, 31.2304, 121.4737)
  let dc = cgcs.distance(39.9087, 116.3975, 31.2304, 121.4737)
  assert_true((dw - dc).abs() > 0.0)
  assert_true((dw - dc).abs() < 1.0)
  let invalid = try {
    ignore(@geodesic.Rhumb::new(0.0, 0.0))
    false
  } catch {
    _ => true
  }
  assert_true(invalid)
}
