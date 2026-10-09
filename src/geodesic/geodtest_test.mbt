///|
fn max_err(errs : Array[Double]) -> Double {
  let mut m = 0.0
  for e in errs {
    if e > m || e.is_nan() {
      m = e
    }
  }
  m
}

///|
fn ang_err(a : Double, b : Double) -> Double {
  // angular difference in degrees, wrapped
  let mut d = (a - b).abs()
  if d > 180.0 {
    d = 360.0 - d
  }
  d
}

///|
test "GeodTest inverse: distance, azimuths, reduced length, area" {
  let g = @geodesic.Geodesic::wgs84()
  let es = []
  let ea = []
  let em = []
  let eS = []
  for row in geodtest {
    let r = g.inverse(row[0], row[1], row[3], row[4])
    es.push((r.s12 - row[6]).abs())
    // Nearly antipodal points (a12 ~ 180) can have several shortest
    // geodesics, so azimuths, m12 and S12 are not unique there.
    if row[7] > 179.5 {
      continue
    }
    em.push((r.m12 - row[8]).abs())
    eS.push((r.area12 - row[9]).abs() / g.ellipsoid_area())
    // azimuth error expressed as a displacement (m) at the far end
    if row[0].abs() < 89.9999 && row[3].abs() < 89.9999 {
      ea.push(ang_err(r.azi1, row[2]) * row[6] * 3.14159 / 180.0)
      ea.push(ang_err(r.azi2, row[5]) * row[6] * 3.14159 / 180.0)
    }
  }
  println(
    "inverse max |ds12| m = \{max_err(es)}, max azi err (m at distance) = \{max_err(ea)}, max |dm12| = \{max_err(em)}, max |dS12|/area0 = \{max_err(eS)}",
  )
  assert_true(max_err(es) < 2.0e-8) // 20 nm
  assert_true(max_err(ea) < 1.0e-6) // 1 micron; worst cases are near-antipodal
  assert_true(max_err(em) < 2.0e-8)
  assert_true(max_err(eS) < 1.0e-12)
}

///|
test "GeodTest direct: position, azimuth, area" {
  let g = @geodesic.Geodesic::wgs84()
  let ep = []
  let ea = []
  let eS = []
  for row in geodtest {
    let r = g.direct(row[0], row[1], row[2], row[6], long_unroll=true)
    let dlat = (r.lat2 - row[3]).abs() * 111320.0
    let dlon = (r.lon2 - row[4]).abs() *
      111320.0 *
      @math.cos(row[3] * 3.141592653589793 / 180.0)
    ep.push(dlat + dlon)
    ea.push(ang_err(r.azi2, row[5]))
    eS.push((r.area12 - row[9]).abs() / g.ellipsoid_area())
  }
  println(
    "direct max pos err m = \{max_err(ep)}, max azi2 err deg = \{max_err(ea)}, max |dS12|/area0 = \{max_err(eS)}",
  )
  assert_true(max_err(ep) < 2.0e-8)
  assert_true(max_err(ea) < 1.0e-8)
  assert_true(max_err(eS) < 1.0e-11)
}
