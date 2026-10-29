///|
test "GeodesicExact inverse and direct at f = 0.1" {
  let g = @geodesic.GeodesicExact::new(6378137.0, 0.1)
  let r = g.inverse(40.6, -73.8, 1.5, 0.0)
  assert_true((r.s12 - 8323146.9846059354).abs() < 1.0e-6)
  assert_true((r.azi1 - 94.973407482443889).abs() < 1.0e-11)
  assert_true((r.azi2 - 127.90469767058218).abs() < 1.0e-11)
  assert_true((r.a12 - 81.390825504822601).abs() < 1.0e-11)
  assert_true((r.m12 - 5882704.3517672326).abs() < 1.0e-5)
  assert_true((r.scale12 - 0.20329884227443962).abs() < 1.0e-12)
  assert_true((r.scale21 - 0.15107485075491633).abs() < 1.0e-12)
  assert_true(
    (r.area12 - 19678368251345.023).abs() / g.ellipsoid_area() < 1.0e-12,
  )
  let d = g.direct(40.6, -73.8, 45.0, 10000000.0)
  assert_true((d.lat2 - 33.830172877019528).abs() < 1.0e-11)
  assert_true((d.lon2 - 45.621635867229685).abs() < 1.0e-11)
  assert_true((d.azi2 - 139.16745368446479).abs() < 1.0e-11)
  assert_true((d.a12 - 93.948275977567477).abs() < 1.0e-11)
  assert_true((d.m12 - 6544871.7763396213).abs() < 1.0e-4)
  assert_true((d.scale12 - 0.035541045029230806).abs() < 1.0e-12)
  assert_true((d.scale21 - 0.01760011296392848).abs() < 1.0e-12)
  assert_true(
    (d.area12 - 59075487728653.039).abs() / g.ellipsoid_area() < 1.0e-12,
  )
  let back = g.inverse(40.6, -73.8, d.lat2, d.lon2)
  assert_true((back.s12 - 10000000.0).abs() < 1.0e-6)
}

///|
test "GeodesicExact rejects a non-positive axis" {
  let bad = try {
    ignore(@geodesic.GeodesicExact::new(0.0, 0.1))
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let flat = try {
    ignore(@geodesic.GeodesicExact::new(6378137.0, 1.0))
    false
  } catch {
    _ => true
  }
  assert_true(flat)
}

///|
fn exact_geod(a : Double, f : Double) -> @geodesic.GeodesicExact {
  @geodesic.GeodesicExact::new(a, f) catch {
    _ => abort("bad ellipsoid")
  }
}

///|
fn dmax(a : Double, b : Double) -> Double {
  if b > a {
    b
  } else {
    a
  }
}

///|
test "GeodesicExact matches GeodSolve -E for |f| > 0.02" {
  let mut ds = 0.0
  let mut da = 0.0
  let mut dm = 0.0
  let mut dS = 0.0
  let mut d_scale = 0.0
  let mut dp = 0.0
  let mut n_inv = 0
  let mut n_dir = 0
  for row in exact_cases {
    let g = exact_geod(row[0], row[1])
    let area0 = g.ellipsoid_area()
    if row[2] == 0.0 {
      n_inv = n_inv + 1
      let r = g.inverse(row[3], row[4], row[5], row[6])
      ds = dmax(ds, (r.s12 - row[13]).abs())
      dm = dmax(dm, (r.m12 - row[15]).abs())
      d_scale = dmax(d_scale, (r.scale12 - row[16]).abs())
      d_scale = dmax(d_scale, (r.scale21 - row[17]).abs())
      dS = dmax(dS, (r.area12 - row[18]).abs() / area0)
      if row[14] <= 179.5 && row[3].abs() < 89.999 && row[5].abs() < 89.999 {
        da = dmax(da, ang_err(r.azi1, row[9]))
        da = dmax(da, ang_err(r.azi2, row[12]))
      }
    } else {
      n_dir = n_dir + 1
      let r = g.direct(row[3], row[4], row[5], row[6])
      dp = dmax(dp, (r.lat2 - row[10]).abs())
      dp = dmax(dp, (r.lon2 - row[11]).abs())
      da = dmax(da, ang_err(r.azi2, row[12]))
      dm = dmax(dm, (r.m12 - row[15]).abs())
      d_scale = dmax(d_scale, (r.scale12 - row[16]).abs())
      d_scale = dmax(d_scale, (r.scale21 - row[17]).abs())
      dS = dmax(dS, (r.area12 - row[18]).abs() / area0)
      ds = dmax(ds, (r.s12 - row[13]).abs())
    }
  }
  println(
    "exact n_inv=\{n_inv} n_dir=\{n_dir} max |ds| m=\{ds} max |dazi| deg=\{da} max |dm| m=\{dm} max |dM|=\{d_scale} max |dS|/area=\{dS} max |dpos| deg=\{dp}",
  )
  assert_true(n_inv == 105)
  assert_true(n_dir == 6)
  assert_true(ds < 2.0e-8) // 20 nm
  assert_true(da < 1.0e-12)
  assert_true(dm < 2.0e-8)
  assert_true(d_scale < 1.0e-13)
  assert_true(dS < 1.0e-15)
  assert_true(dp < 1.0e-12)
}

///|
test "GeodesicExact line, arc, and sphere" {
  let g = @geodesic.GeodesicExact::new(6378137.0, 0.1)
  let line = g.inverse_line(10.0, 20.0, -30.0, 140.0)
  let end = line.position(line.distance())
  assert_true((end.lat2 - -30.0).abs() < 1.0e-9)
  assert_true((end.lon2 - 140.0).abs() < 1.0e-9)
  let pts = line.densify(2.0e6)
  assert_true(pts.length() >= 2)
  assert_true((pts[0].0 - 10.0).abs() < 1.0e-12)
  let arc = g.arc_direct(0.0, 0.0, 90.0, 10.0)
  let back = g.inverse(0.0, 0.0, arc.lat2, arc.lon2)
  assert_true((back.a12 - 10.0).abs() < 1.0e-9)
  let sphere = @geodesic.GeodesicExact::new(6371000.0, 0.0)
  let q = sphere.inverse(0.0, 0.0, 90.0, 0.0)
  let quarter = 6371000.0 * 3.141592653589793 / 2.0
  assert_true((q.s12 - quarter).abs() < 1.0e-6)
  assert_true((q.azi1 - 0.0).abs() < 1.0e-12)
}

///|
test "series geodesic is less accurate when |f| = 0.2" {
  let exact = @geodesic.GeodesicExact::new(6378137.0, 0.2)
  let series = @geodesic.Geodesic::new(6378137.0, 0.2)
  let gap = (exact.distance(10.0, 20.0, -30.0, 140.0) -
  series.distance(10.0, 20.0, -30.0, 140.0)).abs()
  println("series gap at f=0.2, 10 20 to -30 140: \{gap} m")
  // The series distance is centimetres off. Exact matches GeodSolve -E.
  assert_true(gap > 0.01)
}
