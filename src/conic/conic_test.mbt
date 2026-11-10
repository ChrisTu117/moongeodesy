///|
let deg : Double = 3.141592653589793 / 180.0

///|
fn grid_of(
  s : FixedArray[Double],
  fe~ : Double = 0.0,
  fn_~ : Double = 0.0,
) -> @conic.ConicGrid {
  try {
    match s[0].to_int() {
      0 =>
        @conic.ConicGrid::lambert_2sp(
          s[1],
          s[2],
          lat_origin=s[3],
          lon_origin=s[4],
          stdlat1=s[5],
          stdlat2=s[6],
          false_easting=fe,
          false_northing=fn_,
        )
      1 =>
        @conic.ConicGrid::lambert_1sp(
          s[1],
          s[2],
          lat_origin=s[3],
          lon_origin=s[4],
          k0=s[7],
          false_easting=fe,
          false_northing=fn_,
        )
      _ =>
        @conic.ConicGrid::albers(
          s[1],
          s[2],
          lat_origin=s[3],
          lon_origin=s[4],
          stdlat1=s[5],
          stdlat2=s[6],
          false_easting=fe,
          false_northing=fn_,
        )
    }
  } catch {
    e => abort("\{e.to_repr()}")
  }
}

///|
fn maxv(a : Double, b : Double) -> Double {
  if b.is_nan() || b > a {
    b
  } else {
    a
  }
}

///|
/// Horizontal distance (m, approximate) between two geographic points.
fn ground_err(
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> Double {
  let mut dl = (lon1 - lon2).abs()
  if dl > 180.0 {
    dl = 360.0 - dl
  }
  let dy = (lat1 - lat2).abs() * 111320.0
  let dx = dl * 111320.0 * @math.cos(lat1 * deg)
  maxv(dx, dy)
}

///|
test "matches a 50-digit evaluation of the EPSG / Snyder formulas" {
  // tools/gen_conic.py: 18 projections (LCC 2SP and 1SP, Albers; China
  // CGCS2000 and Krasovsky, US, Europe, Australia, southern hemisphere,
  // parallels 1e-4 degrees apart, polar, near the equator, one parallel,
  // a sphere), 1792 points from 89°S to 89°N and up to 170° from the
  // central meridian.
  let mut dxy = 0.0
  let mut dgamma = 0.0
  let mut dk = 0.0
  let mut rt = 0.0
  let mut dkr = 0.0
  for r in conic_mp {
    let g = grid_of(conic_setups[r[0].to_int()])
    let p = g.forward(r[1], r[2])
    dxy = maxv(dxy, maxv((p.easting - r[3]).abs(), (p.northing - r[4]).abs()))
    dgamma = maxv(dgamma, (p.gamma - r[5]).abs())
    dk = maxv(dk, (p.k - r[6]).abs() / r[6])
    let b = g.reverse(r[3], r[4])
    rt = maxv(rt, ground_err(r[1], r[2], b.lat, b.lon))
    // The scale from reverse equals the scale from forward.
    dkr = maxv(dkr, (b.k - r[6]).abs() / r[6])
  }
  println(
    "conic vs 50-digit: max |dx|,|dy| = \{dxy} m, max |dgamma| = \{dgamma} deg, max |dk|/k = \{dk}, reverse of the exact x, y: \{rt} m, k from reverse: \{dkr}",
  )
  assert_true(dxy < 2.5e-8)
  assert_true(dgamma < 1.0e-12)
  assert_true(dk < 1.0e-13)
  assert_true(rt < 2.5e-7)
  assert_true(dkr < 2.0e-12)
}

///|
test "agrees with PROJ" {
  let mut dxy = 0.0
  let mut dxy_close = 0.0
  let mut dgamma = 0.0
  let mut dk = 0.0
  let mut dh = 0.0
  let mut inv = 0.0
  for r in conic_proj {
    let i = r[0].to_int()
    let s = conic_setups[i]
    let g = grid_of(s)
    let p = g.forward(r[1], r[2])
    let d = maxv((p.easting - r[3]).abs(), (p.northing - r[4]).abs())
    // Setup 10 has standard parallels 40 and 40.0001 degrees: PROJ's
    // n = log(m1/m2) / log(t1/t2) loses about 8 digits to cancellation.
    // The 50-digit reference above shows our values are right there.
    if i == 10 {
      dxy_close = maxv(dxy_close, d)
    } else {
      dxy = maxv(dxy, d)
    }
    dgamma = maxv(dgamma, (p.gamma - r[5]).abs())
    dk = maxv(dk, (p.k - r[6]).abs())
    // Conformal: the meridional scale equals k. Equal area: it is 1 / k.
    let h = if s[0].to_int() == 2 { 1.0 / p.k } else { p.k }
    dh = maxv(dh, (h - r[7]).abs() / r[7])
    inv = maxv(inv, ground_err(r[1], r[2], r[8], r[9]))
  }
  println(
    "conic vs PROJ: max |dx|,|dy| = \{dxy} m (close parallels: \{dxy_close} m), max |dgamma| = \{dgamma} deg, max |dk| = \{dk}, max |dh|/h = \{dh}; PROJ's own inverse is within \{inv} m",
  )
  assert_true(dxy < 2.5e-7)
  assert_true(dxy_close < 0.01)
  // PROJ's convergence and scale factors come from numerical derivatives.
  assert_true(dgamma < 5.0e-8)
  assert_true(dk < 1.0e-6)
  assert_true(dh < 1.0e-5)
}

///|
test "IOGP GIGS 5102, 5103, 5109" {
  // The GIGS tolerances are 0.03 m (LCC) and 0.05 m (Albers) for
  // coordinates, 3e-7 and 6e-7 degree for latitude and longitude.
  for i, s in gigs_sets {
    let g = grid_of(s, fe=s[8], fn_=s[9])
    let lu = s[10]
    let au = s[11]
    let mut de = 0.0
    let mut dl = 0.0
    for r in gigs_points {
      if r[0].to_int() != i {
        continue
      }
      let p = g.forward(r[2] * au, r[3] * au)
      de = maxv(
        de,
        maxv((p.easting - r[4] * lu).abs(), (p.northing - r[5] * lu).abs()),
      )
      let b = g.reverse(r[4] * lu, r[5] * lu)
      dl = maxv(dl, maxv((b.lat - r[2] * au).abs(), (b.lon - r[3] * au).abs()))
    }
    println(
      "GIGS set \{i}: max |dE|,|dN| = \{de} m, max |dlat|,|dlon| = \{dl} deg",
    )
    // The files give 2 or 3 decimals. The residuals (up to 0.8 mm for
    // Lambert, 4.5 mm on GIGS-5109-09 for Albers) are a little above that
    // rounding, but PROJ shows the same differences at the same points, so
    // they are in the GIGS values. GIGS itself allows 3-5 cm.
    let (etol, ltol) = if s[0].to_int() == 2 {
      (0.005, 7.0e-8)
    } else {
      (0.001, 1.0e-8)
    }
    assert_true(de < etol)
    assert_true(dl < ltol)
  }
}

///|
test "EPSG Guidance Note 7-2 worked examples" {
  // Clarke 1866: a = 6378206.4 m, b = 6356583.8 m.
  let a = 6378206.4
  let f = 1.0 - 6356583.8 / a
  let ftus = 1200.0 / 3937.0
  // NAD27 / Texas South Central (EPSG 9802), units US survey feet.
  let texas = @conic.ConicGrid::lambert_2sp(
    a,
    f,
    lat_origin=27.0 + 50.0 / 60.0,
    lon_origin=-99.0,
    stdlat1=28.0 + 23.0 / 60.0,
    stdlat2=30.0 + 17.0 / 60.0,
    false_easting=2000000.0 * ftus,
  )
  let p = texas.forward(28.5, -96.0)
  assert_true((p.easting / ftus - 2963503.91).abs() < 0.005)
  assert_true((p.northing / ftus - 254759.80).abs() < 0.005)
  let b = texas.reverse(2963503.91 * ftus, 254759.80 * ftus)
  assert_true((b.lat - 28.5).abs() < 2.0e-8 && (b.lon - -96.0).abs() < 2.0e-8)
  // JAD69 / Jamaica National Grid (EPSG 9801).
  let jamaica = @conic.ConicGrid::lambert_1sp(
    a,
    f,
    lat_origin=18.0,
    lon_origin=-77.0,
    false_easting=250000.0,
    false_northing=150000.0,
  )
  let lat = 17.0 + 55.0 / 60.0 + 55.80 / 3600.0
  let lon = -(76.0 + 56.0 / 60.0 + 37.26 / 3600.0)
  let q = jamaica.forward(lat, lon)
  assert_true((q.easting - 255966.58).abs() < 0.005)
  assert_true((q.northing - 142493.51).abs() < 0.005)
  let c = jamaica.reverse(255966.58, 142493.51)
  assert_true((c.lat - lat).abs() < 1.0e-7 && (c.lon - lon).abs() < 1.0e-7)
}

///|
test "Snyder (USGS PP 1395) numerical examples" {
  // Appendix A: Clarke 1866, standard parallels 29.5/45.5 (Albers) and
  // 33/45 (Lambert), origin 23°N 96°W, point 35°N 75°W. Snyder carries
  // about 8 significant digits by hand, so agree to 0.1 m.
  let a = 6378206.4
  let f = 1.0 - (1.0 - 0.00676866).sqrt()
  let albers = @conic.ConicGrid::albers(
    a,
    f,
    lat_origin=23.0,
    lon_origin=-96.0,
    stdlat1=29.5,
    stdlat2=45.5,
  )
  let p = albers.forward(35.0, -75.0)
  assert_true((p.easting - 1885472.7).abs() < 0.1)
  assert_true((p.northing - 1535925.0).abs() < 0.1)
  let lcc = @conic.ConicGrid::lambert_2sp(
    a,
    f,
    lat_origin=23.0,
    lon_origin=-96.0,
    stdlat1=33.0,
    stdlat2=45.0,
  )
  let q = lcc.forward(35.0, -75.0)
  assert_true((q.easting - 1894410.9).abs() < 0.1)
  assert_true((q.northing - 1564649.5).abs() < 0.1)
}

///|
test "Lambert is conformal: measured scale is k in every direction" {
  let g = @geodesic.Geodesic::wgs84()
  let lcc = @conic.ConicGrid::lambert_2sp(
    6378137.0,
    1.0 / 298.257222101,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  for pt in [(20.0, 100.0), (36.0, 116.0), (50.0, 80.0), (30.0, 125.0)] {
    let p = lcc.forward(pt.0, pt.1)
    // Steps of 1e-4 degree (about 11 m) north, east, and north-east.
    for d in [(1.0e-4, 0.0), (0.0, 1.0e-4), (1.0e-4, 1.0e-4)] {
      let q = lcc.forward(pt.0 + d.0, pt.1 + d.1)
      let grid = @math.hypot(q.easting - p.easting, q.northing - p.northing)
      let ground = g.inverse(pt.0, pt.1, pt.0 + d.0, pt.1 + d.1).s12
      // Scale at the midpoint, to second order.
      let mid = lcc.forward(pt.0 + d.0 / 2.0, pt.1 + d.1 / 2.0)
      assert_true((grid / ground - mid.k).abs() < 1.0e-9)
    }
  }
}

///|
test "Albers is equal-area: projected areas equal ellipsoidal areas" {
  let g = @geodesic.Geodesic::wgs84()
  let aea = @conic.ConicGrid::albers(
    6378137.0,
    1.0 / 298.257223563,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  // A 30° x 50° box over China, edges sampled every 0.01°: the geodesic
  // polygon and the straight-edged polygon in the plane then differ by far
  // less than the bound below.
  let pts : Array[(Double, Double)] = []
  let n = 5000
  for i in 0..<n {
    pts.push((20.0, 80.0 + 50.0 * i.to_double() / n.to_double()))
  }
  for i in 0..<3000 {
    pts.push((20.0 + 30.0 * i.to_double() / 3000.0, 130.0))
  }
  for i in 0..<n {
    pts.push((50.0, 130.0 - 50.0 * i.to_double() / n.to_double()))
  }
  for i in 0..<3000 {
    pts.push((50.0 - 30.0 * i.to_double() / 3000.0, 80.0))
  }
  let ellipsoidal = g.polygon_area(pts).area.abs()
  let mut planar = 0.0
  let m = pts.length()
  for i in 0..<m {
    let p = aea.forward(pts[i].0, pts[i].1)
    let q = aea.forward(pts[(i + 1) % m].0, pts[(i + 1) % m].1)
    planar = planar + (p.easting * q.northing - q.easting * p.northing)
  }
  planar = planar.abs() / 2.0
  println(
    "Albers area check: ellipsoidal \{ellipsoidal} m^2, planar \{planar} m^2, relative difference \{(planar - ellipsoidal).abs() / ellipsoidal}",
  )
  assert_true((planar - ellipsoidal).abs() / ellipsoidal < 1.0e-6)
  // k is the parallel scale and 1/k the meridian scale: their product is 1.
  for lat in [-60.0, 0.0, 25.0, 36.0, 47.0, 80.0] {
    let p = aea.forward(lat, 120.0)
    let dn = aea.forward(lat + 1.0e-3, 120.0)
    let ground = g.inverse(lat, 120.0, lat + 1.0e-3, 120.0).s12
    let h = @math.hypot(dn.easting - p.easting, dn.northing - p.northing) /
      ground
    assert_true((h * aea.forward(lat + 0.5e-3, 120.0).k - 1.0).abs() < 1.0e-8)
  }
}

///|
test "scale is 1 on the standard parallels and smallest at the origin latitude" {
  let lcc = @conic.LambertConformalConic::new(
    6378137.0,
    1.0 / 298.257222101,
    25.0,
    47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  let aea = @conic.AlbersEqualArea::new(
    6378137.0,
    1.0 / 298.257222101,
    25.0,
    47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  for lat in [25.0, 47.0] {
    assert_true((lcc.forward(105.0, lat, 110.0).k - 1.0).abs() < 1.0e-14)
    assert_true((aea.forward(105.0, lat, 110.0).k - 1.0).abs() < 1.0e-14)
  }
  let lat0 = lcc.origin_latitude()
  assert_true(lat0 > 25.0 && lat0 < 47.0)
  let k0 = lcc.central_scale()
  assert_true((lcc.forward(105.0, lat0, 105.0).k - k0).abs() < 1.0e-15)
  for d in [-1.0, -0.01, 0.01, 1.0] {
    assert_true(lcc.forward(105.0, lat0 + d, 105.0).k > k0)
  }
  // The cone constant is sin(origin latitude).
  assert_true((lcc.cone_constant() - @math.sin(lat0 * deg)).abs() < 1.0e-15)
  assert_true(aea.origin_latitude() > 25.0 && aea.origin_latitude() < 47.0)
}

///|
test "limits: Mercator, polar stereographic, and close parallels" {
  let a = 6378137.0
  let f = 1.0 / 298.257223563
  // LCC with the standard parallel on the equator is Mercator.
  let merc = @conic.LambertConformalConic::new_1sp(a, f, 0.0) catch {
    e => abort("\{e.to_repr()}")
  }
  let p = merc.forward(0.0, 40.0, 30.0)
  assert_true((p.x - a * 30.0 * deg).abs() < 1.0e-8)
  assert_eq(p.gamma, 0.0)
  // Mercator northing: a * asinh(tan chi), chi the conformal latitude.
  let es = (f * (2.0 - f)).sqrt()
  let t = @math.tan(40.0 * deg)
  let tchi = t *
    @math.hypot(1.0, @math.sinh(es * @math.atanh(es * @math.sin(40.0 * deg)))) -
    @math.sinh(es * @math.atanh(es * @math.sin(40.0 * deg))) *
    @math.hypot(1.0, t)
  assert_true((p.y - a * @math.asinh(tchi)).abs() < 1.0e-7)
  // LCC with the standard parallel at the pole is polar stereographic:
  // the pole maps to the origin and parallels are circles.
  let ps = @conic.LambertConformalConic::new_1sp(a, f, 90.0, k0=0.994) catch {
    e => abort("\{e.to_repr()}")
  }
  let pole = ps.forward(0.0, 90.0, 0.0)
  assert_true(pole.x.abs() < 1.0e-9 && pole.y.abs() < 1.0e-9)
  assert_true((pole.k - 0.994).abs() < 1.0e-15)
  let r1 = ps.forward(0.0, 80.0, 0.0)
  let r2 = ps.forward(0.0, 80.0, 123.0)
  assert_true(
    (@math.hypot(r1.x, r1.y) - @math.hypot(r2.x, r2.y)).abs() < 1.0e-8,
  )
  // Same as this library's UPS (polar stereographic, k0 = 0.994).
  let ups = @utm.PolarStereographic::ups()
  for pt in [(80.0, 0.0), (84.5, 123.0), (60.0, -45.0), (89.99, 179.0)] {
    let c = ps.forward(0.0, pt.0, pt.1)
    let u = ups.forward(true, pt.0, pt.1)
    assert_true((c.x - u.x).abs() < 1.0e-8 && (c.y - u.y).abs() < 1.0e-8)
    assert_true((c.k - u.k).abs() < 1.0e-14)
  }
  // Parallels 1e-9 degrees apart define a slightly different projection
  // from the one-parallel one. The difference, a few micrometres, is
  // computed here from two nearly equal quantities; it must match a 50-digit
  // evaluation (tools/gen_conic.py formulas), which naive double-precision
  // formulas cannot do. Both grids use the false origin 40°N 0°E.
  let one = @conic.ConicGrid::lambert_1sp(a, f, lat_origin=40.0, lon_origin=0.0) catch {
    e => abort("\{e.to_repr()}")
  }
  let two = @conic.ConicGrid::lambert_2sp(
    a,
    f,
    lat_origin=40.0,
    lon_origin=0.0,
    stdlat1=40.0,
    stdlat2=40.0 + 1.0e-9,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  for
    c in [
      (55.0, 20.0, 3.66808904812e-6, 9.65725533851e-7),
      (10.0, -60.0, 1.57930680721e-5, -2.54672644642e-5),
      (40.0, 0.0, 0.0, 0.0),
    ] {
    let u = one.forward(c.0, c.1)
    let v = two.forward(c.0, c.1)
    assert_true(
      (u.easting - v.easting - c.2).abs() < 2.0e-8 &&
      (u.northing - v.northing - c.3).abs() < 2.0e-8,
    )
  }
}

///|
test "longitudes wrap and reverse normalizes" {
  let g = @conic.ConicGrid::albers(
    6378137.0,
    1.0 / 298.257223563,
    lat_origin=0.0,
    lon_origin=170.0,
    stdlat1=20.0,
    stdlat2=50.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  let a = g.forward(40.0, -175.0)
  let b = g.forward(40.0, 185.0)
  assert_true((a.easting - b.easting).abs() < 1.0e-9)
  let back = g.reverse(a.easting, a.northing)
  assert_true((back.lon - -175.0).abs() < 1.0e-9)
}

///|
fn conic_err(f : () -> Unit raise @conic.ConicError) -> @conic.ConicError? {
  try f() catch {
    e => Some(e)
  } noraise {
    _ => None
  }
}

///|
test "invalid parameters raise ConicError" {
  let a = 6378137.0
  let f = 1.0 / 298.257223563
  assert_true(
    conic_err(() => {
      ignore(@conic.LambertConformalConic::new(0.0, f, 30.0, 60.0))
    })
    is Some(@conic.InvalidEllipsoid(_, _)),
  )
  assert_true(
    conic_err(() => ignore(@conic.AlbersEqualArea::new(a, 1.0, 30.0, 60.0)))
    is Some(@conic.InvalidEllipsoid(_, _)),
  )
  assert_true(
    conic_err(() => ignore(@conic.LambertConformalConic::new(a, f, 91.0, 60.0)))
    is Some(@conic.InvalidStandardLatitude(_)),
  )
  assert_true(
    conic_err(() => {
      ignore(@conic.AlbersEqualArea::new(a, f, 30.0, @double.not_a_number))
    })
    is Some(@conic.InvalidStandardLatitude(_)),
  )
  assert_true(
    conic_err(() => ignore(@conic.LambertConformalConic::new(a, f, 90.0, 60.0)))
    is Some(@conic.IncompatibleStandardLatitudes(_, _)),
  )
  assert_true(
    conic_err(() => ignore(@conic.AlbersEqualArea::new(a, f, 90.0, -90.0)))
    is Some(@conic.IncompatibleStandardLatitudes(_, _)),
  )
  assert_true(
    conic_err(() => {
      ignore(@conic.LambertConformalConic::new_1sp(a, f, 40.0, k0=0.0))
    })
    is Some(@conic.InvalidScale(_)),
  )
  assert_true(
    conic_err(() => {
      ignore(
        @conic.ConicGrid::albers(
          a,
          f,
          lat_origin=95.0,
          lon_origin=0.0,
          stdlat1=20.0,
          stdlat2=50.0,
        ),
      )
    })
    is Some(@conic.InvalidOriginLatitude(_)),
  )
  // A sphere and a prolate ellipsoid are accepted.
  ignore(@conic.AlbersEqualArea::new(6371000.0, 0.0, 30.0, 60.0))
  ignore(@conic.LambertConformalConic::new(6371000.0, -0.01, 30.0, 60.0))
}

///|
test "Albers: poles round-trip exactly; near-pole points stay accurate" {
  let a = 6378137.0
  let f = 1.0 / 298.257223563
  for
    s in [
      (25.0, 47.0, 105.0),
      (-18.0, -36.0, 132.0),
      (30.0, -30.0, 0.0),
      (60.0, 85.0, 0.0),
      (89.0, 89.0, 0.0),
      (90.0, 90.0, 0.0),
      (-90.0, -90.0, 0.0),
    ] {
    let p = @conic.AlbersEqualArea::new(a, f, s.0, s.1) catch {
      e => abort("\{e.to_repr()}")
    }
    for pole in [90.0, -90.0] {
      // The pole of a polar projection's far side is the outer circle; skip
      // the azimuthal case's antipodal pole, which is a whole circle.
      if (s.0 == 90.0 && pole < 0.0) || (s.0 == -90.0 && pole > 0.0) {
        continue
      }
      for lon in [s.2, s.2 + 30.0, s.2 - 150.0] {
        let fw = p.forward(s.2, pole, lon)
        let back = p.reverse(s.2, fw.x, fw.y)
        assert_eq(back.lat, pole)
        // 1e-4 degree (11 m) from the pole the map still resolves latitude.
        let near = if pole > 0.0 { pole - 1.0e-4 } else { pole + 1.0e-4 }
        let g = p.forward(s.2, near, lon)
        let b = p.reverse(s.2, g.x, g.y)
        assert_true((b.lat - near).abs() < 1.0e-7)
        // Within about 1e-6 degree of the pole the answer is the pole or
        // close to it: never worse than 4e-6 degree (45 cm).
        let close = if pole > 0.0 { pole - 1.0e-6 } else { pole + 1.0e-6 }
        let c = p.forward(s.2, close, lon)
        let d = p.reverse(s.2, c.x, c.y)
        assert_true((d.lat - close).abs() < 4.0e-6)
      }
    }
  }
}

///|
test "Lambert: both poles round-trip; the far pole is at a huge distance" {
  let lcc = @conic.LambertConformalConic::new(
    6378137.0,
    1.0 / 298.257223563,
    25.0,
    47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  let apex = lcc.forward(105.0, 90.0, 130.0)
  assert_eq(lcc.reverse(105.0, apex.x, apex.y).lat, 90.0)
  let far = lcc.forward(105.0, -90.0, 130.0)
  assert_true(@math.hypot(far.x, far.y) > 1.0e20)
  assert_eq(lcc.reverse(105.0, far.x, far.y).lat, -90.0)
}

///|
test "prolate ellipsoid: conformal and equal-area, measured with geodesics" {
  // No reference implementation handles f < 0 here (PROJ does not), so
  // check the defining properties with this library's geodesics, which are
  // tested on prolate ellipsoids against GeodSolve.
  let a = 6378137.0
  let f = -0.01
  let g = @geodesic.Geodesic::new(a, f)
  let lcc = @conic.LambertConformalConic::new(a, f, 25.0, 47.0)
  let aea = @conic.AlbersEqualArea::new(a, f, 25.0, 47.0)
  for pt in [(20.0, 100.0), (36.0, 116.0), (50.0, 80.0)] {
    let p = lcc.forward(105.0, pt.0, pt.1)
    for d in [(1.0e-4, 0.0), (0.0, 1.0e-4), (1.0e-4, 1.0e-4)] {
      let q = lcc.forward(105.0, pt.0 + d.0, pt.1 + d.1)
      let grid = @math.hypot(q.x - p.x, q.y - p.y)
      let ground = g.inverse(pt.0, pt.1, pt.0 + d.0, pt.1 + d.1).s12
      let mid = lcc.forward(105.0, pt.0 + d.0 / 2.0, pt.1 + d.1 / 2.0)
      assert_true((grid / ground - mid.k).abs() < 1.0e-9)
    }
    let back = lcc.reverse(105.0, p.x, p.y)
    assert_true((back.lat - pt.0).abs() < 1.0e-12)
    let r = aea.forward(105.0, pt.0, pt.1)
    let rb = aea.reverse(105.0, r.x, r.y)
    assert_true(
      (rb.lat - pt.0).abs() < 1.0e-12 && (rb.lon - pt.1).abs() < 1.0e-12,
    )
  }
  // Equal area: a 10° x 10° box, edges sampled every 0.01°.
  let ring : Array[(Double, Double)] = []
  for i in 0..<1000 {
    ring.push((30.0, 100.0 + i.to_double() * 0.01))
  }
  for i in 0..<1000 {
    ring.push((30.0 + i.to_double() * 0.01, 110.0))
  }
  for i in 0..<1000 {
    ring.push((40.0, 110.0 - i.to_double() * 0.01))
  }
  for i in 0..<1000 {
    ring.push((40.0 - i.to_double() * 0.01, 100.0))
  }
  let exact = g.polygon_area(ring).area.abs()
  let mut s = 0.0
  for i in 0..<ring.length() {
    let p = aea.forward(105.0, ring[i].0, ring[i].1)
    let q = aea.forward(
      105.0,
      ring[(i + 1) % ring.length()].0,
      ring[(i + 1) % ring.length()].1,
    )
    s = s + (p.x * q.y - q.x * p.y)
  }
  assert_true(((s.abs() / 2.0 - exact) / exact).abs() < 1.0e-6)
}

///|
test "latitudes beyond 90 degrees give NaN" {
  let a = 6378137.0
  let f = 1.0 / 298.257223563
  let lcc = @conic.LambertConformalConic::new(a, f, 25.0, 47.0)
  let aea = @conic.AlbersEqualArea::new(a, f, 25.0, 47.0)
  for lat in [90.5, -91.0, @double.not_a_number] {
    let p = lcc.forward(105.0, lat, 110.0)
    let q = aea.forward(105.0, lat, 110.0)
    assert_true(p.x.is_nan() && p.y.is_nan() && q.x.is_nan() && q.y.is_nan())
  }
}
