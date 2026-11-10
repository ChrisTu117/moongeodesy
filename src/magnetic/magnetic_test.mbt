///|
fn max_abs(errs : Array[Double]) -> Double {
  let mut m = 0.0
  for e in errs {
    if e.is_nan() || e.abs() > m {
      m = e.abs()
    }
  }
  m
}

///|
/// Difference of two angles in degrees, reduced to [-180, 180].
fn ang_err(a : Double, b : Double) -> Double {
  let mut d = (a - b) % 360.0
  if d > 180.0 {
    d = d - 360.0
  } else if d < -180.0 {
    d = d + 360.0
  }
  d.abs()
}

///|
test "WMM2025 matches the 100 official double-precision test values" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let dyz = [] // Y, Z (nT)
  let dxhf = [] // X, H, F (nT)
  let drate = [] // dX/dt ... dF/dt (nT/yr)
  let dang = [] // D, I (deg), printed to 2 decimals
  let dangrate = [] // dD/dt, dI/dt (deg/yr)
  for row in wmm2025_precise {
    let e = wmm.field(row[2], row[3], row[1] * 1000.0, row[0])
    dxhf.push(e.h - row[6])
    dxhf.push(e.x - row[7])
    dyz.push(e.y - row[8])
    dyz.push(e.z - row[9])
    dxhf.push(e.f - row[10])
    dang.push(ang_err(e.decl, row[4]))
    dang.push(e.incl - row[5])
    dangrate.push(e.decl_dot - row[11])
    dangrate.push(e.incl_dot - row[12])
    drate.push(e.hdot - row[13])
    drate.push(e.xdot - row[14])
    drate.push(e.ydot - row[15])
    drate.push(e.zdot - row[16])
    drate.push(e.fdot - row[17])
  }
  println(
    "WMM2025 vs 100 official values: max |dY|,|dZ| = \{max_abs(dyz)} nT, max |dX|,|dH|,|dF| = \{max_abs(dxhf)} nT, max |d rate| = \{max_abs(drate)} nT/yr, max |dD|,|dI| = \{max_abs(dang)} deg, max |dD/dt|,|dI/dt| = \{max_abs(dangrate)} deg/yr",
  )
  // The file has 6 decimals (half a unit is 5e-7). Y and Z agree to that
  // rounding plus a few units of double rounding.
  assert_true(max_abs(dyz) < 2.5e-6)
  // X (and so H and F) differ by up to 7.2e-4 nT, 1.3e-8 of the field. The
  // independent 40-digit evaluation below reproduces our X, not the file's,
  // and the file agrees with both on the equator, so the difference is in
  // the file. It is far below the 0.1 nT the report quotes for
  // single-precision software and the model's ~100 nT accuracy.
  assert_true(max_abs(dxhf) < 1.0e-3)
  assert_true(max_abs(drate) < 2.0e-6)
  assert_true(max_abs(dangrate) <= 5.1e-7)
  // D and I are printed with 2 decimals.
  assert_true(max_abs(dang) <= 0.005 + 1.0e-9)
}

///|
test "WMM2025 matches an independent 40-digit evaluation (119 points)" {
  // tools/gen_wmm_independent.py: exact Legendre polynomials, numerical
  // gradient of the potential in ECEF, projection on the geodetic frame.
  // Includes both poles on three meridians, -1 km and 850 km, 2025.0 and
  // 2030.0, and nine Chinese cities.
  let wmm = @magnetic.MagneticModel::wmm2025()
  let d = []
  let drate = []
  for row in wmm2025_independent {
    let e = wmm.field(row[2], row[3], row[1] * 1000.0, row[0])
    d.push(e.x - row[4])
    d.push(e.y - row[5])
    d.push(e.z - row[6])
    drate.push(e.xdot - row[7])
    drate.push(e.ydot - row[8])
    drate.push(e.zdot - row[9])
  }
  println(
    "WMM2025 vs independent evaluation: max |dX|,|dY|,|dZ| = \{max_abs(d)} nT, max |d rate| = \{max_abs(drate)} nT/yr",
  )
  assert_true(max_abs(d) < 1.0e-9)
  assert_true(max_abs(drate) < 1.0e-12)
}

///|
test "WMM2025 matches report Table 6 (rounded) including grid variation" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let dnt = []
  let ddeg = []
  for row in wmm2025_table6 {
    let e = wmm.field(row[2], row[3], row[1] * 1000.0, row[0])
    for
      pair in [
        (e.x, row[4]),
        (e.y, row[5]),
        (e.z, row[6]),
        (e.h, row[7]),
        (e.f, row[8]),
        (e.xdot, row[12]),
        (e.ydot, row[13]),
        (e.zdot, row[14]),
        (e.hdot, row[15]),
        (e.fdot, row[16]),
      ] {
      dnt.push(pair.0 - pair.1)
    }
    ddeg.push(e.incl - row[9])
    ddeg.push(ang_err(e.decl, row[10]))
    if row[11].is_nan() {
      assert_true(e.gv.is_nan())
    } else {
      ddeg.push(ang_err(e.gv, row[11]))
    }
    ddeg.push(e.incl_dot - row[17])
    ddeg.push(e.decl_dot - row[18])
  }
  // 0.1 nT and 0.01 degree rounding.
  assert_true(max_abs(dnt) <= 0.05 + 1.0e-9)
  assert_true(max_abs(ddeg) <= 0.005 + 1.0e-9)
}

///|
test "report Table 3b: high-precision numerical example" {
  // 2027.5, 100 km, 80°S, 240°E. Values copied from Table 3b of the
  // WMM2025 technical report (10 significant decimals).
  let e = @magnetic.MagneticModel::wmm2025().field(
    -80.0, 240.0, 100000.0, 2027.5,
  )
  let rad = 180.0 / 3.141592653589793
  let checks = [
    (e.x, 5983.9760496518),
    (e.y, 14760.1359757868),
    (e.z, -49317.6706154255),
    (e.xdot, 30.5553533530),
    (e.ydot, -8.0494228995),
    (e.zdot, 89.2174069382),
    (e.f, 51825.6907172314),
    (e.h, 15927.0079860130),
    (e.fdot, -83.6643506802),
    (e.hdot, 4.0203361601),
  ]
  // Half a unit in the 10th decimal.
  for c in checks {
    assert_true((c.0 - c.1).abs() <= 5.1e-11)
  }
  // Angles are given in radians to 10 decimals (half a unit is 5e-11 rad).
  assert_true((e.decl / rad - 1.1856308407).abs() <= 5.1e-11)
  assert_true((e.incl / rad - -1.2584221541).abs() <= 5.1e-11)
  assert_true((e.decl_dot / rad - -0.0019677910).abs() <= 5.1e-11)
  assert_true((e.incl_dot / rad - 0.0006028663).abs() <= 5.1e-11)
}

///|
test "report Table 3d: grivation on UPS and UTM grids from our own utm package" {
  // Table 3d: TrueN-to-MagN 67.9316432316°; GridN-to-MagN for UPS north,
  // UPS south, UTM zone 10 and zone 11.
  let lat = -80.0
  let lon = 240.0
  let e = @magnetic.MagneticModel::wmm2025().field(lat, lon, 100000.0, 2027.5)
  assert_true((e.decl - 67.9316432316).abs() < 1.0e-9)
  // UPS: convergence is +lon (north grid) or -lon (south grid).
  assert_true(
    ang_err(@magnetic.grivation(e.decl, lon), -172.0683567684) < 1.0e-9,
  )
  assert_true(
    ang_err(@magnetic.grivation(e.decl, -lon), 307.9316432316) < 1.0e-9,
  )
  // The WMM polar grid variation is the UPS south value here.
  assert_true(ang_err(e.gv, 307.9316432316) < 1.0e-9)
  // UTM zones 10 and 11: convergence from our transverse Mercator.
  let z10 = @utm.from_latlon(lat, lon - 360.0, zone=10)
  let z11 = @utm.from_latlon(lat, lon - 360.0, zone=11)
  assert_true((z10.gamma - -2.9545046801).abs() < 1.0e-9)
  assert_true((z11.gamma - 2.9545046801).abs() < 1.0e-9)
  assert_true(
    ang_err(@magnetic.grivation(e.decl, z10.gamma), 70.8861479117) < 1.0e-9,
  )
  assert_true(
    ang_err(@magnetic.grivation(e.decl, z11.gamma), 64.9771385515) < 1.0e-9,
  )
}

///|
test "field is finite and continuous at the geographic poles" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  for lat in [90.0, -90.0] {
    for lon in [0.0, 90.0, -135.0] {
      let p = wmm.field(lat, lon, 0.0, 2026.0)
      let near = wmm.field(lat - lat.signum() * 1.0e-7, lon, 0.0, 2026.0)
      for v in [p.x, p.y, p.z, p.h, p.f, p.decl, p.incl, p.xdot, p.ydot, p.zdot] {
        assert_true(!v.is_nan() && !v.is_inf())
      }
      // 1e-7 degree is about 1 cm; the field changes by far less than 1e-3 nT.
      assert_true((p.x - near.x).abs() < 1.0e-3)
      assert_true((p.y - near.y).abs() < 1.0e-3)
      assert_true((p.z - near.z).abs() < 1.0e-3)
    }
  }
  // At a pole the vector does not depend on longitude, only its north/east
  // split does: H, Z, and F are the same for every meridian.
  let a = wmm.field(90.0, 0.0, 0.0, 2026.0)
  let b = wmm.field(90.0, 123.0, 0.0, 2026.0)
  assert_true((a.h - b.h).abs() < 1.0e-6)
  assert_true((a.z - b.z).abs() < 1.0e-6)
  // Rotating the meridian by 123 degrees rotates D by the same angle.
  assert_true(ang_err(b.decl - a.decl, 123.0) < 1.0e-9)
}

///|
test "longitude is periodic" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let a = wmm.field(31.2304, 121.4737, 0.0, 2026.5)
  let b = wmm.field(31.2304, 121.4737 - 360.0, 0.0, 2026.5)
  let c = wmm.field(31.2304, 121.4737 + 720.0, 0.0, 2026.5)
  for p in [(a.x, b.x), (a.y, b.y), (a.z, b.z), (a.x, c.x), (a.y, c.y)] {
    assert_true((p.0 - p.1).abs() < 1.0e-8)
  }
}

///|
test "rates agree with finite differences in time" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let dt = 0.25
  let a = wmm.field(39.9087, 116.3975, 50.0, 2026.0 - dt)
  let b = wmm.field(39.9087, 116.3975, 50.0, 2026.0 + dt)
  let m = wmm.field(39.9087, 116.3975, 50.0, 2026.0)
  // X, Y, Z are linear in time, so the central difference is exact.
  assert_true(((b.x - a.x) / (2.0 * dt) - m.xdot).abs() < 1.0e-8)
  assert_true(((b.y - a.y) / (2.0 * dt) - m.ydot).abs() < 1.0e-8)
  assert_true(((b.z - a.z) / (2.0 * dt) - m.zdot).abs() < 1.0e-8)
  // D, I, H, F are not linear; agree to second order.
  assert_true(((b.decl - a.decl) / (2.0 * dt) - m.decl_dot).abs() < 1.0e-5)
  assert_true(((b.incl - a.incl) / (2.0 * dt) - m.incl_dot).abs() < 1.0e-5)
  assert_true(((b.h - a.h) / (2.0 * dt) - m.hdot).abs() < 1.0e-3)
  assert_true(((b.f - a.f) / (2.0 * dt) - m.fdot).abs() < 1.0e-3)
}

///|
test "derived elements are consistent with X, Y, Z" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let e = wmm.field(22.5431, 114.0579, 0.0, 2026.0)
  assert_true((e.h * e.h - e.x * e.x - e.y * e.y).abs() / (e.h * e.h) < 1.0e-14)
  assert_true((e.f * e.f - e.h * e.h - e.z * e.z).abs() / (e.f * e.f) < 1.0e-14)
  let rad = 3.141592653589793 / 180.0
  assert_true(
    (e.h * @math.cos(e.decl * rad) - e.x).abs() < 1.0e-9 &&
    (e.h * @math.sin(e.decl * rad) - e.y).abs() < 1.0e-9,
  )
  assert_true((e.f * @math.sin(e.incl * rad) - e.z).abs() < 1.0e-9)
}

///|
fn field_err(
  wmm : @magnetic.MagneticModel,
  lat : Double,
  lon : Double,
  h : Double,
  t : Double,
) -> @magnetic.MagneticError? {
  try wmm.field(lat, lon, h, t) catch {
    e => Some(e)
  } noraise {
    _ => None
  }
}

///|
fn date_err(y : Int32, m : Int32, d : Int32) -> @magnetic.MagneticError? {
  try @magnetic.decimal_year(y, m, d) catch {
    e => Some(e)
  } noraise {
    _ => None
  }
}

///|
test "domain errors" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let bad : Array[(Double, Double, Double, Double)] = [
    (90.5, 0.0, 0.0, 2026.0),
    (@double.not_a_number, 0.0, 0.0, 2026.0),
    (0.0, @double.infinity, 0.0, 2026.0),
    (0.0, 0.0, -1000.5, 2026.0),
    (0.0, 0.0, 850000.5, 2026.0),
    (0.0, 0.0, 0.0, 2024.999),
    (0.0, 0.0, 0.0, 2030.001),
    (0.0, 0.0, @double.not_a_number, 2026.0),
    (0.0, 0.0, 0.0, @double.not_a_number),
  ]
  for b in bad {
    assert_true(field_err(wmm, b.0, b.1, b.2, b.3) is Some(_))
  }
  assert_true(
    field_err(wmm, 91.0, 0.0, 0.0, 2026.0) is Some(@magnetic.InvalidLatitude(_)),
  )
  assert_true(
    field_err(wmm, 0.0, 0.0, 0.0, 2031.0) is Some(@magnetic.DateOutOfRange(_)),
  )
  assert_true(
    field_err(wmm, 0.0, 0.0, 900000.0, 2026.0)
    is Some(@magnetic.HeightOutOfRange(_)),
  )
  assert_true(
    field_err(wmm, 0.0, @double.neg_infinity, 0.0, 2026.0)
    is Some(@magnetic.InvalidLongitude(_)),
  )
  // The ends of the domain are accepted.
  let _ = wmm.field(-90.0, 0.0, -1000.0, 2025.0)
  let _ = wmm.field(90.0, 0.0, 850000.0, 2030.0)
}

///|
test "decimal_year follows the WMM software convention" {
  assert_eq(@magnetic.decimal_year(2025, 1, 1), 2025.0)
  assert_eq(@magnetic.decimal_year(2025, 7, 2), 2025.0 + 182.0 / 365.0)
  assert_eq(@magnetic.decimal_year(2028, 12, 31), 2028.0 + 365.0 / 366.0)
  assert_eq(@magnetic.decimal_year(2028, 3, 1), 2028.0 + 60.0 / 366.0)
  assert_eq(@magnetic.decimal_year(2100, 3, 1), 2100.0 + 59.0 / 365.0)
  assert_eq(@magnetic.decimal_year(2000, 3, 1), 2000.0 + 60.0 / 366.0)
  for
    d in [
      (2025, 2, 29),
      (2025, 13, 1),
      (2025, 4, 31),
      (2025, 0, 1),
      (2025, 1, 0),
    ] {
    assert_true(date_err(d.0, d.1, d.2) is Some(@magnetic.InvalidDate(_, _, _)))
  }
}

///|
test "compass zones use the WMM2025 Blackout and Caution thresholds" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  // First official test point: H = 1504.3 nT near the north dip pole.
  let pole = wmm.field(89.0, -121.0, 28000.0, 2025.0)
  assert_eq(pole.zone(), @magnetic.Blackout)
  // H = 2164.3 nT.
  let caution = wmm.field(80.0, -96.0, 48000.0, 2025.0)
  assert_eq(caution.zone(), @magnetic.Caution)
  let beijing = wmm.field(39.9087, 116.3975, 0.0, 2026.0)
  assert_eq(beijing.zone(), @magnetic.Normal)
}

///|
test "bearing conversions" {
  // With D = -7.5 (west), true 0 shows as magnetic 7.5.
  assert_eq(@magnetic.magnetic_bearing(0.0, -7.5), 7.5)
  assert_eq(@magnetic.true_bearing(7.5, -7.5), 0.0)
  assert_eq(@magnetic.magnetic_bearing(355.0, 10.0), 345.0)
  assert_eq(@magnetic.true_bearing(355.0, 10.0), 5.0)
  for b in [0.0, 12.5, 180.0, 359.9] {
    for d in [-20.0, -0.1, 3.0, 179.0] {
      let m = @magnetic.magnetic_bearing(b, d)
      assert_true(m >= 0.0 && m < 360.0)
      assert_true(ang_err(@magnetic.true_bearing(m, d), b) < 1.0e-12)
    }
  }
  assert_eq(@magnetic.grivation(10.0, -175.0), -175.0)
  assert_eq(@magnetic.grivation(10.0, -170.0), -180.0)
}
