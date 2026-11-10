///|
test "report Table 3b intermediate quantities" {
  // 2027.5, 100 km, 80°S, 240°E; values from Table 3b of the WMM2025
  // technical report (10 decimals; half a unit is 5e-11).
  let wmm = MagneticModel::wmm2025()
  let (phip, r) = geocentric(-80.0, 100000.0)
  assert_true((phip - -1.3951289589).abs() <= 5.1e-11)
  // r is given to 10 decimals, below one unit in the last place of a double
  // (9.3e-10 m here), so allow two units.
  assert_true((r - 6457402.3484473705).abs() <= 1.9e-9)
  let (g, h) = wmm.coefficients_at(2027.5)
  let k = (n : Int32, m : Int32) => n * (n + 1) / 2 + m
  for
    c in [
      (g[k(1, 0)], -29321.8),
      (g[k(1, 1)], -1386.55),
      (g[k(2, 0)], -2585.6),
      (g[k(2, 1)], 2938.1),
      (g[k(2, 2)], 1629.3),
      (h[k(1, 1)], 4491.65),
      (h[k(2, 1)], -3202.85),
      (h[k(2, 2)], -845.35),
    ] {
    assert_true((c.0 - c.1).abs() <= 1.0e-9)
  }
  let sp = @math.sin(phip)
  let cp = @math.cos(phip)
  let (t, dt) = legendre_poly(wmm.nmax, sp)
  let ratio = wmm.radius / r
  let (xp, yp, zp) = wmm.sum_field(g, h, t, dt, sp, cp, ratio, 240.0)
  // The geocentric components agree within one unit of the 10th decimal
  // (X' differs by 5.4e-11, just over half a unit).
  assert_true((xp - 5928.0241392588).abs() <= 1.0e-10)
  assert_true((yp - 14760.1359757868).abs() <= 1.0e-10)
  assert_true((zp - -49324.4273570284).abs() <= 1.0e-10)
  let (xd, yd, zd) = wmm.sum_field(
    wmm.gdot,
    wmm.hdot,
    t,
    dt,
    sp,
    cp,
    ratio,
    240.0,
  )
  assert_true((xd - 30.6565457063).abs() <= 1.0e-10)
  assert_true((yd - -8.0494228995).abs() <= 1.0e-10)
  assert_true((zd - 89.1826862382).abs() <= 1.0e-10)
}
