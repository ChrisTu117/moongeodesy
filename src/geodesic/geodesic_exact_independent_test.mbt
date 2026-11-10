// GeodesicExact against tools/gen_exact_independent.py: quadrature of the
// classical distance and longitude integrals in mpmath, not GeographicLib.

///|
/// Metres per radian along the meridian and along the parallel at `lat`.
fn exact_indep_radii(a : Double, f : Double, lat : Double) -> (Double, Double) {
  let e2 = f * (2.0 - f)
  let phi = lat * 3.141592653589793 / 180.0
  let w = 1.0 - e2 * @math.sin(phi) * @math.sin(phi)
  let m = a * (1.0 - e2) / (w * w.sqrt())
  let n = a / w.sqrt()
  (m, n * @math.cos(phi))
}

///|
test "GeodesicExact direct matches mpmath quadrature, f from -0.5 to 0.5" {
  let a = 6378137.0
  let mut dp = 0.0
  let mut da = 0.0
  let mut dS = 0.0
  let mut da_pole = 0.0
  let mut dS_pole = 0.0
  let mut n = 0
  for row in exact_indep {
    if row[1] != 1.0 {
      continue
    }
    let f = exact_indep_f[row[0].to_int()]
    let g = exact_geod(a, f)
    let r = g.direct(row[2], row[3], row[4], row[5])
    n = n + 1
    let (mr, pr) = exact_indep_radii(a, f, row[6])
    let rad = 3.141592653589793 / 180.0
    let pos = (r.lat2 - row[6]).abs() * rad * mr +
      indep_ang(r.lon2, row[7]) * rad * pr
    dp = indep_max(dp, pos)
    // Within 1 degree of a pole a nanometre of position is a large change
    // of longitude, which moves azi2 and S12 with it.
    let near_pole = @cmp.maximum(row[2].abs(), row[6].abs()) >= 89.0
    let dazi = indep_ang(r.azi2, row[8])
    let dS12 = if row[9].is_nan() {
      0.0
    } else {
      (r.area12 - row[9]).abs() / g.ellipsoid_area()
    }
    if near_pole {
      da_pole = indep_max(da_pole, dazi)
      dS_pole = indep_max(dS_pole, dS12)
    } else {
      da = indep_max(da, dazi)
      dS = indep_max(dS, dS12)
    }
  }
  println(
    "exact direct vs mpmath: n=\{n} max position err m=\{dp} max |dazi2| deg=\{da} max |dS12|/area=\{dS}; within 1 deg of a pole: |dazi2| deg=\{da_pole} |dS12|/area=\{dS_pole}",
  )
  assert_eq(n, 209)
  assert_true(dp < 2.5e-8) // README: position < 25 nm
  assert_true(da < 1.0e-12) // README: azi2 < 1e-12 deg
  assert_true(dS < 1.0e-15) // README: area < 1e-15 of the ellipsoid
  assert_true(da_pole < 1.0e-10) // README: near a pole, azi2 < 1e-10 deg
  assert_true(dS_pole < 1.0e-13) // README: and area < 1e-13
}

///|
test "GeodesicExact inverse matches mpmath quadrature, f from -0.5 to 0.5" {
  let a = 6378137.0
  let mut ds = 0.0
  let mut da = 0.0
  let mut dS = 0.0
  let mut n = 0
  for row in exact_indep {
    if row[1] != 0.0 {
      continue
    }
    let f = exact_indep_f[row[0].to_int()]
    let g = exact_geod(a, f)
    let r = g.inverse(row[2], row[3], row[4], row[5])
    n = n + 1
    ds = indep_max(ds, (r.s12 - row[6]).abs())
    // Azimuth errors as sideways offsets at the far end of the line.
    let rad = 3.141592653589793 / 180.0
    da = indep_max(da, indep_ang(r.azi1, row[7]) * rad * row[6])
    da = indep_max(da, indep_ang(r.azi2, row[8]) * rad * row[6])
    dS = indep_max(dS, (r.area12 - row[9]).abs() / g.ellipsoid_area())
  }
  println(
    "exact inverse vs mpmath: n=\{n} max |ds12| m=\{ds} max azimuth offset m=\{da} max |dS12|/area=\{dS}",
  )
  assert_eq(n, 165)
  assert_true(ds < 2.0e-8) // README: distance < 20 nm
  assert_true(da < 1.0e-8) // README: azimuth < 10 nm at the far end
  assert_true(dS < 1.0e-15)
}
