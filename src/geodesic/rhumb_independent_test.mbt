// Rhumb lines against tools/gen_rhumb_independent.py: textbook formulas in
// mpmath (isometric latitude, meridian-arc quadrature, Green's theorem for
// the area), not GeographicLib.

///|
fn indep_rhumb(code : Double) -> @geodesic.Rhumb {
  if code == 0.0 {
    @geodesic.Rhumb::wgs84()
  } else {
    @geodesic.Rhumb::new(6378137.0, -1.0 / 298.257223563) catch {
      _ => abort("bad ellipsoid")
    }
  }
}

///|
/// Metres per degree of arc used to turn angle errors into distances.
let indep_m_per_deg : Double = 6378137.0 * 3.141592653589793 / 180.0

///|
/// Running maximum that keeps a NaN, so a NaN error fails the bound.
fn indep_max(m : Double, e : Double) -> Double {
  if m.is_nan() || e.is_nan() {
    @double.not_a_number
  } else if e > m {
    e
  } else {
    m
  }
}

///|
fn indep_ang(a : Double, b : Double) -> Double {
  let mut d = (a - b).abs()
  d = d - 360.0 * (d / 360.0).round()
  d.abs()
}

///|
test "rhumb inverse matches independent mpmath formulas" {
  let mut ds = 0.0
  let mut da = 0.0
  let mut dS = 0.0
  let mut n = 0
  for row in rhumb_indep_inverse {
    let rh = indep_rhumb(row[0])
    let r = rh.inverse(row[1], row[2], row[3], row[4])
    n = n + 1
    ds = indep_max(ds, (r.s12 - row[6]).abs())
    dS = indep_max(dS, (r.area12 - row[7]).abs() / rh.ellipsoid_area())
    if row[5].is_nan() {
      assert_true(
        r.azi12.is_nan(),
        msg="azimuth defined at \{row[1]} \{row[2]} \{row[3]} \{row[4]}",
      )
    } else {
      // Azimuth error as a sideways offset at the far end of the line.
      let off = indep_ang(r.azi12, row[5]) * 3.141592653589793 / 180.0 * row[6]
      da = indep_max(da, off)
    }
  }
  println(
    "rhumb inverse vs mpmath: n=\{n} max |ds12| m=\{ds} max azimuth offset m=\{da} max |dS12|/area=\{dS}",
  )
  assert_eq(n, 303)
  assert_true(ds < 1.2e-8) // README: distance < 12 nm
  assert_true(da < 7.0e-9) // README: azimuth < 7 nm at the far end
  assert_true(dS < 2.0e-16) // README: area < 2e-16 of the ellipsoid
}

///|
test "rhumb direct matches independent mpmath formulas" {
  let mut dp = 0.0
  let mut dS = 0.0
  let mut near_pos = 0.0
  let mut near_S = 0.0
  let mut n = 0
  let mut n_pole = 0
  let mut n_near = 0
  for row in rhumb_indep_direct {
    let rh = indep_rhumb(row[0])
    let r = rh.direct(row[2], row[3], row[4], row[5], long_unroll=row[1] == 1.0)
    n = n + 1
    let dlat = (r.lat2 - row[6]).abs() * indep_m_per_deg
    if row[7].is_nan() {
      // The line passes a pole: only the latitude is defined.
      n_pole = n_pole + 1
      assert_true(r.lon2.is_nan(), msg="lon2 at \{row[2]} \{row[4]} \{row[5]}")
      assert_true(
        r.area12.is_nan(),
        msg="area12 at \{row[2]} \{row[4]} \{row[5]}",
      )
      dp = indep_max(dp, dlat)
    } else if @cmp.maximum(row[2].abs(), row[6].abs()) >= 89.9 {
      // The line runs within 0.1 degree (11 km) of a pole, where it may
      // circle the pole thousands of times (2 900 turns 1 km from it).
      // Errors are compared with the distance travelled and the area swept.
      n_near = n_near + 1
      let u = rh.direct(row[2], row[3], row[4], row[5], long_unroll=true)
      let coslat = @math.cos(row[6] * 3.141592653589793 / 180.0)
      let dlon12 = (u.lon2 - u.lon1 - row[9]).abs()
      let dpos = dlat + dlon12 * coslat * indep_m_per_deg
      // Allowed: 11 nm plus 4e-13 of the distance travelled.
      near_pos = indep_max(near_pos, dpos - 4.0e-13 * row[5].abs())
      let swept = @cmp.maximum(row[8].abs(), rh.ellipsoid_area())
      near_S = indep_max(near_S, (r.area12 - row[8]).abs() / swept)
    } else {
      let coslat = @math.cos(row[6] * 3.141592653589793 / 180.0)
      let dlon = if row[1] == 1.0 {
        (r.lon2 - row[7]).abs()
      } else {
        indep_ang(r.lon2, row[7])
      }
      dp = indep_max(dp, dlat + dlon * coslat * indep_m_per_deg)
      dS = indep_max(dS, (r.area12 - row[8]).abs() / rh.ellipsoid_area())
    }
  }
  println(
    "rhumb direct vs mpmath: n=\{n} past a pole=\{n_pole} near a pole=\{n_near} max position err m=\{dp} max |dS12|/area=\{dS} near-pole position err - 4e-13 s12 m=\{near_pos} near-pole |dS12|/max(|S12|, area)=\{near_S}",
  )
  assert_eq(n, 179)
  assert_eq(n_pole, 37)
  assert_true(dp < 1.1e-8) // README: position < 11 nm
  assert_true(dS < 2.0e-16) // README: area < 2e-16 of the ellipsoid
  assert_eq(n_near, 19)
  assert_true(near_pos < 1.1e-8) // README: near a pole, 11 nm + 4e-13 s12
  assert_true(near_S < 4.0e-13) // README: area 4e-13 of the area swept or ellipsoid
}
