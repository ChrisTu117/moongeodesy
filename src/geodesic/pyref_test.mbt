///|
/// |a - b| scaled: absolute below 1, relative above.
fn rel(a : Double, b : Double) -> Double {
  if a.is_nan() && b.is_nan() {
    return 0.0
  }
  let d = (a - b).abs()
  let s = b.abs()
  if s > 1.0 {
    d / s
  } else {
    d
  }
}

///|
fn angd(a : Double, b : Double) -> Double {
  let mut d = (a - b).abs()
  if d > 180.0 {
    d = 360.0 - d
  }
  d
}

///|
test "inverse agrees with reference Python geographiclib" {
  let g = @geodesic.Geodesic::wgs84()
  let mut worst = 0.0
  let mut worst_ang = 0.0
  for v in pyref_inverse {
    let r = g.inverse(v[0], v[1], v[2], v[3])
    for
      p in [
        (r.s12, v[4]),
        (r.a12, v[7]),
        (r.m12, v[8]),
        (r.scale12, v[9]),
        (r.scale21, v[10]),
      ] {
      let e = if v[7] > 179.5 && p.1 != v[4] { 0.0 } else { rel(p.0, p.1) }
      if e > worst {
        worst = e
      }
    }
    // Azimuths: compare as displacement at the far end (m); short lines are
    // ill-conditioned in angle. Skip nearly antipodal (non-unique) cases.
    let ea = if v[7] > 179.5 {
      0.0
    } else {
      @cmp.maximum(angd(r.azi1, v[5]), angd(r.azi2, v[6])) *
      @cmp.maximum(v[4], 1.0) *
      0.0174533
    }
    if ea > worst_ang {
      worst_ang = ea
    }
    // area: relative to the ellipsoid's total area
    let es = if v[7] > 179.5 {
      0.0
    } else {
      (r.area12 - v[11]).abs() / g.ellipsoid_area()
    }
    if es > worst {
      worst = es
    }
  }
  println(
    "inverse vs python: worst rel \{worst}, worst azimuth displacement \{worst_ang} m",
  )
  assert_true(worst < 1.0e-11)
  assert_true(worst_ang < 1.0e-7)
}

///|
test "direct agrees with reference Python geographiclib" {
  let g = @geodesic.Geodesic::wgs84()
  let mut worst = 0.0
  let mut worst_ang = 0.0
  for v in pyref_direct {
    let r = g.direct(v[0], v[1], v[2], v[3], long_unroll=true)
    for
      p in [(r.a12, v[7]), (r.m12, v[8]), (r.scale12, v[9]), (r.scale21, v[10])] {
      let e = rel(p.0, p.1)
      if e > worst {
        worst = e
      }
    }
    for p in [(r.lat2, v[4]), (r.lon2, v[5]), (r.azi2, v[6])] {
      let e = angd(p.0, p.1)
      if e > worst_ang {
        worst_ang = e
      }
    }
    let es = (r.area12 - v[11]).abs() / g.ellipsoid_area()
    if es > worst {
      worst = es
    }
  }
  println("direct vs python: worst rel \{worst}, worst angle \{worst_ang} deg")
  assert_true(worst < 1.0e-12)
  assert_true(worst_ang < 1.0e-12)
}
