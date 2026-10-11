// Port of GeographicLib's GeodesicLine (Charles Karney, MIT/X11 License).

///|
/// A geodesic line starting at (lat1, lon1) with azimuth azi1; positions
/// along it can be computed cheaply once it has been constructed.
pub struct GeodesicLine {
  a : Double
  f : Double
  b : Double
  c2 : Double
  f1 : Double
  lat1 : Double
  lon1 : Double
  azi1 : Double
  salp1 : Double
  calp1 : Double
  dn1 : Double
  salp0 : Double
  calp0 : Double
  ssig1 : Double
  csig1 : Double
  somg1 : Double
  comg1 : Double
  k2 : Double
  a1m1 : Double
  c1a : FixedArray[Double]
  b11 : Double
  stau1 : Double
  ctau1 : Double
  c1pa : FixedArray[Double]
  a2m1 : Double
  c2a : FixedArray[Double]
  b21 : Double
  c3a : FixedArray[Double]
  a3c : Double
  b31 : Double
  c4a : FixedArray[Double]
  a4 : Double
  b41 : Double
  mut s13 : Double
  mut a13 : Double
}

///|
fn GeodesicLine::make(
  geod : Geodesic,
  lat1 : Double,
  lon1 : Double,
  azi1_in : Double,
  salp1_in : Double,
  calp1_in : Double,
) -> GeodesicLine {
  let lat1 = @geomath.lat_fix(lat1)
  let (azi1, salp1, calp1) = if salp1_in.is_nan() || calp1_in.is_nan() {
    let (s, c) = @geomath.sincosd(@geomath.ang_round(azi1_in))
    (@geomath.ang_normalize(azi1_in), s, c)
  } else {
    (azi1_in, salp1_in, calp1_in)
  }
  let (sb, cb) = @geomath.sincosd(@geomath.ang_round(lat1))
  let (sbet1, cb1) = @geomath.norm2(sb * geod.f1, cb)
  let cbet1 = @geomath.fmax(tiny, cb1)
  let dn1 = (1.0 + geod.ep2 * @geomath.sq(sbet1)).sqrt()
  let salp0 = salp1 * cbet1
  let calp0 = @geomath.hypot(calp1, salp1 * sbet1)
  let somg1 = salp0 * sbet1
  let c0 = if sbet1 != 0.0 || calp1 != 0.0 { cbet1 * calp1 } else { 1.0 }
  let comg1 = c0
  let (ssig1, csig1) = @geomath.norm2(sbet1, c0)
  let k2 = @geomath.sq(calp0) * geod.ep2
  let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2)
  let a1m1 = a1m1f(eps)
  let c1a = FixedArray::make(order + 1, 0.0)
  c1f(eps, c1a)
  let b11 = sin_cos_series(true, ssig1, csig1, c1a)
  let s = @math.sin(b11)
  let c = @math.cos(b11)
  let stau1 = ssig1 * c + csig1 * s
  let ctau1 = csig1 * c - ssig1 * s
  let c1pa = FixedArray::make(order + 1, 0.0)
  c1pf(eps, c1pa)
  let a2m1 = a2m1f(eps)
  let c2a = FixedArray::make(order + 1, 0.0)
  c2f(eps, c2a)
  let b21 = sin_cos_series(true, ssig1, csig1, c2a)
  let c3a = FixedArray::make(order, 0.0)
  geod.c3f(eps, c3a)
  let a3c = -geod.f * salp0 * geod.a3f(eps)
  let b31 = sin_cos_series(true, ssig1, csig1, c3a)
  let c4a = FixedArray::make(order, 0.0)
  geod.c4f(eps, c4a)
  let a4 = @geomath.sq(geod.a) * calp0 * salp0 * geod.e2
  let b41 = sin_cos_series(false, ssig1, csig1, c4a)
  {
    a: geod.a,
    f: geod.f,
    b: geod.b,
    c2: geod.c2,
    f1: geod.f1,
    lat1,
    lon1,
    azi1,
    salp1,
    calp1,
    dn1,
    salp0,
    calp0,
    ssig1,
    csig1,
    somg1,
    comg1,
    k2,
    a1m1,
    c1a,
    b11,
    stau1,
    ctau1,
    c1pa,
    a2m1,
    c2a,
    b21,
    c3a,
    a3c,
    b31,
    c4a,
    a4,
    b41,
    s13: @double.not_a_number,
    a13: @double.not_a_number,
  }
}

///|
/// Result of a direct problem / a position on a geodesic line.
pub(all) struct DirectResult {
  /// Latitude of point 1 (degrees).
  lat1 : Double
  /// Longitude of point 1 (degrees).
  lon1 : Double
  /// Azimuth at point 1 (degrees clockwise from north).
  azi1 : Double
  /// Latitude of point 2 (degrees). NaN when point 1 was outside [-90, 90].
  lat2 : Double
  /// Longitude of point 2 (degrees).
  lon2 : Double
  /// Forward azimuth at point 2 (degrees clockwise from north).
  azi2 : Double
  /// Distance from point 1 to point 2 (m).
  s12 : Double
  /// Arc length on the auxiliary sphere (degrees).
  a12 : Double
  /// Reduced length (m).
  m12 : Double
  /// Geodesic scale M12 (dimensionless).
  scale12 : Double
  /// Geodesic scale M21 (dimensionless).
  scale21 : Double
  /// Area between the geodesic and the equator (m²).
  area12 : Double
} derive(Debug)

///|
pub extend DirectResult with @debug.Debug::{to_repr}

///|
fn GeodesicLine::gen_position(
  self : GeodesicLine,
  arcmode : Bool,
  s12_a12 : Double,
  long_unroll : Bool,
) -> DirectResult {
  let mut b12 = 0.0
  let mut ab1 = 0.0
  let mut sig12 = 0.0
  let mut ssig12 = 0.0
  let mut csig12 = 0.0
  if arcmode {
    sig12 = @geomath.radians(s12_a12)
    let (s, c) = @geomath.sincosd(s12_a12)
    ssig12 = s
    csig12 = c
  } else {
    let mut tau12 = s12_a12 / (self.b * (1.0 + self.a1m1))
    if !@geomath.is_finite(tau12) {
      tau12 = @double.not_a_number
    }
    let s = @math.sin(tau12)
    let c = @math.cos(tau12)
    b12 = -sin_cos_series(
      true,
      self.stau1 * c + self.ctau1 * s,
      self.ctau1 * c - self.stau1 * s,
      self.c1pa,
    )
    sig12 = tau12 - (b12 - self.b11)
    ssig12 = @math.sin(sig12)
    csig12 = @math.cos(sig12)
    if self.f.abs() > 0.01 {
      let ssig2 = self.ssig1 * csig12 + self.csig1 * ssig12
      let csig2 = self.csig1 * csig12 - self.ssig1 * ssig12
      b12 = sin_cos_series(true, ssig2, csig2, self.c1a)
      let serr = (1.0 + self.a1m1) * (sig12 + (b12 - self.b11)) -
        s12_a12 / self.b
      sig12 = sig12 - serr / (1.0 + self.k2 * @geomath.sq(ssig2)).sqrt()
      ssig12 = @math.sin(sig12)
      csig12 = @math.cos(sig12)
    }
  }
  let ssig2 = self.ssig1 * csig12 + self.csig1 * ssig12
  let mut csig2 = self.csig1 * csig12 - self.ssig1 * ssig12
  let dn2 = (1.0 + self.k2 * @geomath.sq(ssig2)).sqrt()
  if arcmode || self.f.abs() > 0.01 {
    b12 = sin_cos_series(true, ssig2, csig2, self.c1a)
  }
  ab1 = (1.0 + self.a1m1) * (b12 - self.b11)
  let sbet2 = self.calp0 * ssig2
  let mut cbet2 = @geomath.hypot(self.salp0, self.calp0 * csig2)
  if cbet2 == 0.0 {
    cbet2 = tiny
    csig2 = tiny
  }
  let salp2 = self.salp0
  let calp2 = self.calp0 * csig2
  let s12 = if arcmode {
    self.b * ((1.0 + self.a1m1) * sig12 + ab1)
  } else {
    s12_a12
  }
  // longitude
  let somg2 = self.salp0 * ssig2
  let comg2 = csig2
  let e = @geomath.copysign(1.0, self.salp0)
  let omg12 = if long_unroll {
    e *
    (
      sig12 -
      (@math.atan2(ssig2, csig2) - @math.atan2(self.ssig1, self.csig1)) +
      (@math.atan2(e * somg2, comg2) - @math.atan2(e * self.somg1, self.comg1))
    )
  } else {
    @math.atan2(
      somg2 * self.comg1 - comg2 * self.somg1,
      comg2 * self.comg1 + somg2 * self.somg1,
    )
  }
  let lam12 = omg12 +
    self.a3c *
    (sig12 + (sin_cos_series(true, ssig2, csig2, self.c3a) - self.b31))
  let lon12 = @geomath.degrees(lam12)
  let lon2 = if long_unroll {
    self.lon1 + lon12
  } else {
    @geomath.ang_normalize(
      @geomath.ang_normalize(self.lon1) + @geomath.ang_normalize(lon12),
    )
  }
  let lat2 = @geomath.atan2d(sbet2, self.f1 * cbet2)
  let azi2 = @geomath.atan2d(salp2, calp2)
  // reduced length and geodesic scale
  let b22 = sin_cos_series(true, ssig2, csig2, self.c2a)
  let ab2 = (1.0 + self.a2m1) * (b22 - self.b21)
  let j12 = (self.a1m1 - self.a2m1) * sig12 + (ab1 - ab2)
  let m12 = self.b *
    (
      dn2 * (self.csig1 * ssig2) -
      self.dn1 * (self.ssig1 * csig2) -
      self.csig1 * csig2 * j12
    )
  let t = self.k2 *
    (ssig2 - self.ssig1) *
    (ssig2 + self.ssig1) /
    (self.dn1 + dn2)
  let mm12 = csig12 + (t * ssig2 - csig2 * j12) * self.ssig1 / self.dn1
  let mm21 = csig12 - (t * self.ssig1 - self.csig1 * j12) * ssig2 / dn2
  // area
  let b42 = sin_cos_series(false, ssig2, csig2, self.c4a)
  let (salp12, calp12) = if self.calp0 == 0.0 || self.salp0 == 0.0 {
    (
      salp2 * self.calp1 - calp2 * self.salp1,
      calp2 * self.calp1 + salp2 * self.salp1,
    )
  } else {
    (
      self.calp0 *
      self.salp0 *
      (if csig12 <= 0.0 {
        self.csig1 * (1.0 - csig12) + ssig12 * self.ssig1
      } else {
        ssig12 * (self.csig1 * ssig12 / (1.0 + csig12) + self.ssig1)
      }),
      @geomath.sq(self.salp0) + @geomath.sq(self.calp0) * self.csig1 * csig2,
    )
  }
  let area12 = self.c2 * @math.atan2(salp12, calp12) +
    self.a4 * (b42 - self.b41)
  let a12 = if arcmode { s12_a12 } else { @geomath.degrees(sig12) }
  {
    lat1: self.lat1,
    lon1: if long_unroll {
      self.lon1
    } else {
      @geomath.ang_normalize(self.lon1)
    },
    azi1: self.azi1,
    lat2,
    lon2,
    azi2,
    s12,
    a12,
    m12,
    scale12: mm12,
    scale21: mm21,
    area12,
  }
}

///|
/// Position at distance `s12` metres from point 1.
/// See `Geodesic::inverse` for units, poles, and `long_unroll`.
pub fn GeodesicLine::position(
  self : GeodesicLine,
  s12 : Double,
  long_unroll? : Bool = false,
) -> DirectResult {
  self.gen_position(false, s12, long_unroll)
}

///|
/// Position at arc length `a12` @geomath.degrees from point 1.
pub fn GeodesicLine::arc_position(
  self : GeodesicLine,
  a12 : Double,
  long_unroll? : Bool = false,
) -> DirectResult {
  self.gen_position(true, a12, long_unroll)
}

///|
/// Distance to the reference point 3 (set by `inverse_line` / `direct_line`).
pub fn GeodesicLine::distance(self : GeodesicLine) -> Double {
  self.s13
}

///|
/// Arc length to the reference point 3 in @geomath.degrees.
pub fn GeodesicLine::arc(self : GeodesicLine) -> Double {
  self.a13
}

///|
/// Set reference point 3 to distance `s13` metres from point 1.
/// `distance` and `densify` use this point.
pub fn GeodesicLine::set_distance(self : GeodesicLine, s13 : Double) -> Unit {
  self.s13 = s13
  self.a13 = self.gen_position(false, s13, false).a12
}

///|
/// Set reference point 3 to arc length `a13` degrees from point 1.
pub fn GeodesicLine::set_arc(self : GeodesicLine, a13 : Double) -> Unit {
  self.a13 = a13
  self.s13 = self.gen_position(true, a13, false).s12
}

///|
/// Points spaced at most `max_step` metres apart along the line from
/// point 1 to point 3 (inclusive), for example to draw a route.
/// Each point is `(lat, lon)` in degrees.
/// Point 3 must already be set (`inverse_line`, `direct_line`,
/// `set_distance`, or `set_arc`). Otherwise the distance is NaN.
pub fn GeodesicLine::densify(
  self : GeodesicLine,
  max_step : Double,
) -> Array[(Double, Double)] {
  let total = self.s13
  let n = @geomath.fmax(1.0, (total / max_step).ceil()).to_int()
  let pts = []
  for i in 0..<=n {
    let p = self.position(total * i.to_double() / n.to_double())
    pts.push((p.lat2, p.lon2))
  }
  pts
}

///|
/// A geodesic line starting at `(lat1, lon1)` with azimuth `azi1` (degrees).
/// `distance()` is NaN until `set_distance`, `set_arc`, or `direct_line`.
pub fn Geodesic::line(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
) -> GeodesicLine {
  GeodesicLine::make(
    self, lat1, lon1, azi1, @double.not_a_number, @double.not_a_number,
  )
}

///|
/// The geodesic line between two points; `distance()` returns their separation.
pub fn Geodesic::inverse_line(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> GeodesicLine {
  let r = self.gen_inverse(lat1, lon1, lat2, lon2)
  let azi1 = @geomath.atan2d(r.salp1, r.calp1)
  let line = GeodesicLine::make(self, lat1, lon1, azi1, r.salp1, r.calp1)
  line.set_arc(r.a12)
  line
}

///|
/// A geodesic line with point 3 at distance `s12` from point 1.
pub fn Geodesic::direct_line(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  s12 : Double,
) -> GeodesicLine {
  let line = self.line(lat1, lon1, azi1)
  line.set_distance(s12)
  line
}

///|
/// Solve the direct geodesic problem: travel `s12` metres from (lat1, lon1)
/// with initial azimuth `azi1` (degrees, clockwise from north).
/// See `inverse` for poles, NaN latitudes, and `long_unroll`.
pub fn Geodesic::direct(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  s12 : Double,
  long_unroll? : Bool = false,
) -> DirectResult {
  self.line(lat1, lon1, azi1).gen_position(false, s12, long_unroll)
}

///|
/// Direct problem with the distance given as an arc length in @geomath.degrees.
pub fn Geodesic::arc_direct(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  a12 : Double,
  long_unroll? : Bool = false,
) -> DirectResult {
  self.line(lat1, lon1, azi1).gen_position(true, a12, long_unroll)
}
