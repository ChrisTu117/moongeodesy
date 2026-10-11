// Transverse Mercator projection using Krüger's series to order n^6, after
// C. F. F. Karney, "Transverse Mercator with an accuracy of a few
// nanometers", J. Geodesy 85, 475-485 (2011),
// https://doi.org/10.1007/s00190-011-0445-3 — following the structure of
// GeographicLib::TransverseMercator (MIT/X11 License). Accurate to a few
// nanometres within ~3900 km of the central meridian.

///|
/// e * atanh(e * x) (or the oblate/prolate analogue) for signed `es`.
fn eatanhe(x : Double, es : Double) -> Double {
  if es > 0.0 {
    es * @math.atanh(es * x)
  } else {
    -es * @math.atan(es * x)
  }
}

///|
/// tan(chi) where chi is the conformal latitude, given tau = tan(phi).
fn taupf(tau : Double, es : Double) -> Double {
  let tau1 = @geomath.hypot(1.0, tau)
  let sig = @math.sinh(eatanhe(tau / tau1, es))
  @geomath.hypot(1.0, sig) * tau - sig * tau1
}

///|
/// Inverse of `taupf` by Newton's method.
fn tauf(taup : Double, es : Double) -> Double {
  let tol = @geomath.dbl_epsilon.sqrt() / 10.0
  let taumax = 2.0 / @geomath.dbl_epsilon.sqrt()
  let e2m = 1.0 - es * es
  let mut tau = if taup.abs() > 70.0 {
    taup * @math.exp(eatanhe(1.0, es))
  } else {
    taup / e2m
  }
  let stol = tol * @geomath.fmax(1.0, taup.abs())
  if !(tau.abs() < taumax) {
    return tau
  }
  for _ in 0..<5 {
    let taupa = taupf(tau, es)
    let dtau = (taup - taupa) *
      (1.0 + e2m * tau * tau) /
      (e2m * @geomath.hypot(1.0, tau) * @geomath.hypot(1.0, taupa))
    tau = tau + dtau
    if !(dtau.abs() >= stol) {
      break
    }
  }
  tau
}

///|
/// A transverse Mercator projection on an ellipsoid.
pub struct TransverseMercator {
  a : Double
  f : Double
  k0 : Double
  e2 : Double
  es : Double
  e2m : Double
  c : Double
  n : Double
  a1 : Double
  b1 : Double
  alp : FixedArray[Double]
  bet : FixedArray[Double]
}

///|
/// Result of a forward projection.
pub(all) struct Projected {
  /// Easting from the origin (m). No false easting is applied.
  x : Double
  /// Northing from the origin (m). No false northing is applied.
  y : Double
  /// Meridian convergence (degrees).
  gamma : Double
  /// Point scale factor (dimensionless).
  k : Double
} derive(Debug)

///|
pub extend Projected with @debug.Debug::{to_repr}

///|
/// Result of a reverse projection.
pub(all) struct Geographic {
  /// Latitude (degrees).
  lat : Double
  /// Longitude (degrees).
  lon : Double
  /// Meridian convergence (degrees).
  gamma : Double
  /// Point scale factor (dimensionless).
  k : Double
} derive(Debug)

///|
pub extend Geographic with @debug.Debug::{to_repr}

///|
/// Create a transverse Mercator projection for an ellipsoid with equatorial
/// radius `a`, flattening `f` and central scale factor `k0`.
pub fn TransverseMercator::new(
  a : Double,
  f : Double,
  k0 : Double,
) -> TransverseMercator {
  let e2 = f * (2.0 - f)
  let es = (if f < 0.0 { -1.0 } else { 1.0 }) * e2.abs().sqrt()
  let e2m = 1.0 - e2
  let c = e2m.sqrt() * @math.exp(eatanhe(1.0, es))
  let n = f / (2.0 - f)
  let n2 = n * n
  let b1 = (1.0 + n2 * (1.0 / 4.0 + n2 * (1.0 / 64.0 + n2 / 256.0))) / (1.0 + n)
  let a1 = b1 * a
  // Krüger series coefficients (Karney 2011, eqs. 35-36), order n^6.
  let n3 = n2 * n
  let n4 = n3 * n
  let n5 = n4 * n
  let n6 = n5 * n
  let alp : FixedArray[Double] = [
    0.0,
    n / 2.0 -
    2.0 * n2 / 3.0 +
    5.0 * n3 / 16.0 +
    41.0 * n4 / 180.0 -
    127.0 * n5 / 288.0 +
    7891.0 * n6 / 37800.0,
    13.0 * n2 / 48.0 -
    3.0 * n3 / 5.0 +
    557.0 * n4 / 1440.0 +
    281.0 * n5 / 630.0 -
    1983433.0 * n6 / 1935360.0,
    61.0 * n3 / 240.0 -
    103.0 * n4 / 140.0 +
    15061.0 * n5 / 26880.0 +
    167603.0 * n6 / 181440.0,
    49561.0 * n4 / 161280.0 - 179.0 * n5 / 168.0 + 6601661.0 * n6 / 7257600.0,
    34729.0 * n5 / 80640.0 - 3418889.0 * n6 / 1995840.0,
    212378941.0 * n6 / 319334400.0,
  ]
  let bet : FixedArray[Double] = [
    0.0,
    n / 2.0 -
    2.0 * n2 / 3.0 +
    37.0 * n3 / 96.0 -
    n4 / 360.0 -
    81.0 * n5 / 512.0 +
    96199.0 * n6 / 604800.0,
    n2 / 48.0 +
    n3 / 15.0 -
    437.0 * n4 / 1440.0 +
    46.0 * n5 / 105.0 -
    1118711.0 * n6 / 3870720.0,
    17.0 * n3 / 480.0 -
    37.0 * n4 / 840.0 -
    209.0 * n5 / 4480.0 +
    5569.0 * n6 / 90720.0,
    4397.0 * n4 / 161280.0 - 11.0 * n5 / 504.0 - 830251.0 * n6 / 7257600.0,
    4583.0 * n5 / 161280.0 - 108847.0 * n6 / 3991680.0,
    20648693.0 * n6 / 638668800.0,
  ]
  { a, f, k0, e2, es, e2m, c, n, a1, b1, alp, bet, }
}

///|
/// The UTM projection's transverse Mercator (WGS-84, k0 = 0.9996).
pub fn TransverseMercator::utm() -> TransverseMercator {
  utm_tm
}

///|
let utm_tm : TransverseMercator = TransverseMercator::new(
  6378137.0,
  1.0 / 298.257223563,
  0.9996,
)

///|
/// Project `(lat, lon)` with central meridian `lon0`. All angles are degrees.
/// `x` is easting and `y` is northing from the central meridian and the
/// equator (m). No false offsets are applied.
/// The 6th-order series stays within a few nanometres for points within
/// about 3900 km of the central meridian. Latitude ±90 is valid.
pub fn TransverseMercator::forward(
  self : TransverseMercator,
  lon0 : Double,
  lat_in : Double,
  lon_in : Double,
) -> Projected {
  let lat0 = @geomath.lat_fix(lat_in)
  let (lon_d, _) = @geomath.ang_diff(lon0, lon_in)
  let mut latsign = if lat0 < 0.0 ||
    (lat0 == 0.0 && @geomath.copysign(1.0, lat0) < 0.0) {
    -1.0
  } else {
    1.0
  }
  let lonsign = if @geomath.copysign(1.0, lon_d) < 0.0 { -1.0 } else { 1.0 }
  let lat = lat0 * latsign
  let mut lon = lon_d * lonsign
  let backside = lon > 90.0
  if backside {
    if lat == 0.0 {
      latsign = -1.0
    }
    lon = 180.0 - lon
  }
  let (sphi, cphi) = @geomath.sincosd(lat)
  let (slam, clam) = @geomath.sincosd(lon)
  let mut xip = 0.0
  let mut etap = 0.0
  let mut gamma = 0.0
  let mut k = 0.0
  if lat != 90.0 {
    let tau = sphi / cphi
    let taup = taupf(tau, self.es)
    xip = @math.atan2(taup, clam)
    etap = @math.asinh(slam / @geomath.hypot(taup, clam))
    gamma = @geomath.atan2d(slam * taup, clam * @geomath.hypot(1.0, taup))
    k = (self.e2m + self.e2 * cphi * cphi).sqrt() *
      @geomath.hypot(1.0, tau) /
      @geomath.hypot(taup, clam)
  } else {
    xip = @geomath.pi / 2.0
    etap = 0.0
    gamma = lon
    k = self.c
  }
  // xi + i eta = zeta' + sum alp_j sin(2 j zeta')
  let mut xi = xip
  let mut eta = etap
  let mut p = 1.0
  let mut q = 0.0
  for j in 1..<=6 {
    let tj = 2.0 * j.to_double()
    let s = @math.sin(tj * xip)
    let c = @math.cos(tj * xip)
    let sh = @math.sinh(tj * etap)
    let ch = @math.cosh(tj * etap)
    xi = xi + self.alp[j] * s * ch
    eta = eta + self.alp[j] * c * sh
    p = p + tj * self.alp[j] * c * ch
    q = q + tj * self.alp[j] * s * sh
  }
  gamma = gamma + @geomath.atan2d(q, p)
  k = k * self.b1 * @geomath.hypot(p, q)
  let y = self.a1 *
    self.k0 *
    (if backside { @geomath.pi - xi } else { xi }) *
    latsign
  let x = self.a1 * self.k0 * eta * lonsign
  if backside {
    gamma = 180.0 - gamma
  }
  gamma = gamma * latsign * lonsign
  gamma = @geomath.ang_normalize(gamma)
  { x, y, gamma, k: k * self.k0, }
}

///|
/// Inverse projection: (x, y) relative to central meridian `lon0`.
pub fn TransverseMercator::reverse(
  self : TransverseMercator,
  lon0 : Double,
  x : Double,
  y : Double,
) -> Geographic {
  let mut xi = y / (self.a1 * self.k0)
  let mut eta = x / (self.a1 * self.k0)
  let xisign = if @geomath.copysign(1.0, xi) < 0.0 { -1.0 } else { 1.0 }
  let etasign = if @geomath.copysign(1.0, eta) < 0.0 { -1.0 } else { 1.0 }
  xi = xi * xisign
  eta = eta * etasign
  let backside = xi > @geomath.pi / 2.0
  if backside {
    xi = @geomath.pi - xi
  }
  let mut xip = xi
  let mut etap = eta
  let mut p = 1.0
  let mut q = 0.0
  for j in 1..<=6 {
    let tj = 2.0 * j.to_double()
    let s = @math.sin(tj * xi)
    let c = @math.cos(tj * xi)
    let sh = @math.sinh(tj * eta)
    let ch = @math.cosh(tj * eta)
    xip = xip - self.bet[j] * s * ch
    etap = etap - self.bet[j] * c * sh
    p = p - tj * self.bet[j] * c * ch
    q = q + tj * self.bet[j] * s * sh
  }
  let mut gamma = @geomath.atan2d(q, p)
  let mut k = self.b1 / @geomath.hypot(p, q)
  let s = @math.sinh(etap)
  let c = @geomath.fmax(0.0, @math.cos(xip))
  let r = @geomath.hypot(s, c)
  let mut lat = 0.0
  let mut lon = 0.0
  if r != 0.0 {
    lon = @geomath.atan2d(s, c)
    let sxip = @math.sin(xip)
    let tau = tauf(sxip / r, self.es)
    gamma = gamma + @geomath.atan2d(sxip * @math.tanh(etap), c)
    lat = @geomath.atan2d(tau, 1.0)
    k = k *
      (self.e2m + self.e2 / (1.0 + tau * tau)).sqrt() *
      @geomath.hypot(1.0, tau) *
      r
  } else {
    lat = 90.0
    lon = 0.0
    k = k * self.c
  }
  lat = lat * xisign
  if backside {
    lon = 180.0 - lon
  }
  lon = lon * etasign
  lon = @geomath.ang_normalize(lon + lon0)
  if backside {
    gamma = 180.0 - gamma
  }
  gamma = @geomath.ang_normalize(gamma * xisign * etasign)
  { lat, lon, gamma, k: k * self.k0, }
}
