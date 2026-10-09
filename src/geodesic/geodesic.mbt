// MoonBit port of GeographicLib's Geodesic class (Karney 2013,
// "Algorithms for geodesics", J. Geodesy 87, 43-55,
// https://doi.org/10.1007/s00190-012-0578-z).
// Original: Charles Karney, MIT/X11 License, https://geographiclib.sourceforge.io/
// This is a close transcription of the Python implementation (geographiclib 2.x).

///|
let order : Int = 6

///|
let n_c3x : Int = 15 // order * (order - 1) / 2

///|
let n_c4x : Int = 21 // order * (order + 1) / 2

///|
let maxit1 : Int = 20

///|
let maxit2 : Int = 83 // maxit1 + 53 + 10

///|
let tiny : Double = dbl_min.sqrt()

///|
let tol0 : Double = dbl_epsilon

///|
let tol1 : Double = 200.0 * dbl_epsilon

///|
let tol2 : Double = dbl_epsilon.sqrt()

///|
let tolb : Double = dbl_epsilon

///|
let xthresh : Double = 1000.0 * dbl_epsilon.sqrt()

// Output masks (subset of GeodesicCapability used internally).

///|
let m_distance : Int = 1 << 10

///|
let m_reducedlength : Int = 1 << 12

///|
let m_geodesicscale : Int = 1 << 13

///|
let m_area : Int = 1 << 14

///|
let m_all : Int = m_distance | m_reducedlength | m_geodesicscale | m_area

///|
/// Error raised for an invalid ellipsoid.
pub suberror GeodesicError {
  InvalidEllipsoid(String)
} derive(Debug)

///|
pub extend GeodesicError with @debug.Debug::{to_repr}

///|
/// An ellipsoid of revolution on which geodesic problems are solved.
pub struct Geodesic {
  /// Equatorial radius in metres.
  a : Double
  /// Flattening.
  f : Double
  f1 : Double
  e2 : Double
  ep2 : Double
  n : Double
  b : Double
  c2 : Double
  etol2 : Double
  a3x : FixedArray[Double]
  c3x : FixedArray[Double]
  c4x : FixedArray[Double]
}

///|
/// Construct an ellipsoid with equatorial radius `a` (metres) and
/// flattening `f`. Raises if `a` or the polar semi-axis is not positive.
pub fn Geodesic::new(a : Double, f : Double) -> Geodesic raise GeodesicError {
  let f1 = 1.0 - f
  let e2 = f * (2.0 - f)
  let ep2 = e2 / sq(f1)
  let n = f / (2.0 - f)
  let b = a * f1
  let c2 = (
      sq(a) +
      sq(b) *
      (if e2 == 0.0 {
        1.0
      } else {
        (if e2 > 0.0 {
          @math.atanh(e2.sqrt())
        } else {
          @math.atan((-e2).sqrt())
        }) /
        e2.abs().sqrt()
      })
    ) /
    2.0
  let etol2 = 0.1 *
    tol2 /
    (fmax(0.001, f.abs()) * fmin(1.0, 1.0 - f / 2.0) / 2.0).sqrt()
  if !(is_finite(a) && a > 0.0) {
    raise InvalidEllipsoid("Equatorial radius is not positive")
  }
  if !(is_finite(b) && b > 0.0) {
    raise InvalidEllipsoid("Polar semi-axis is not positive")
  }
  let g : Geodesic = {
    a,
    f,
    f1,
    e2,
    ep2,
    n,
    b,
    c2,
    etol2,
    a3x: FixedArray::make(order, 0.0),
    c3x: FixedArray::make(n_c3x, 0.0),
    c4x: FixedArray::make(n_c4x, 0.0),
  }
  g.a3coeff()
  g.c3coeff()
  g.c4coeff()
  g
}

///|
/// WGS-84 semi-major axis (m).
pub let wgs84_a : Double = 6378137.0

///|
/// WGS-84 flattening.
pub let wgs84_f : Double = 1.0 / 298.257223563

///|
let wgs84_geod : Geodesic = Geodesic::new(wgs84_a, wgs84_f) catch {
  _ => abort("unreachable")
}

///|
/// The WGS-84 ellipsoid.
pub fn Geodesic::wgs84() -> Geodesic {
  wgs84_geod
}

///|
pub fn Geodesic::equatorial_radius(self : Geodesic) -> Double {
  self.a
}

///|
pub fn Geodesic::flattening(self : Geodesic) -> Double {
  self.f
}

///|
/// Total area of the ellipsoid in m².
pub fn Geodesic::ellipsoid_area(self : Geodesic) -> Double {
  4.0 * pi * self.c2
}

///|
fn sin_cos_series(
  sinp : Bool,
  sinx : Double,
  cosx : Double,
  c : FixedArray[Double],
) -> Double {
  let mut k = c.length()
  let mut n = k - (if sinp { 1 } else { 0 })
  let ar = 2.0 * (cosx - sinx) * (cosx + sinx)
  let mut y1 = 0.0
  let mut y0 = 0.0
  if (n & 1) != 0 {
    k -= 1
    y0 = c[k]
  }
  n = n / 2
  while n > 0 {
    n -= 1
    k -= 1
    y1 = ar * y0 - y1 + c[k]
    k -= 1
    y0 = ar * y1 - y0 + c[k]
  }
  if sinp {
    2.0 * sinx * cosx * y0
  } else {
    cosx * (y0 - y1)
  }
}

///|
fn astroid(x : Double, y : Double) -> Double {
  let p = sq(x)
  let q = sq(y)
  let r = (p + q - 1.0) / 6.0
  if !(q == 0.0 && r <= 0.0) {
    let s = p * q / 4.0
    let r2 = sq(r)
    let r3 = r * r2
    let disc = s * (s + 2.0 * r3)
    let mut u = r
    if disc >= 0.0 {
      let mut t3 = s + r3
      t3 = t3 + (if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() })
      let t = @math.cbrt(t3)
      u = u + t + (if t != 0.0 { r2 / t } else { 0.0 })
    } else {
      let ang = @math.atan2((-disc).sqrt(), -(s + r3))
      u = u + 2.0 * r * @math.cos(ang / 3.0)
    }
    let v = (sq(u) + q).sqrt()
    let uv = if u < 0.0 { q / (v - u) } else { u + v }
    let w = (uv - q) / (2.0 * v)
    uv / ((uv + sq(w)).sqrt() + w)
  } else {
    0.0
  }
}

///|
let coeff_a1 : FixedArray[Double] = [1.0, 4.0, 64.0, 0.0, 256.0]

///|
fn a1m1f(eps : Double) -> Double {
  let m = order / 2
  let t = polyval(m, coeff_a1, 0, sq(eps)) / coeff_a1[m + 1]
  (t + eps) / (1.0 - eps)
}

///|
let coeff_c1 : FixedArray[Double] = [
  -1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0,
  512.0, -7.0, 1280.0, -7.0, 2048.0,
]

///|
fn series_eps2(
  eps : Double,
  coeff : FixedArray[Double],
  c : FixedArray[Double],
) -> Unit {
  let eps2 = sq(eps)
  let mut d = eps
  let mut o = 0
  for l in 1..<=order {
    let m = (order - l) / 2
    c[l] = d * polyval(m, coeff, o, eps2) / coeff[o + m + 1]
    o += m + 2
    d = d * eps
  }
}

///|
fn c1f(eps : Double, c : FixedArray[Double]) -> Unit {
  series_eps2(eps, coeff_c1, c)
}

///|
let coeff_c1p : FixedArray[Double] = [
  205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0,
  384.0, -7173.0, 2695.0, 7680.0, 3467.0, 7680.0, 38081.0, 61440.0,
]

///|
fn c1pf(eps : Double, c : FixedArray[Double]) -> Unit {
  series_eps2(eps, coeff_c1p, c)
}

///|
let coeff_a2 : FixedArray[Double] = [-11.0, -28.0, -192.0, 0.0, 256.0]

///|
fn a2m1f(eps : Double) -> Double {
  let m = order / 2
  let t = polyval(m, coeff_a2, 0, sq(eps)) / coeff_a2[m + 1]
  (t - eps) / (1.0 + eps)
}

///|
let coeff_c2 : FixedArray[Double] = [
  1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0,
  512.0, 63.0, 1280.0, 77.0, 2048.0,
]

///|
fn c2f(eps : Double, c : FixedArray[Double]) -> Unit {
  series_eps2(eps, coeff_c2, c)
}

///|
let coeff_a3 : FixedArray[Double] = [
  -3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0,
  -1.0, 2.0, 1.0, 1.0,
]

///|
fn Geodesic::a3coeff(self : Geodesic) -> Unit {
  let mut o = 0
  let mut k = 0
  for j = order - 1; j >= 0; j = j - 1 {
    let m = imin(order - j - 1, j)
    self.a3x[k] = polyval(m, coeff_a3, o, self.n) / coeff_a3[o + m + 1]
    k += 1
    o += m + 2
  }
}

///|
let coeff_c3 : FixedArray[Double] = [
  3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0,
  4.0, 5.0, 256.0, 1.0, 3.0, 128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0,
  7.0, 512.0, -10.0, 9.0, 384.0, 5.0, -9.0, 5.0, 192.0, 7.0, 512.0, -14.0, 7.0, 512.0,
  21.0, 2560.0,
]

///|
fn Geodesic::c3coeff(self : Geodesic) -> Unit {
  let mut o = 0
  let mut k = 0
  for l in 1..<order {
    for j = order - 1; j >= l; j = j - 1 {
      let m = imin(order - j - 1, j)
      self.c3x[k] = polyval(m, coeff_c3, o, self.n) / coeff_c3[o + m + 1]
      k += 1
      o += m + 2
    }
  }
}

///|
let coeff_c4 : FixedArray[Double] = [
  97.0, 15015.0, 1088.0, 156.0, 45045.0, -224.0, -4784.0, 1573.0, 45045.0, -10656.0,
  14144.0, -4576.0, -858.0, 45045.0, 64.0, 624.0, -4576.0, 6864.0, -3003.0, 15015.0,
  100.0, 208.0, 572.0, 3432.0, -12012.0, 30030.0, 45045.0, 1.0, 9009.0, -2944.0,
  468.0, 135135.0, 5792.0, 1040.0, -1287.0, 135135.0, 5952.0, -11648.0, 9152.0, -2574.0,
  135135.0, -64.0, -624.0, 4576.0, -6864.0, 3003.0, 135135.0, 8.0, 10725.0, 1856.0,
  -936.0, 225225.0, -8448.0, 4992.0, -1144.0, 225225.0, -1440.0, 4160.0, -4576.0,
  1716.0, 225225.0, -136.0, 63063.0, 1024.0, -208.0, 105105.0, 3584.0, -3328.0, 1144.0,
  315315.0, -128.0, 135135.0, -2560.0, 832.0, 405405.0, 128.0, 99099.0,
]

///|
fn Geodesic::c4coeff(self : Geodesic) -> Unit {
  let mut o = 0
  let mut k = 0
  for l in 0..<order {
    for j = order - 1; j >= l; j = j - 1 {
      let m = order - j - 1
      self.c4x[k] = polyval(m, coeff_c4, o, self.n) / coeff_c4[o + m + 1]
      k += 1
      o += m + 2
    }
  }
}

///|
fn Geodesic::a3f(self : Geodesic, eps : Double) -> Double {
  polyval(order - 1, self.a3x, 0, eps)
}

///|
fn Geodesic::c3f(self : Geodesic, eps : Double, c : FixedArray[Double]) -> Unit {
  let mut mult = 1.0
  let mut o = 0
  for l in 1..<order {
    let m = order - l - 1
    mult = mult * eps
    c[l] = mult * polyval(m, self.c3x, o, eps)
    o += m + 1
  }
}

///|
fn Geodesic::c4f(self : Geodesic, eps : Double, c : FixedArray[Double]) -> Unit {
  let mut mult = 1.0
  let mut o = 0
  for l in 0..<order {
    let m = order - l - 1
    c[l] = mult * polyval(m, self.c4x, o, eps)
    o += m + 1
    mult = mult * eps
  }
}

///|
/// Returns (s12b, m12b, m0, M12, M21).
fn Geodesic::lengths(
  self : Geodesic,
  eps : Double,
  sig12 : Double,
  ssig1 : Double,
  csig1 : Double,
  dn1 : Double,
  ssig2 : Double,
  csig2 : Double,
  dn2 : Double,
  cbet1 : Double,
  cbet2 : Double,
  outmask : Int,
  c1a : FixedArray[Double],
  c2a : FixedArray[Double],
) -> (Double, Double, Double, Double, Double) {
  let nan = @double.not_a_number
  let mut s12b = nan
  let mut m12b = nan
  let mut m0 = nan
  let mut m0x = nan
  let mut j12 = nan
  let mut mm12 = nan
  let mut mm21 = nan
  let mut a1 = nan
  let mut a2 = nan
  let red_or_scale = m_reducedlength | m_geodesicscale
  if (outmask & (m_distance | red_or_scale)) != 0 {
    a1 = a1m1f(eps)
    c1f(eps, c1a)
    if (outmask & red_or_scale) != 0 {
      a2 = a2m1f(eps)
      c2f(eps, c2a)
      m0x = a1 - a2
      a2 = 1.0 + a2
    }
    a1 = 1.0 + a1
  }
  if (outmask & m_distance) != 0 {
    let b1 = sin_cos_series(true, ssig2, csig2, c1a) -
      sin_cos_series(true, ssig1, csig1, c1a)
    s12b = a1 * (sig12 + b1)
    if (outmask & red_or_scale) != 0 {
      let b2 = sin_cos_series(true, ssig2, csig2, c2a) -
        sin_cos_series(true, ssig1, csig1, c2a)
      j12 = m0x * sig12 + (a1 * b1 - a2 * b2)
    }
  } else if (outmask & red_or_scale) != 0 {
    for l in 1..<order {
      c2a[l] = a1 * c1a[l] - a2 * c2a[l]
    }
    j12 = m0x * sig12 +
      (
        sin_cos_series(true, ssig2, csig2, c2a) -
        sin_cos_series(true, ssig1, csig1, c2a)
      )
  }
  if (outmask & m_reducedlength) != 0 {
    m0 = m0x
    m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * j12
  }
  if (outmask & m_geodesicscale) != 0 {
    let csig12 = csig1 * csig2 + ssig1 * ssig2
    let t = self.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
    mm12 = csig12 + (t * ssig2 - csig2 * j12) * ssig1 / dn1
    mm21 = csig12 - (t * ssig1 - csig1 * j12) * ssig2 / dn2
  }
  (s12b, m12b, m0, mm12, mm21)
}

///|
/// Returns (sig12, salp1, calp1, salp2, calp2, dnm).
fn Geodesic::inverse_start(
  self : Geodesic,
  sbet1 : Double,
  cbet1 : Double,
  dn1 : Double,
  sbet2 : Double,
  cbet2 : Double,
  dn2 : Double,
  lam12 : Double,
  slam12 : Double,
  clam12 : Double,
  c1a : FixedArray[Double],
  c2a : FixedArray[Double],
) -> (Double, Double, Double, Double, Double, Double) {
  let nan = @double.not_a_number
  let mut sig12 = -1.0
  let mut salp2 = nan
  let mut calp2 = nan
  let mut dnm = nan
  let sbet12 = sbet2 * cbet1 - cbet2 * sbet1
  let cbet12 = cbet2 * cbet1 + sbet2 * sbet1
  let mut sbet12a = sbet2 * cbet1
  sbet12a = sbet12a + cbet2 * sbet1
  let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5
  let mut somg12 = 0.0
  let mut comg12 = 0.0
  if shortline {
    let mut sbetm2 = sq(sbet1 + sbet2)
    sbetm2 = sbetm2 / (sbetm2 + sq(cbet1 + cbet2))
    dnm = (1.0 + self.ep2 * sbetm2).sqrt()
    let omg12 = lam12 / (self.f1 * dnm)
    somg12 = @math.sin(omg12)
    comg12 = @math.cos(omg12)
  } else {
    somg12 = slam12
    comg12 = clam12
  }
  let mut salp1 = cbet2 * somg12
  let mut calp1 = if comg12 >= 0.0 {
    sbet12 + cbet2 * sbet1 * sq(somg12) / (1.0 + comg12)
  } else {
    sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12)
  }
  let ssig12 = hypot(salp1, calp1)
  let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12
  if shortline && ssig12 < self.etol2 {
    salp2 = cbet1 * somg12
    calp2 = sbet12 -
      cbet1 *
      sbet2 *
      (if comg12 >= 0.0 { sq(somg12) / (1.0 + comg12) } else { 1.0 - comg12 })
    let (s, c) = norm2(salp2, calp2)
    salp2 = s
    calp2 = c
    sig12 = @math.atan2(ssig12, csig12)
  } else if self.n.abs() >= 0.1 ||
    csig12 >= 0.0 ||
    ssig12 >= 6.0 * self.n.abs() * pi * sq(cbet1) {
    ()
  } else {
    let lam12x = @math.atan2(-slam12, -clam12)
    let mut x = 0.0
    let mut y = 0.0
    let mut lamscale = 0.0
    let mut betscale = 0.0
    if self.f >= 0.0 {
      let k2 = sq(sbet1) * self.ep2
      let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2)
      lamscale = self.f * cbet1 * self.a3f(eps) * pi
      betscale = lamscale * cbet1
      x = lam12x / lamscale
      y = sbet12a / betscale
    } else {
      let cbet12a = cbet2 * cbet1 - sbet2 * sbet1
      let bet12a = @math.atan2(sbet12a, cbet12a)
      let (_, m12b, m0, _, _) = self.lengths(
        self.n,
        pi + bet12a,
        sbet1,
        -cbet1,
        dn1,
        sbet2,
        cbet2,
        dn2,
        cbet1,
        cbet2,
        m_reducedlength,
        c1a,
        c2a,
      )
      x = -1.0 + m12b / (cbet1 * cbet2 * m0 * pi)
      betscale = if x < -0.01 { sbet12a / x } else { -self.f * sq(cbet1) * pi }
      lamscale = betscale / cbet1
      y = lam12x / lamscale
    }
    if y > -tol1 && x > -1.0 - xthresh {
      if self.f >= 0.0 {
        salp1 = fmin(1.0, -x)
        calp1 = -(1.0 - sq(salp1)).sqrt()
      } else {
        calp1 = fmax(if x > -tol1 { 0.0 } else { -1.0 }, x)
        salp1 = (1.0 - sq(calp1)).sqrt()
      }
    } else {
      let k = astroid(x, y)
      let omg12a = lamscale *
        (if self.f >= 0.0 { -x * k / (1.0 + k) } else { -y * (1.0 + k) / k })
      somg12 = @math.sin(omg12a)
      comg12 = -@math.cos(omg12a)
      salp1 = cbet2 * somg12
      calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12)
    }
  }
  if !(salp1 <= 0.0) {
    let (s, c) = norm2(salp1, calp1)
    salp1 = s
    calp1 = c
  } else {
    salp1 = 1.0
    calp1 = 0.0
  }
  (sig12, salp1, calp1, salp2, calp2, dnm)
}

///|
priv struct Lambda12Result {
  lam12 : Double
  salp2 : Double
  calp2 : Double
  sig12 : Double
  ssig1 : Double
  csig1 : Double
  ssig2 : Double
  csig2 : Double
  eps : Double
  domg12 : Double
  dlam12 : Double
}

///|
fn Geodesic::lambda12(
  self : Geodesic,
  sbet1 : Double,
  cbet1 : Double,
  dn1 : Double,
  sbet2 : Double,
  cbet2 : Double,
  dn2 : Double,
  salp1 : Double,
  calp1_in : Double,
  slam120 : Double,
  clam120 : Double,
  diffp : Bool,
  c1a : FixedArray[Double],
  c2a : FixedArray[Double],
  c3a : FixedArray[Double],
) -> Lambda12Result {
  let calp1 = if sbet1 == 0.0 && calp1_in == 0.0 { -tiny } else { calp1_in }
  let salp0 = salp1 * cbet1
  let calp0 = hypot(calp1, salp1 * sbet1)
  let somg1 = salp0 * sbet1
  let comg1 = calp1 * cbet1
  let (ssig1, csig1) = norm2(sbet1, calp1 * cbet1)
  let salp2 = if cbet2 != cbet1 { salp0 / cbet2 } else { salp1 }
  let calp2 = if cbet2 != cbet1 || sbet2.abs() != -sbet1 {
    (sq(calp1 * cbet1) +
    (if cbet1 < -sbet1 {
      (cbet2 - cbet1) * (cbet1 + cbet2)
    } else {
      (sbet1 - sbet2) * (sbet1 + sbet2)
    })).sqrt() /
    cbet2
  } else {
    calp1.abs()
  }
  let somg2 = salp0 * sbet2
  let comg2 = calp2 * cbet2
  let (ssig2, csig2) = norm2(sbet2, calp2 * cbet2)
  let sig12 = @math.atan2(
    fmax(0.0, csig1 * ssig2 - ssig1 * csig2) + 0.0,
    csig1 * csig2 + ssig1 * ssig2,
  )
  let somg12 = fmax(0.0, comg1 * somg2 - somg1 * comg2) + 0.0
  let comg12 = comg1 * comg2 + somg1 * somg2
  let eta = @math.atan2(
    somg12 * clam120 - comg12 * slam120,
    comg12 * clam120 + somg12 * slam120,
  )
  let k2 = sq(calp0) * self.ep2
  let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2)
  self.c3f(eps, c3a)
  let b312 = sin_cos_series(true, ssig2, csig2, c3a) -
    sin_cos_series(true, ssig1, csig1, c3a)
  let domg12 = -self.f * self.a3f(eps) * salp0 * (sig12 + b312)
  let lam12 = eta + domg12
  let dlam12 = if diffp {
    if calp2 == 0.0 {
      -2.0 * self.f1 * dn1 / sbet1
    } else {
      let (_, d, _, _, _) = self.lengths(
        eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, m_reducedlength,
        c1a, c2a,
      )
      d * self.f1 / (calp2 * cbet2)
    }
  } else {
    @double.not_a_number
  }
  {
    lam12,
    salp2,
    calp2,
    sig12,
    ssig1,
    csig1,
    ssig2,
    csig2,
    eps,
    domg12,
    dlam12,
  }
}

///|
priv struct GenInverse {
  a12 : Double
  s12 : Double
  salp1 : Double
  calp1 : Double
  salp2 : Double
  calp2 : Double
  m12 : Double
  mm12 : Double
  mm21 : Double
  ss12 : Double
}

///|
fn Geodesic::gen_inverse(
  self : Geodesic,
  lat1_in : Double,
  lon1 : Double,
  lat2_in : Double,
  lon2 : Double,
) -> GenInverse {
  let nan = @double.not_a_number
  let outmask = m_all
  let mut a12 = nan
  let mut s12 = nan
  let mut m12 = nan
  let mut mm12 = nan
  let mut mm21 = nan
  let mut ss12 = nan
  let (lon12_0, lon12s_0) = ang_diff(lon1, lon2)
  let mut lonsign = copysign(1.0, lon12_0)
  let lon12 = lonsign * lon12_0
  let mut lon12s = lonsign * lon12s_0
  let lam12 = radians(lon12)
  let (slam12, clam12) = sincosde(lon12, lon12s)
  lon12s = 180.0 - lon12 - lon12s
  let mut lat1 = ang_round(lat_fix(lat1_in))
  let mut lat2 = ang_round(lat_fix(lat2_in))
  let swapp = if lat1.abs() < lat2.abs() || lat2.is_nan() { -1.0 } else { 1.0 }
  if swapp < 0.0 {
    lonsign = -lonsign
    let t = lat1
    lat1 = lat2
    lat2 = t
  }
  let latsign = copysign(1.0, -lat1)
  lat1 = lat1 * latsign
  lat2 = lat2 * latsign
  let (sb1, cb1) = sincosd(lat1)
  let (sbet1, cbet1_0) = norm2(sb1 * self.f1, cb1)
  let cbet1 = fmax(tiny, cbet1_0)
  let (sb2, cb2) = sincosd(lat2)
  let (sbet2_0, cbet2_0) = norm2(sb2 * self.f1, cb2)
  let mut sbet2 = sbet2_0
  let mut cbet2 = fmax(tiny, cbet2_0)
  if cbet1 < -sbet1 {
    if cbet2 == cbet1 {
      sbet2 = copysign(sbet1, sbet2)
    }
  } else if sbet2.abs() == -sbet1 {
    cbet2 = cbet1
  }
  let dn1 = (1.0 + self.ep2 * sq(sbet1)).sqrt()
  let dn2 = (1.0 + self.ep2 * sq(sbet2)).sqrt()
  let c1a = FixedArray::make(order + 1, 0.0)
  let c2a = FixedArray::make(order + 1, 0.0)
  let c3a = FixedArray::make(order, 0.0)
  let mut meridian = lat1 == -90.0 || slam12 == 0.0
  let mut calp1 = 0.0
  let mut salp1 = 0.0
  let mut calp2 = 0.0
  let mut salp2 = 0.0
  let mut sig12 = 0.0
  let mut s12x = 0.0
  let mut m12x = 0.0
  if meridian {
    calp1 = clam12
    salp1 = slam12
    calp2 = 1.0
    salp2 = 0.0
    let ssig1 = sbet1
    let csig1 = calp1 * cbet1
    let ssig2 = sbet2
    let csig2 = calp2 * cbet2
    sig12 = @math.atan2(
      fmax(0.0, csig1 * ssig2 - ssig1 * csig2) + 0.0,
      csig1 * csig2 + ssig1 * ssig2,
    )
    let (sx, mx, _, m1, m2) = self.lengths(
      self.n,
      sig12,
      ssig1,
      csig1,
      dn1,
      ssig2,
      csig2,
      dn2,
      cbet1,
      cbet2,
      outmask | m_distance | m_reducedlength,
      c1a,
      c2a,
    )
    s12x = sx
    m12x = mx
    mm12 = m1
    mm21 = m2
    if sig12 < tol2 || m12x >= 0.0 {
      if sig12 < 3.0 * tiny || (sig12 < tol0 && (s12x < 0.0 || m12x < 0.0)) {
        sig12 = 0.0
        m12x = 0.0
        s12x = 0.0
      }
      m12x = m12x * self.b
      s12x = s12x * self.b
      a12 = degrees(sig12)
    } else {
      meridian = false
    }
  }
  let mut somg12 = 2.0
  let mut comg12 = 0.0
  let mut omg12 = 0.0
  if !meridian && sbet1 == 0.0 && (self.f <= 0.0 || lon12s >= self.f * 180.0) {
    calp1 = 0.0
    calp2 = 0.0
    salp1 = 1.0
    salp2 = 1.0
    s12x = self.a * lam12
    sig12 = lam12 / self.f1
    omg12 = sig12
    m12x = self.b * @math.sin(sig12)
    mm12 = @math.cos(sig12)
    mm21 = mm12
    a12 = lon12 / self.f1
  } else if !meridian {
    let (sg, sa1, ca1, sa2, ca2, dnm) = self.inverse_start(
      sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, c1a, c2a,
    )
    sig12 = sg
    salp1 = sa1
    calp1 = ca1
    salp2 = sa2
    calp2 = ca2
    if sig12 >= 0.0 {
      s12x = sig12 * self.b * dnm
      m12x = sq(dnm) * self.b * @math.sin(sig12 / dnm)
      mm12 = @math.cos(sig12 / dnm)
      mm21 = mm12
      a12 = degrees(sig12)
      omg12 = lam12 / (self.f1 * dnm)
    } else {
      let mut numit = 0
      let mut tripn = false
      let mut tripb = false
      let mut salp1a = tiny
      let mut calp1a = 1.0
      let mut salp1b = tiny
      let mut calp1b = -1.0
      let mut r : Lambda12Result = {
        lam12: 0.0,
        salp2: 0.0,
        calp2: 0.0,
        sig12: 0.0,
        ssig1: 0.0,
        csig1: 0.0,
        ssig2: 0.0,
        csig2: 0.0,
        eps: 0.0,
        domg12: 0.0,
        dlam12: 0.0,
      }
      while true {
        r = self.lambda12(
          sbet1,
          cbet1,
          dn1,
          sbet2,
          cbet2,
          dn2,
          salp1,
          calp1,
          slam12,
          clam12,
          numit < maxit1,
          c1a,
          c2a,
          c3a,
        )
        let v = r.lam12
        let dv = r.dlam12
        if tripb ||
          !(v.abs() >= (if tripn { 8.0 } else { 1.0 }) * tol0) ||
          numit == maxit2 {
          break
        }
        if v > 0.0 && (numit > maxit1 || calp1 / salp1 > calp1b / salp1b) {
          salp1b = salp1
          calp1b = calp1
        } else if v < 0.0 && (numit > maxit1 || calp1 / salp1 < calp1a / salp1a) {
          salp1a = salp1
          calp1a = calp1
        }
        numit += 1
        if numit < maxit1 && dv > 0.0 {
          let dalp1 = -v / dv
          if dalp1.abs() < pi {
            let sdalp1 = @math.sin(dalp1)
            let cdalp1 = @math.cos(dalp1)
            let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1
            if nsalp1 > 0.0 {
              calp1 = calp1 * cdalp1 - salp1 * sdalp1
              salp1 = nsalp1
              let (s, c) = norm2(salp1, calp1)
              salp1 = s
              calp1 = c
              tripn = v.abs() <= 16.0 * tol0
              continue
            }
          }
        }
        salp1 = (salp1a + salp1b) / 2.0
        calp1 = (calp1a + calp1b) / 2.0
        let (s, c) = norm2(salp1, calp1)
        salp1 = s
        calp1 = c
        tripn = false
        tripb = (salp1a - salp1).abs() + (calp1a - calp1) < tolb ||
          (salp1 - salp1b).abs() + (calp1 - calp1b) < tolb
      }
      salp2 = r.salp2
      calp2 = r.calp2
      sig12 = r.sig12
      let (sx, mx, _, m1, m2) = self.lengths(
        r.eps,
        sig12,
        r.ssig1,
        r.csig1,
        dn1,
        r.ssig2,
        r.csig2,
        dn2,
        cbet1,
        cbet2,
        outmask | m_distance,
        c1a,
        c2a,
      )
      mm12 = m1
      mm21 = m2
      m12x = mx * self.b
      s12x = sx * self.b
      a12 = degrees(sig12)
      let sdomg12 = @math.sin(r.domg12)
      let cdomg12 = @math.cos(r.domg12)
      somg12 = slam12 * cdomg12 - clam12 * sdomg12
      comg12 = clam12 * cdomg12 + slam12 * sdomg12
    }
  }
  s12 = 0.0 + s12x
  m12 = 0.0 + m12x
  // Area
  let salp0 = salp1 * cbet1
  let calp0 = hypot(calp1, salp1 * sbet1)
  if calp0 != 0.0 && salp0 != 0.0 {
    let (ssig1, csig1) = norm2(sbet1, calp1 * cbet1)
    let (ssig2, csig2) = norm2(sbet2, calp2 * cbet2)
    let k2 = sq(calp0) * self.ep2
    let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2)
    let a4 = sq(self.a) * calp0 * salp0 * self.e2
    let c4a = FixedArray::make(order, 0.0)
    self.c4f(eps, c4a)
    let b41 = sin_cos_series(false, ssig1, csig1, c4a)
    let b42 = sin_cos_series(false, ssig2, csig2, c4a)
    ss12 = a4 * (b42 - b41)
  } else {
    ss12 = 0.0
  }
  if !meridian && somg12 == 2.0 {
    somg12 = @math.sin(omg12)
    comg12 = @math.cos(omg12)
  }
  let alp12 = if !meridian && comg12 > -0.7071 && sbet2 - sbet1 < 1.75 {
    let domg12 = 1.0 + comg12
    let dbet1 = 1.0 + cbet1
    let dbet2 = 1.0 + cbet2
    2.0 *
    @math.atan2(
      somg12 * (sbet1 * dbet2 + sbet2 * dbet1),
      domg12 * (sbet1 * sbet2 + dbet1 * dbet2),
    )
  } else {
    let mut salp12 = salp2 * calp1 - calp2 * salp1
    let mut calp12 = calp2 * calp1 + salp2 * salp1
    if salp12 == 0.0 && calp12 < 0.0 {
      salp12 = tiny * calp1
      calp12 = -1.0
    }
    @math.atan2(salp12, calp12)
  }
  ss12 = ss12 + self.c2 * alp12
  ss12 = ss12 * swapp * lonsign * latsign
  ss12 = ss12 + 0.0
  if swapp < 0.0 {
    let t1 = salp1
    salp1 = salp2
    salp2 = t1
    let t2 = calp1
    calp1 = calp2
    calp2 = t2
    let t3 = mm12
    mm12 = mm21
    mm21 = t3
  }
  salp1 = salp1 * swapp * lonsign
  calp1 = calp1 * swapp * latsign
  salp2 = salp2 * swapp * lonsign
  calp2 = calp2 * swapp * latsign
  { a12, s12, salp1, calp1, salp2, calp2, m12, mm12, mm21, ss12, }
}

///|
/// Result of the inverse geodesic problem.
pub(all) struct InverseResult {
  lat1 : Double
  lon1 : Double
  lat2 : Double
  lon2 : Double
  /// Distance between the points (m).
  s12 : Double
  /// Azimuth at point 1 (degrees clockwise from north).
  azi1 : Double
  /// Forward azimuth at point 2 (degrees).
  azi2 : Double
  /// Arc length on the auxiliary sphere (degrees).
  a12 : Double
  /// Reduced length m12 (m).
  m12 : Double
  /// Geodesic scale M12.
  scale12 : Double
  /// Geodesic scale M21.
  scale21 : Double
  /// Area between the geodesic and the equator, S12 (m²).
  area12 : Double
} derive(Debug)

///|
pub extend InverseResult with @debug.Debug::{to_repr}

///|
/// Solve the inverse geodesic problem: the shortest path between two points.
/// With `long_unroll=true`, `lon2` is reported as `lon1` plus the signed
/// longitude difference instead of being reduced to [-180, 180].
pub fn Geodesic::inverse(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
  long_unroll? : Bool = false,
) -> InverseResult {
  let r = self.gen_inverse(lat1, lon1, lat2, lon2)
  let lon2o = if long_unroll {
    let (lon12, e) = ang_diff(lon1, lon2)
    lon1 + lon12 + e
  } else {
    ang_normalize(lon2)
  }
  {
    lat1: lat_fix(lat1),
    lon1: if long_unroll {
      lon1
    } else {
      ang_normalize(lon1)
    },
    lat2: lat_fix(lat2),
    lon2: lon2o,
    s12: r.s12,
    azi1: atan2d(r.salp1, r.calp1),
    azi2: atan2d(r.salp2, r.calp2),
    a12: r.a12,
    m12: r.m12,
    scale12: r.mm12,
    scale21: r.mm21,
    area12: r.ss12,
  }
}

///|
/// Shortest distance between two points in metres.
pub fn Geodesic::distance(
  self : Geodesic,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> Double {
  self.gen_inverse(lat1, lon1, lat2, lon2).s12
}
