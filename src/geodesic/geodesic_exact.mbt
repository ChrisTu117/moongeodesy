// MoonBit port of GeographicLib::GeodesicExact and GeodesicLineExact
// (Charles Karney, MIT/X11 License). Elliptic-integral geodesics, accurate
// for any flattening with a > 0 and b > 0, including |f| > 1/50 where the
// series in `Geodesic` lose digits. Transcription of GeodesicExact.cpp and
// GeodesicLineExact.cpp from GeographicLib 2.5.

///|
fn dn_beta(
  f : Double,
  f1 : Double,
  ep2 : Double,
  e2 : Double,
  sbet : Double,
  cbet : Double,
) -> Double {
  if f >= 0.0 {
    (1.0 + ep2 * @geomath.sq(sbet)).sqrt()
  } else {
    (1.0 - e2 * @geomath.sq(cbet)).sqrt() / f1
  }
}

///|
fn i4_asinhsqrt(x : Double) -> Double {
  if x == 0.0 {
    1.0
  } else if x > 0.0 {
    @math.asinh(x.sqrt()) / x.sqrt()
  } else {
    @math.asin((-x).sqrt()) / (-x).sqrt()
  }
}

///|
fn i4_t(x : Double) -> Double {
  x + ((1.0 + x).sqrt() * i4_asinhsqrt(x) - 1.0)
}

///|
fn i4_td(x : Double) -> Double {
  if x == 0.0 {
    4.0 / 3.0
  } else {
    1.0 + (1.0 - i4_asinhsqrt(x) / (1.0 + x).sqrt()) / (2.0 * x)
  }
}

///|
priv struct I4 {
  x : Double
  tx : Double
  tdx : Double
  k2 : Double
  sx : Double
  sx1 : Double
  sxx1 : Double
  asinhsx : Double
}

///|
fn I4::new(ep2 : Double, k2 : Double) -> I4 {
  let sx = ep2.abs().sqrt()
  let sx1 = (1.0 + ep2).sqrt()
  {
    x: ep2,
    tx: i4_t(ep2),
    tdx: i4_td(ep2),
    k2,
    sx,
    sx1,
    sxx1: sx * sx1,
    asinhsx: if ep2 > 0.0 {
      @math.asinh(sx)
    } else {
      @math.asin(sx)
    },
  }
}

///|
fn I4::dtx(self : I4, y : Double) -> Double {
  if self.x == y {
    self.tdx
  } else if self.x * y <= 0.0 {
    (self.tx - i4_t(y)) / (self.x - y)
  } else {
    let sy = y.abs().sqrt()
    let sy1 = (1.0 + y).sqrt()
    let z = (self.x - y) / (self.sx * sy1 + sy * self.sx1)
    let d1 = 2.0 * self.sx * sy
    let d2 = 2.0 * (self.x * sy * sy1 + y * self.sxx1)
    if self.x > 0.0 {
      1.0 + @math.asinh(z) / z / d1 - (self.asinhsx + @math.asinh(sy)) / d2
    } else {
      1.0 - @math.asin(z) / z / d1 - (self.asinhsx + @math.asin(sy)) / d2
    }
  }
}

///|
fn I4::eval(self : I4, sig : Double) -> Double {
  let ssig = @math.sin(sig)
  -self.dtx(self.k2 * @geomath.sq(ssig)) * ssig / 2.0
}

///|
fn i4_coeffs(ep2 : Double, k2 : Double, n_c4 : Int32) -> Array[Double] {
  let i4 = I4::new(ep2, k2)
  dst_transform(n_c4, fn(sig) { i4.eval(sig) })
}

///|
/// An ellipsoid on which geodesic problems are solved with elliptic integrals.
/// Use this when `|f|` is larger than about 1/50. `Geodesic` is the series
/// solution and is the one to use for WGS-84, GRS80, and CGCS2000.
pub struct GeodesicExact {
  a : Double
  f : Double
  f1 : Double
  e2 : Double
  ep2 : Double
  n : Double
  b : Double
  c2 : Double
  etol2 : Double
  n_c4 : Int32
}

///|
/// Equatorial radius `a` (metres) and flattening `f`. Raises
/// `InvalidEllipsoid` when `a` or the polar semi-axis `a * (1 - f)` is not
/// finite and positive. The elliptic integrals stay accurate for large `|f|`.
pub fn GeodesicExact::new(
  a : Double,
  f : Double,
) -> GeodesicExact raise GeodesicError {
  let f1 = 1.0 - f
  let e2 = f * (2.0 - f)
  let ep2 = e2 / @geomath.sq(f1)
  let n = f / (2.0 - f)
  let b = a * f1
  let c2 = if f == 0.0 {
    @geomath.sq(a)
  } else {
    let q = if f > 0.0 {
      @math.asinh(ep2.sqrt())
    } else {
      @math.atan((-e2).sqrt())
    }
    (@geomath.sq(a) + @geomath.sq(b) * (q / e2.abs().sqrt())) / 2.0
  }
  let etol2 = 0.1 *
    tol2 /
    (@geomath.fmax(0.001, f.abs()) * @geomath.fmin(1.0, 1.0 - f / 2.0) / 2.0).sqrt()
  if !(@geomath.is_finite(a) && a > 0.0) {
    raise InvalidEllipsoid("Equatorial radius is not positive")
  }
  if !(@geomath.is_finite(b) && b > 0.0) {
    raise InvalidEllipsoid("Polar semi-axis is not positive")
  }
  { a, f, f1, e2, ep2, n, b, c2, etol2, n_c4: c4_nterms(n), }
}

///|
let wgs84_exact : GeodesicExact = GeodesicExact::new(wgs84_a, wgs84_f) catch {
  _ => abort("unreachable")
}

///|
/// WGS-84 solved with elliptic integrals.
pub fn GeodesicExact::wgs84() -> GeodesicExact {
  wgs84_exact
}

///|
let grs80_exact : GeodesicExact = GeodesicExact::new(
  6378137.0,
  1.0 / 298.257222101,
) catch {
  _ => abort("unreachable")
}

///|
/// The GRS80 ellipsoid.
pub fn GeodesicExact::grs80() -> GeodesicExact {
  grs80_exact
}

///|
/// The CGCS2000 ellipsoid. Same `a` and `f` as GRS80.
pub fn GeodesicExact::cgcs2000() -> GeodesicExact {
  grs80_exact
}

///|
/// Equatorial radius `a` in metres.
pub fn GeodesicExact::equatorial_radius(self : GeodesicExact) -> Double {
  self.a
}

///|
/// Flattening `f`.
pub fn GeodesicExact::flattening(self : GeodesicExact) -> Double {
  self.f
}

///|
/// Total area of the ellipsoid in m².
pub fn GeodesicExact::ellipsoid_area(self : GeodesicExact) -> Double {
  4.0 * @geomath.pi * self.c2
}

///|
priv struct ExactLen {
  s12b : Double
  m12b : Double
  m12 : Double
  m21 : Double
}

///|
fn GeodesicExact::lengths(
  self : GeodesicExact,
  e : Elliptic,
  sig12 : Double,
  ssig1 : Double,
  csig1 : Double,
  dn1 : Double,
  ssig2 : Double,
  csig2 : Double,
  dn2 : Double,
  cbet1 : Double,
  cbet2 : Double,
) -> ExactLen {
  let half = @geomath.pi / 2.0
  let s12b = e.ec /
    half *
    (sig12 + (e.delta_e(ssig2, csig2, dn2) - e.delta_e(ssig1, csig1, dn1)))
  let m0x = -e.k2 * e.dc / half
  let j12 = m0x *
    (sig12 + (e.delta_d(ssig2, csig2, dn2) - e.delta_d(ssig1, csig1, dn1)))
  let m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * j12
  let csig12 = csig1 * csig2 + ssig1 * ssig2
  let t = self.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
  let m12 = csig12 + (t * ssig2 - csig2 * j12) * ssig1 / dn1
  let m21 = csig12 - (t * ssig1 - csig1 * j12) * ssig2 / dn2
  { s12b, m12b, m12, m21, }
}

///|
priv struct ExactStart {
  sig12 : Double
  salp1 : Double
  calp1 : Double
  salp2 : Double
  calp2 : Double
  dnm : Double
}

///|
fn GeodesicExact::inverse_start(
  self : GeodesicExact,
  e0 : Elliptic,
  sbet1 : Double,
  cbet1 : Double,
  dn1 : Double,
  sbet2 : Double,
  cbet2 : Double,
  dn2 : Double,
  lam12 : Double,
  slam12 : Double,
  clam12 : Double,
) -> ExactStart {
  let nan = @double.not_a_number
  let mut sig12 = -1.0
  let mut salp2 = nan
  let mut calp2 = nan
  let mut dnm = nan
  let sbet12 = sbet2 * cbet1 - cbet2 * sbet1
  let cbet12 = cbet2 * cbet1 + sbet2 * sbet1
  let sbet12a = sbet2 * cbet1 + cbet2 * sbet1
  let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5
  let mut somg12 = 0.0
  let mut comg12 = 0.0
  if shortline {
    let mut sbetm2 = @geomath.sq(sbet1 + sbet2)
    sbetm2 = sbetm2 / (sbetm2 + @geomath.sq(cbet1 + cbet2))
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
    sbet12 + cbet2 * sbet1 * @geomath.sq(somg12) / (1.0 + comg12)
  } else {
    sbet12a - cbet2 * sbet1 * @geomath.sq(somg12) / (1.0 - comg12)
  }
  let ssig12 = @geomath.hypot(salp1, calp1)
  let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12
  if shortline && ssig12 < self.etol2 {
    salp2 = cbet1 * somg12
    calp2 = sbet12 -
      cbet1 *
      sbet2 *
      (if comg12 >= 0.0 {
        @geomath.sq(somg12) / (1.0 + comg12)
      } else {
        1.0 - comg12
      })
    let (s, c) = @geomath.norm2(salp2, calp2)
    salp2 = s
    calp2 = c
    sig12 = @math.atan2(ssig12, csig12)
  } else if self.n.abs() > 0.1 ||
    csig12 >= 0.0 ||
    ssig12 >= 6.0 * self.n.abs() * @geomath.pi * @geomath.sq(cbet1) {
    ()
  } else {
    let lam12x = @math.atan2(-slam12, -clam12)
    let mut x = 0.0
    let mut y = 0.0
    let mut lamscale = 0.0
    let mut betscale = 0.0
    if self.f >= 0.0 {
      let k2 = @geomath.sq(sbet1) * self.ep2
      let e = Elliptic::from_parts(-k2, -self.ep2, 1.0 + k2, 1.0 + self.ep2)
      lamscale = self.e2 / self.f1 * cbet1 * 2.0 * e.hc
      betscale = lamscale * cbet1
      x = lam12x / lamscale
      y = sbet12a / betscale
    } else {
      let cbet12a = cbet2 * cbet1 - sbet2 * sbet1
      let bet12a = @math.atan2(sbet12a, cbet12a)
      let len = self.lengths(
        e0,
        @geomath.pi + bet12a,
        sbet1,
        -cbet1,
        dn1,
        sbet2,
        cbet2,
        dn2,
        cbet1,
        cbet2,
      )
      // m0 is the secular coefficient from Lengths: -k2 * D / (pi/2).
      let m0 = -e0.k2 * e0.dc / (@geomath.pi / 2.0)
      x = -1.0 + len.m12b / (cbet1 * cbet2 * m0 * @geomath.pi)
      betscale = if x < -0.01 {
        sbet12a / x
      } else {
        -self.f * @geomath.sq(cbet1) * @geomath.pi
      }
      lamscale = betscale / cbet1
      y = lam12x / lamscale
    }
    if y > -tol1 && x > -1.0 - xthresh {
      if self.f >= 0.0 {
        salp1 = @geomath.fmin(1.0, -x)
        calp1 = -(1.0 - @geomath.sq(salp1)).sqrt()
      } else {
        calp1 = @geomath.fmax(if x > -tol1 { 0.0 } else { -1.0 }, x)
        salp1 = (1.0 - @geomath.sq(calp1)).sqrt()
      }
    } else {
      let k = astroid(x, y)
      let omg12a = lamscale *
        (if self.f >= 0.0 { -x * k / (1.0 + k) } else { -y * (1.0 + k) / k })
      somg12 = @math.sin(omg12a)
      comg12 = -@math.cos(omg12a)
      salp1 = cbet2 * somg12
      calp1 = sbet12a - cbet2 * sbet1 * @geomath.sq(somg12) / (1.0 - comg12)
    }
  }
  if !(salp1 <= 0.0) {
    let (s, c) = @geomath.norm2(salp1, calp1)
    salp1 = s
    calp1 = c
  } else {
    salp1 = 1.0
    calp1 = 0.0
  }
  { sig12, salp1, calp1, salp2, calp2, dnm, }
}

///|
priv struct ExactLam {
  lam12 : Double
  salp2 : Double
  calp2 : Double
  sig12 : Double
  ssig1 : Double
  csig1 : Double
  ssig2 : Double
  csig2 : Double
  domg12 : Double
  dlam12 : Double
}

///|
fn GeodesicExact::lambda12(
  self : GeodesicExact,
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
) -> ExactLam {
  let calp1 = if sbet1 == 0.0 && calp1_in == 0.0 { -tiny } else { calp1_in }
  let salp0 = salp1 * cbet1
  let calp0 = @geomath.hypot(calp1, salp1 * sbet1)
  let somg1 = salp0 * sbet1
  let comg1 = calp1 * cbet1
  let cchi1 = self.f1 * dn1 * comg1
  let (ssig1, csig1) = @geomath.norm2(sbet1, comg1)
  let salp2 = if cbet2 != cbet1 { salp0 / cbet2 } else { salp1 }
  let calp2 = if cbet2 != cbet1 || sbet2.abs() != -sbet1 {
    (@geomath.sq(calp1 * cbet1) +
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
  let cchi2 = self.f1 * dn2 * comg2
  let (ssig2, csig2) = @geomath.norm2(sbet2, comg2)
  let sig12 = @math.atan2(
    @geomath.fmax(0.0, csig1 * ssig2 - ssig1 * csig2),
    csig1 * csig2 + ssig1 * ssig2,
  )
  let somg12 = @geomath.fmax(0.0, comg1 * somg2 - somg1 * comg2)
  let comg12 = comg1 * comg2 + somg1 * somg2
  let k2 = @geomath.sq(calp0) * self.ep2
  let e = Elliptic::from_parts(-k2, -self.ep2, 1.0 + k2, 1.0 + self.ep2)
  let schi12 = @geomath.fmax(0.0, cchi1 * somg2 - somg1 * cchi2)
  let cchi12 = cchi1 * cchi2 + somg1 * somg2
  let eta = @math.atan2(
    schi12 * clam120 - cchi12 * slam120,
    cchi12 * clam120 + schi12 * slam120,
  )
  let half = @geomath.pi / 2.0
  let deta12 = -self.e2 /
    self.f1 *
    salp0 *
    e.hc /
    half *
    (sig12 + (e.delta_h(ssig2, csig2, dn2) - e.delta_h(ssig1, csig1, dn1)))
  let lam12 = eta + deta12
  let domg12 = deta12 +
    @math.atan2(
      schi12 * comg12 - cchi12 * somg12,
      cchi12 * comg12 + schi12 * somg12,
    )
  let dlam12 = if calp2 == 0.0 {
    -2.0 * self.f1 * dn1 / sbet1
  } else {
    let len = self.lengths(
      e, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
    )
    len.m12b * self.f1 / (calp2 * cbet2)
  }
  { lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, domg12, dlam12, }
}

///|
fn GeodesicExact::gen_inverse(
  self : GeodesicExact,
  lat1_in : Double,
  lon1 : Double,
  lat2_in : Double,
  lon2 : Double,
) -> GenInverse {
  let nan = @double.not_a_number
  let mut a12 = nan
  let mut s12 = nan
  let mut m12 = nan
  let mut mm12 = nan
  let mut mm21 = nan
  let mut ss12 = nan
  let (lon12_0, lon12s_0) = @geomath.ang_diff(lon1, lon2)
  let mut lonsign = @geomath.copysign(1.0, lon12_0)
  let lon12 = lonsign * lon12_0
  let mut lon12s = lonsign * lon12s_0
  let lam12 = @geomath.radians(lon12)
  let (slam12, clam12) = @geomath.sincosde(lon12, lon12s)
  lon12s = 180.0 - lon12 - lon12s
  let mut lat1 = @geomath.ang_round(@geomath.lat_fix(lat1_in))
  let mut lat2 = @geomath.ang_round(@geomath.lat_fix(lat2_in))
  let swapp = if lat1.abs() < lat2.abs() || lat2.is_nan() { -1.0 } else { 1.0 }
  if swapp < 0.0 {
    lonsign = -lonsign
    let t = lat1
    lat1 = lat2
    lat2 = t
  }
  let latsign = @geomath.copysign(1.0, -lat1)
  lat1 = lat1 * latsign
  lat2 = lat2 * latsign
  let (sb1, cb1) = @geomath.sincosd(lat1)
  let (sbet1, cbet1_0) = @geomath.norm2(sb1 * self.f1, cb1)
  let cbet1 = @geomath.fmax(tiny, cbet1_0)
  let (sb2, cb2) = @geomath.sincosd(lat2)
  let (sbet2_0, cbet2_0) = @geomath.norm2(sb2 * self.f1, cb2)
  let mut sbet2 = sbet2_0
  let mut cbet2 = @geomath.fmax(tiny, cbet2_0)
  if cbet1 < -sbet1 {
    if cbet2 == cbet1 {
      sbet2 = @geomath.copysign(sbet1, sbet2)
    }
  } else if sbet2.abs() == -sbet1 {
    cbet2 = cbet1
  }
  let dn1 = dn_beta(self.f, self.f1, self.ep2, self.e2, sbet1, cbet1)
  let dn2 = dn_beta(self.f, self.f1, self.ep2, self.e2, sbet2, cbet2)
  let e0 = Elliptic::from_k2(-self.ep2)
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
      @geomath.fmax(0.0, csig1 * ssig2 - ssig1 * csig2),
      csig1 * csig2 + ssig1 * ssig2,
    )
    let len = self.lengths(
      e0, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
    )
    s12x = len.s12b
    m12x = len.m12b
    mm12 = len.m12
    mm21 = len.m21
    if sig12 < 1.0 || m12x >= 0.0 {
      if sig12 < 3.0 * tiny || (sig12 < tol0 && (s12x < 0.0 || m12x < 0.0)) {
        sig12 = 0.0
        m12x = 0.0
        s12x = 0.0
      }
      m12x = m12x * self.b
      s12x = s12x * self.b
      a12 = @geomath.degrees(sig12)
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
    let st = self.inverse_start(
      e0, sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12,
    )
    sig12 = st.sig12
    salp1 = st.salp1
    calp1 = st.calp1
    salp2 = st.salp2
    calp2 = st.calp2
    if sig12 >= 0.0 {
      s12x = sig12 * self.b * st.dnm
      m12x = @geomath.sq(st.dnm) * self.b * @math.sin(sig12 / st.dnm)
      mm12 = @math.cos(sig12 / st.dnm)
      mm21 = mm12
      a12 = @geomath.degrees(sig12)
      omg12 = lam12 / (self.f1 * st.dnm)
    } else {
      let mut numit : Int32 = 0
      let mut tripn = false
      let mut tripb = false
      let mut salp1a = tiny
      let mut calp1a = 1.0
      let mut salp1b = tiny
      let mut calp1b = -1.0
      let mut r : ExactLam = {
        lam12: 0.0,
        salp2: 0.0,
        calp2: 0.0,
        sig12: 0.0,
        ssig1: 0.0,
        csig1: 0.0,
        ssig2: 0.0,
        csig2: 0.0,
        domg12: 0.0,
        dlam12: 0.0,
      }
      while true {
        r = self.lambda12(
          sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, slam12, clam12,
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
        if numit < maxit1 && dv > 0.0 {
          let dalp1 = -v / dv
          if dalp1.abs() < @geomath.pi {
            let sdalp1 = @math.sin(dalp1)
            let cdalp1 = @math.cos(dalp1)
            let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1
            if nsalp1 > 0.0 {
              calp1 = calp1 * cdalp1 - salp1 * sdalp1
              salp1 = nsalp1
              let (s, c) = @geomath.norm2(salp1, calp1)
              salp1 = s
              calp1 = c
              tripn = v.abs() <= 16.0 * tol0
              numit = numit + 1
              continue
            }
          }
        }
        salp1 = (salp1a + salp1b) / 2.0
        calp1 = (calp1a + calp1b) / 2.0
        let (s, c) = @geomath.norm2(salp1, calp1)
        salp1 = s
        calp1 = c
        tripn = false
        tripb = (salp1a - salp1).abs() + (calp1a - calp1) < tolb ||
          (salp1 - salp1b).abs() + (calp1 - calp1b) < tolb
        numit = numit + 1
      }
      salp2 = r.salp2
      calp2 = r.calp2
      sig12 = r.sig12
      let k2 = @geomath.sq(@geomath.hypot(calp1, salp1 * sbet1)) * self.ep2
      let e = Elliptic::from_parts(-k2, -self.ep2, 1.0 + k2, 1.0 + self.ep2)
      let len = self.lengths(
        e,
        sig12,
        r.ssig1,
        r.csig1,
        dn1,
        r.ssig2,
        r.csig2,
        dn2,
        cbet1,
        cbet2,
      )
      mm12 = len.m12
      mm21 = len.m21
      m12x = len.m12b * self.b
      s12x = len.s12b * self.b
      a12 = @geomath.degrees(sig12)
      let sdomg12 = @math.sin(r.domg12)
      let cdomg12 = @math.cos(r.domg12)
      somg12 = slam12 * cdomg12 - clam12 * sdomg12
      comg12 = clam12 * cdomg12 + slam12 * sdomg12
    }
  }
  s12 = 0.0 + s12x
  m12 = 0.0 + m12x
  let salp0 = salp1 * cbet1
  let calp0 = @geomath.hypot(calp1, salp1 * sbet1)
  let a4 = @geomath.sq(self.a) * calp0 * salp0 * self.e2
  ss12 = if a4 != 0.0 {
    let (ssig1, csig1) = @geomath.norm2(sbet1, calp1 * cbet1)
    let (ssig2, csig2) = @geomath.norm2(sbet2, calp2 * cbet2)
    let k2 = @geomath.sq(calp0) * self.ep2
    let c4a = i4_coeffs(self.ep2, k2, self.n_c4)
    a4 * dst_integral_diff(ssig1, csig1, ssig2, csig2, c4a)
  } else {
    0.0
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
/// Solve the inverse problem with elliptic integrals.
/// See `Geodesic::inverse` for units, poles, and `long_unroll`.
pub fn GeodesicExact::inverse(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
  long_unroll~ : Bool = false,
) -> InverseResult {
  let r = self.gen_inverse(lat1, lon1, lat2, lon2)
  let lon2o = if long_unroll {
    let (lon12, e) = @geomath.ang_diff(lon1, lon2)
    lon1 + lon12 + e
  } else {
    @geomath.ang_normalize(lon2)
  }
  {
    lat1: @geomath.lat_fix(lat1),
    lon1: if long_unroll {
      lon1
    } else {
      @geomath.ang_normalize(lon1)
    },
    lat2: @geomath.lat_fix(lat2),
    lon2: lon2o,
    s12: r.s12,
    azi1: @geomath.atan2d(r.salp1, r.calp1),
    azi2: @geomath.atan2d(r.salp2, r.calp2),
    a12: r.a12,
    m12: r.m12,
    scale12: r.mm12,
    scale21: r.mm21,
    area12: r.ss12,
  }
}

///|
/// Shortest distance in metres, from elliptic integrals.
pub fn GeodesicExact::distance(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> Double {
  self.gen_inverse(lat1, lon1, lat2, lon2).s12
}

///|
/// A geodesic line evaluated with elliptic integrals.
pub struct GeodesicLineExact {
  b : Double
  c2 : Double
  f1 : Double
  e2 : Double
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
  cchi1 : Double
  k2 : Double
  ek2 : Double
  ekp2 : Double
  ealphap2 : Double
  eeps : Double
  eec : Double
  edc : Double
  ehc : Double
  ee0 : Double
  ee1 : Double
  stau1 : Double
  ctau1 : Double
  dd0 : Double
  dd1 : Double
  hh0 : Double
  hh1 : Double
  a4 : Double
  b41 : Double
  c4a : Array[Double]
  mut s13 : Double
  mut a13 : Double
}

///|
fn GeodesicLineExact::ell(self : GeodesicLineExact) -> Elliptic {
  {
    k2: self.ek2,
    kp2: self.ekp2,
    alphap2: self.ealphap2,
    eps: self.eeps,
    ec: self.eec,
    dc: self.edc,
    hc: self.ehc,
  }
}

///|
fn GeodesicLineExact::make(
  geod : GeodesicExact,
  lat1_in : Double,
  lon1 : Double,
  azi1_in : Double,
  salp1_in : Double,
  calp1_in : Double,
) -> GeodesicLineExact {
  let lat1 = @geomath.lat_fix(lat1_in)
  let (azi1, salp1, calp1) = if salp1_in.is_nan() || calp1_in.is_nan() {
    let (s, c) = @geomath.sincosd(@geomath.ang_round(azi1_in))
    (@geomath.ang_normalize(azi1_in), s, c)
  } else {
    (azi1_in, salp1_in, calp1_in)
  }
  let (sb, cb) = @geomath.sincosd(@geomath.ang_round(lat1))
  let (sbet1, cb1) = @geomath.norm2(sb * geod.f1, cb)
  let cbet1 = @geomath.fmax(tiny, cb1)
  let dn1 = dn_beta(geod.f, geod.f1, geod.ep2, geod.e2, sbet1, cbet1)
  let salp0 = salp1 * cbet1
  let calp0 = @geomath.hypot(calp1, salp1 * sbet1)
  let somg1 = salp0 * sbet1
  let c0 = if sbet1 != 0.0 || calp1 != 0.0 { cbet1 * calp1 } else { 1.0 }
  let (ssig1, csig1) = @geomath.norm2(sbet1, c0)
  let cchi1 = geod.f1 * dn1 * c0
  let k2 = @geomath.sq(calp0) * geod.ep2
  let e = Elliptic::from_parts(-k2, -geod.ep2, 1.0 + k2, 1.0 + geod.ep2)
  let half = @geomath.pi / 2.0
  let ee0 = e.ec / half
  let ee1 = e.delta_e(ssig1, csig1, dn1)
  let s = @math.sin(ee1)
  let c = @math.cos(ee1)
  let stau1 = ssig1 * c + csig1 * s
  let ctau1 = csig1 * c - ssig1 * s
  let dd0 = e.dc / half
  let dd1 = e.delta_d(ssig1, csig1, dn1)
  let hh0 = e.hc / half
  let hh1 = e.delta_h(ssig1, csig1, dn1)
  let a4 = @geomath.sq(geod.a) * calp0 * salp0 * geod.e2
  let (b41, c4a) = if a4 == 0.0 {
    (0.0, [])
  } else {
    let c4a = i4_coeffs(geod.ep2, k2, geod.n_c4)
    (dst_integral(ssig1, csig1, c4a), c4a)
  }
  {
    b: geod.b,
    c2: geod.c2,
    f1: geod.f1,
    e2: geod.e2,
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
    cchi1,
    k2,
    ek2: e.k2,
    ekp2: e.kp2,
    ealphap2: e.alphap2,
    eeps: e.eps,
    eec: e.ec,
    edc: e.dc,
    ehc: e.hc,
    ee0,
    ee1,
    stau1,
    ctau1,
    dd0,
    dd1,
    hh0,
    hh1,
    a4,
    b41,
    c4a,
    s13: @double.not_a_number,
    a13: @double.not_a_number,
  }
}

///|
fn GeodesicLineExact::gen_position(
  self : GeodesicLineExact,
  arcmode : Bool,
  s12_a12 : Double,
  long_unroll : Bool,
) -> DirectResult {
  let mut sig12 = 0.0
  let mut ssig12 = 0.0
  let mut csig12 = 0.0
  let mut e2ell = 0.0
  if arcmode {
    sig12 = @geomath.radians(s12_a12)
    let (s, c) = @geomath.sincosd(s12_a12)
    ssig12 = s
    csig12 = c
  } else {
    let tau12 = s12_a12 / (self.b * self.ee0)
    let s = @math.sin(tau12)
    let c = @math.cos(tau12)
    e2ell = -self
      .ell()
      .delta_einv(
        self.stau1 * c + self.ctau1 * s,
        self.ctau1 * c - self.stau1 * s,
      )
    sig12 = tau12 - (e2ell - self.ee1)
    ssig12 = @math.sin(sig12)
    csig12 = @math.cos(sig12)
  }
  let ssig2 = self.ssig1 * csig12 + self.csig1 * ssig12
  let mut csig2 = self.csig1 * csig12 - self.ssig1 * ssig12
  let ell = self.ell()
  let dn2 = ell.delta_amp(ssig2, csig2)
  if arcmode {
    e2ell = ell.delta_e(ssig2, csig2, dn2)
  }
  let ab1 = self.ee0 * (e2ell - self.ee1)
  let sbet2 = self.calp0 * ssig2
  let mut cbet2 = @geomath.hypot(self.salp0, self.calp0 * csig2)
  if cbet2 == 0.0 {
    cbet2 = tiny
    csig2 = tiny
  }
  let salp2 = self.salp0
  let calp2 = self.calp0 * csig2
  let s12 = if arcmode { self.b * (self.ee0 * sig12 + ab1) } else { s12_a12 }
  let somg2 = self.salp0 * ssig2
  let east = @geomath.copysign(1.0, self.salp0)
  let cchi2 = self.f1 * dn2 * csig2
  let chi12 = if long_unroll {
    east *
    (
      sig12 -
      (@math.atan2(ssig2, csig2) - @math.atan2(self.ssig1, self.csig1)) +
      (
        @math.atan2(east * somg2, cchi2) -
        @math.atan2(east * self.somg1, self.cchi1)
      )
    )
  } else {
    @math.atan2(
      somg2 * self.cchi1 - cchi2 * self.somg1,
      cchi2 * self.cchi1 + somg2 * self.somg1,
    )
  }
  let lam12 = chi12 -
    self.e2 /
    self.f1 *
    self.salp0 *
    self.hh0 *
    (sig12 + (ell.delta_h(ssig2, csig2, dn2) - self.hh1))
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
  let j12 = self.k2 *
    self.dd0 *
    (sig12 + (ell.delta_d(ssig2, csig2, dn2) - self.dd1))
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
  let b42 = if self.a4 == 0.0 {
    0.0
  } else {
    dst_integral(ssig2, csig2, self.c4a)
  }
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
pub fn GeodesicLineExact::position(
  self : GeodesicLineExact,
  s12 : Double,
  long_unroll~ : Bool = false,
) -> DirectResult {
  self.gen_position(false, s12, long_unroll)
}

///|
/// Position at arc length `a12` degrees from point 1.
pub fn GeodesicLineExact::arc_position(
  self : GeodesicLineExact,
  a12 : Double,
  long_unroll~ : Bool = false,
) -> DirectResult {
  self.gen_position(true, a12, long_unroll)
}

///|
/// Distance to reference point 3, in metres.
pub fn GeodesicLineExact::distance(self : GeodesicLineExact) -> Double {
  self.s13
}

///|
/// Arc length to reference point 3, in degrees.
pub fn GeodesicLineExact::arc(self : GeodesicLineExact) -> Double {
  self.a13
}

///|
/// Set reference point 3 to distance `s13` metres from point 1.
pub fn GeodesicLineExact::set_distance(
  self : GeodesicLineExact,
  s13 : Double,
) -> Unit {
  self.s13 = s13
  self.a13 = self.gen_position(false, s13, false).a12
}

///|
/// Set reference point 3 to arc length `a13` degrees from point 1.
pub fn GeodesicLineExact::set_arc(
  self : GeodesicLineExact,
  a13 : Double,
) -> Unit {
  self.a13 = a13
  self.s13 = self.gen_position(true, a13, false).s12
}

///|
/// Points spaced at most `max_step` metres apart from point 1 to point 3.
pub fn GeodesicLineExact::densify(
  self : GeodesicLineExact,
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
pub fn GeodesicExact::line(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
) -> GeodesicLineExact {
  GeodesicLineExact::make(
    self, lat1, lon1, azi1, @double.not_a_number, @double.not_a_number,
  )
}

///|
/// The geodesic line between two points.
pub fn GeodesicExact::inverse_line(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> GeodesicLineExact {
  let r = self.gen_inverse(lat1, lon1, lat2, lon2)
  let azi1 = @geomath.atan2d(r.salp1, r.calp1)
  let line = GeodesicLineExact::make(self, lat1, lon1, azi1, r.salp1, r.calp1)
  line.set_arc(r.a12)
  line
}

///|
/// A geodesic line with point 3 at distance `s12` from point 1.
pub fn GeodesicExact::direct_line(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  s12 : Double,
) -> GeodesicLineExact {
  let line = self.line(lat1, lon1, azi1)
  line.set_distance(s12)
  line
}

///|
/// Solve the direct problem with elliptic integrals.
/// See `Geodesic::inverse` for poles, NaN latitudes, and `long_unroll`.
pub fn GeodesicExact::direct(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  s12 : Double,
  long_unroll~ : Bool = false,
) -> DirectResult {
  self.line(lat1, lon1, azi1).gen_position(false, s12, long_unroll)
}

///|
/// Direct problem with the distance given as an arc length in degrees.
pub fn GeodesicExact::arc_direct(
  self : GeodesicExact,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  a12 : Double,
  long_unroll~ : Bool = false,
) -> DirectResult {
  self.line(lat1, lon1, azi1).gen_position(true, a12, long_unroll)
}
