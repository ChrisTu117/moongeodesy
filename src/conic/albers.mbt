// Albers equal-area conic projection on an ellipsoid, a port of
// GeographicLib::AlbersEqualArea (C. F. F. Karney, MIT/X11 License,
// https://geographiclib.sourceforge.io/). Like the conformal projection it
// is written with divided differences, so it stays accurate for close
// standard parallels and tends to the cylindrical and azimuthal equal-area
// projections in the limits.

///|
/// An Albers equal-area conic projection. Coordinates are relative to the
/// origin at latitude `origin_latitude()` and the central meridian passed
/// to `forward` / `reverse`; `ConicGrid` adds a false origin and false
/// easting/northing as in EPSG.
pub struct AlbersEqualArea {
  priv a : Double
  priv f : Double
  priv fm : Double
  priv e2 : Double
  priv e : Double
  priv e2m : Double
  priv qz : Double
  priv qx : Double
  priv sign : Double
  priv lat0 : Double
  priv k0 : Double
  priv n0 : Double
  priv m02 : Double
  priv nrho0 : Double
  priv k2 : Double
  priv txi0 : Double
  priv scxi0 : Double
  priv sxi0 : Double
}

///|
let tol : Double = @geomath.dbl_epsilon.sqrt()

///|
let tol0 : Double = tol * @geomath.dbl_epsilon.sqrt().sqrt()

///|
/// Ellipsoid quantities used before the projection is built.
priv struct AlbersEll {
  f : Double
  e2 : Double
  e : Double
  e2m : Double
  qz : Double
  qx : Double
}

///|
/// atanh(e x) / e for oblate, atan(sqrt(-e2) x) / sqrt(-e2) for prolate,
/// x for a sphere.
fn AlbersEll::atanhee(self : AlbersEll, x : Double) -> Double {
  if self.f > 0.0 {
    @math.atanh(self.e * x) / self.e
  } else if self.f < 0.0 {
    @math.atan(self.e * x) / self.e
  } else {
    x
  }
}

///|
fn albers_ell(f : Double) -> AlbersEll {
  let e2 = f * (2.0 - f)
  let e = e2.abs().sqrt()
  let e2m = 1.0 - e2
  let ell0 = { f, e2, e, e2m, qz: 0.0, qx: 0.0, }
  let qz = 1.0 + e2m * ell0.atanhee(1.0)
  { f, e2, e, e2m, qz, qx: qz / (2.0 * e2m), }
}

///|
/// (atanhee(x) - atanhee(y)) / (x - y).
fn AlbersEll::datanhee(self : AlbersEll, x : Double, y : Double) -> Double {
  let t = x - y
  let d = 1.0 - self.e2 * x * y
  if t == 0.0 {
    1.0 / d
  } else if x * y < 0.0 {
    (self.atanhee(x) - self.atanhee(y)) / t
  } else {
    self.atanhee(t / d) / t
  }
}

///|
/// (Datanhee(1, y) - Datanhee(1, x)) / (y - x).
fn AlbersEll::ddatanhee(self : AlbersEll, x0 : Double, y0 : Double) -> Double {
  let (x, y) = if y0 < x0 { (y0, x0) } else { (x0, y0) }
  let q1 = self.e2.abs()
  let q2 = (2.0 * self.e / self.e2m * (1.0 - x)).abs()
  if x <= 0.0 || !(@geomath.fmin(q1, q2) < 0.75) {
    self.ddatanhee0(x, y)
  } else if q1 < q2 {
    self.ddatanhee1(x, y)
  } else {
    self.ddatanhee2(x, y)
  }
}

///|
fn AlbersEll::ddatanhee0(self : AlbersEll, x : Double, y : Double) -> Double {
  (self.datanhee(1.0, y) - self.datanhee(x, y)) / (1.0 - x)
}

///|
/// Series in e2, for small e2.
fn AlbersEll::ddatanhee1(self : AlbersEll, x : Double, y : Double) -> Double {
  let mut s = 0.0
  let mut z = 1.0
  let mut k = 1.0
  let mut t = 0.0
  let mut c = 0.0
  let mut en = 1.0
  while true {
    t = y * t + z
    c = c + t
    z = z * x
    t = y * t + z
    c = c + t
    z = z * x
    k = k + 2.0
    en = en * self.e2
    let ds = en * c / k
    s = s + ds
    if !(ds.abs() > s.abs() * eps / 2.0) {
      break
    }
  }
  s
}

///|
/// Series in 1 - x and 1 - y, for x, y close to 1.
fn AlbersEll::ddatanhee2(self : AlbersEll, x : Double, y : Double) -> Double {
  let dx = 1.0 - x
  let dy = 1.0 - y
  let mut xy = 1.0
  let mut yy = 1.0
  let mut ee = self.e2 / sq(self.e2m)
  let mut s = ee
  let mut m = 1
  while true {
    let mut c = (m + 2).to_double()
    let mut t = c
    yy = yy * dy
    xy = dx * xy + yy
    ee = ee / -self.e2m
    if m % 2 == 0 {
      ee = ee * self.e2
    }
    let kmax = (m + 1) / 2
    for k = kmax - 1; k >= 0; k = k - 1 {
      c = c * ((k + 1) * (2 * (k + m - 2 * kmax) + 3)).to_double()
      c = c / ((kmax - k) * (2 * (kmax - k) + 1)).to_double()
      t = self.e2 * t + c
    }
    let ds = t * ee * xy / (m + 2).to_double()
    s = s + ds
    if !(ds.abs() > s.abs() * eps / 2.0) {
      break
    }
    m = m + 1
  }
  s
}

///|
/// atanh(sqrt(x)) / sqrt(x) - 1, accurate for small x.
fn atanhxm1(x : Double) -> Double {
  if x.abs() < 0.5 {
    // Series sum_{j >= 1} x^j / (2j + 1), Horner from the last term; 60
    // terms reach 2^-53 for |x| < 0.5.
    let mut s = 0.0
    let mut j = 60
    while j > 0 {
      s = x * s + 1.0 / (2 * j + 1).to_double()
      j = j - 1
    }
    x * s
  } else {
    let xs = x.abs().sqrt()
    (if x > 0.0 { @math.atanh(xs) } else { @math.atan(xs) }) / xs - 1.0
  }
}

///|
/// tan(xi), xi the authalic latitude, given tan(phi).
fn AlbersEll::txif(self : AlbersEll, tphi : Double) -> Double {
  let cphi = 1.0 / (1.0 + sq(tphi)).sqrt()
  let sphi = tphi * cphi
  let es1 = self.e2 * sphi
  let es2m1 = 1.0 - es1 * sphi
  let es2m1a = self.e2m * es2m1
  (tphi / es2m1 + self.atanhee(sphi) / cphi) /
  (((1.0 + es1) / es2m1a + self.datanhee(1.0, sphi)) *
  ((1.0 - es1) / es2m1a + self.datanhee(1.0, -sphi))).sqrt()
}

///|
/// Inverse of `txif` by Newton's method.
fn AlbersEll::tphif(self : AlbersEll, txi : Double) -> Double {
  let mut tphi = txi
  let stol = tol * @geomath.fmax(1.0, txi.abs())
  for _ in 0..<5 {
    let txia = self.txif(tphi)
    let tphi2 = sq(tphi)
    let scphi2 = 1.0 + tphi2
    let scterm = scphi2 / (1.0 + sq(txia))
    let dtphi = (txi - txia) *
      scterm *
      scterm.sqrt() *
      self.qx *
      sq(1.0 - self.e2 * tphi2 / scphi2)
    tphi = tphi + dtphi
    if !(dtphi.abs() >= stol) {
      break
    }
  }
  tphi
}

///|
/// An Albers equal-area conic projection with two standard parallels
/// `stdlat1` and `stdlat2` (degrees), on which the scale is `k1`. `a` is the
/// equatorial radius (m) and `f` the flattening.
pub fn AlbersEqualArea::new(
  a : Double,
  f : Double,
  stdlat1 : Double,
  stdlat2 : Double,
  k1~ : Double = 1.0,
) -> AlbersEqualArea raise ConicError {
  check_ellipsoid(a, f)
  check_scale(k1)
  check_stdlat(stdlat1)
  check_stdlat(stdlat2)
  let (s1, c1) = @geomath.sincosd(stdlat1)
  let (s2, c2) = @geomath.sincosd(stdlat2)
  if c1 == 0.0 && c2 == 0.0 && s1 * s2 <= 0.0 {
    raise IncompatibleStandardLatitudes(stdlat1, stdlat2)
  }
  albers_init(a, f, s1, c1, s2, c2, k1)
}

///|
/// An Albers equal-area conic projection with one standard parallel
/// `stdlat` (degrees), on which the scale is `k0`.
pub fn AlbersEqualArea::new_1sp(
  a : Double,
  f : Double,
  stdlat : Double,
  k0~ : Double = 1.0,
) -> AlbersEqualArea raise ConicError {
  check_ellipsoid(a, f)
  check_scale(k0)
  check_stdlat(stdlat)
  let (s, c) = @geomath.sincosd(stdlat)
  albers_init(a, f, s, c, s, c, k0)
}

///|
fn albers_init(
  a : Double,
  f : Double,
  sinlat1 : Double,
  coslat1 : Double,
  sinlat2 : Double,
  coslat2 : Double,
  k1 : Double,
) -> AlbersEqualArea {
  let ell = albers_ell(f)
  let fm = 1.0 - f
  let e2 = ell.e2
  let e2m = ell.e2m
  let qz = ell.qz
  let r1 = @geomath.hypot(sinlat1, coslat1)
  let r2 = @geomath.hypot(sinlat2, coslat2)
  let mut sphi1 = sinlat1 / r1
  let mut cphi1 = coslat1 / r1
  let mut sphi2 = sinlat2 / r2
  let mut cphi2 = coslat2 / r2
  let polar = cphi1 == 0.0
  cphi1 = @geomath.fmax(epsx, cphi1)
  cphi2 = @geomath.fmax(epsx, cphi2)
  let sign = if sphi1 + sphi2 >= 0.0 { 1.0 } else { -1.0 }
  sphi1 = sphi1 * sign
  sphi2 = sphi2 * sign
  if sphi1 > sphi2 {
    let ts = sphi1
    sphi1 = sphi2
    sphi2 = ts
    let tc = cphi1
    cphi1 = cphi2
    cphi2 = tc
  }
  let tphi1 = sphi1 / cphi1
  let tphi2 = sphi2 / cphi2
  let mut cc = 1.0
  let tphi0 = if polar || tphi1 == tphi2 {
    tphi2
  } else {
    let tbet1 = fm * tphi1
    let scbet12 = 1.0 + sq(tbet1)
    let tbet2 = fm * tphi2
    let scbet22 = 1.0 + sq(tbet2)
    let txi1 = ell.txif(tphi1)
    let cxi1 = 1.0 / hyp(txi1)
    let sxi1 = txi1 * cxi1
    let txi2 = ell.txif(tphi2)
    let cxi2 = 1.0 / hyp(txi2)
    let sxi2 = txi2 * cxi2
    let dtbet2 = fm * (tbet1 + tbet2)
    let es1 = 1.0 - e2 * sq(sphi1)
    let es2 = 1.0 - e2 * sq(sphi2)
    let dsxi = (
        (1.0 + e2 * sphi1 * sphi2) / (es2 * es1) + ell.datanhee(sphi2, sphi1)
      ) *
      dsn(tphi2, tphi1, sphi2, sphi1) /
      (2.0 * ell.qx)
    let den = (sxi2 + sxi1) * dtbet2 + (scbet22 + scbet12) * dsxi
    let s = 2.0 * dtbet2 / den
    let sm1 = -dsn(tphi2, tphi1, sphi2, sphi1) *
      (
        -((if sphi2 <= 0.0 {
          (1.0 - sxi2) / (1.0 - sphi2)
        } else {
          sq(cxi2 / cphi2) * (1.0 + sphi2) / (1.0 + sxi2)
        }) +
        (if sphi1 <= 0.0 {
          (1.0 - sxi1) / (1.0 - sphi1)
        } else {
          sq(cxi1 / cphi1) * (1.0 + sphi1) / (1.0 + sxi1)
        })) *
        (1.0 + e2 * (sphi1 + sphi2 + sphi1 * sphi2)) /
        (1.0 + (sphi1 + sphi2 + sphi1 * sphi2)) +
        (
          scbet22 *
          (if sphi2 <= 0.0 { 1.0 - sphi2 } else { sq(cphi2) / (1.0 + sphi2) }) +
          scbet12 *
          (if sphi1 <= 0.0 { 1.0 - sphi1 } else { sq(cphi1) / (1.0 + sphi1) })
        ) *
        (
          e2 * (1.0 + sphi1 + sphi2 + e2 * sphi1 * sphi2) / (es1 * es2) +
          e2m * ell.ddatanhee(sphi1, sphi2)
        ) /
        qz
      ) /
      den
    cc = den / (2.0 * scbet12 * scbet22 * dsxi)
    // Newton's method for tan(phi0), the latitude where the scale is k0.
    let mut tphi0 = (tphi2 + tphi1) / 2.0
    let stol = tol0 * @geomath.fmax(1.0, tphi0.abs())
    for _ in 0..<40 {
      let scphi02 = 1.0 + sq(tphi0)
      let scphi0 = scphi02.sqrt()
      let sphi0 = tphi0 / scphi0
      let sphi0m = 1.0 / (scphi0 * (tphi0 + scphi0))
      let g = (1.0 + sq(fm * tphi0)) * sphi0
      let dg = e2m * scphi02 * (1.0 + 2.0 * sq(tphi0)) + e2
      let d = sphi0m *
        (1.0 - e2 * (1.0 + 2.0 * sphi0 * (1.0 + sphi0))) /
        (e2m * (1.0 + sphi0))
      let dd = -2.0 *
        (1.0 - e2 * sq(sphi0) * (2.0 * sphi0 + 3.0)) /
        (e2m * sq(1.0 + sphi0))
      let aa = -e2 *
        sq(sphi0m) *
        (2.0 + (1.0 + e2) * sphi0) /
        (e2m * (1.0 - e2 * sq(sphi0)))
      let bb = sphi0m *
        e2m /
        (1.0 - e2 * sphi0) *
        (atanhxm1(e2 * sq(sphi0m / (1.0 - e2 * sphi0))) - e2 * sphi0m / e2m)
      let dab = 2.0 *
        e2 *
        (2.0 - e2 * (1.0 + sq(sphi0))) /
        (e2m * sq(1.0 - e2 * sq(sphi0)) * scphi02)
      let u = sm1 * g - s / qz * (d - g * (aa + bb))
      let du = sm1 * dg - s / qz * (dd - dg * (aa + bb) - g * dab)
      let dtu = -u / du * (scphi0 * scphi02)
      tphi0 = tphi0 + dtu
      if !(dtu.abs() >= stol) {
        break
      }
    }
    tphi0
  }
  let txi0 = ell.txif(tphi0)
  let scxi0 = hyp(txi0)
  let sxi0 = txi0 / scxi0
  let n0 = tphi0 / hyp(tphi0)
  let m02 = 1.0 / (1.0 + sq(fm * tphi0))
  let nrho0 = if polar { 0.0 } else { a * m02.sqrt() }
  let k0 = (if tphi1 == tphi2 { 1.0 } else { cc / (m02 + n0 * qz * sxi0) }).sqrt() *
    k1
  {
    a,
    f,
    fm,
    e2,
    e: ell.e,
    e2m,
    qz,
    qx: ell.qx,
    sign,
    lat0: sign * @math.atan(tphi0) * @geomath.deg_per_rad,
    k0,
    n0,
    m02,
    nrho0,
    k2: sq(k0),
    txi0,
    scxi0,
    sxi0,
  }
}

///|
fn AlbersEqualArea::ell(self : AlbersEqualArea) -> AlbersEll {
  {
    f: self.f,
    e2: self.e2,
    e: self.e,
    e2m: self.e2m,
    qz: self.qz,
    qx: self.qx,
  }
}

///|
/// Latitude of the origin (degrees), where the scale is `central_scale()`
/// in both directions.
pub fn AlbersEqualArea::origin_latitude(self : AlbersEqualArea) -> Double {
  self.lat0
}

///|
/// Scale factor at the origin latitude.
pub fn AlbersEqualArea::central_scale(self : AlbersEqualArea) -> Double {
  self.k0
}

///|
/// Equatorial radius (m).
pub fn AlbersEqualArea::equatorial_radius(self : AlbersEqualArea) -> Double {
  self.a
}

///|
/// Flattening.
pub fn AlbersEqualArea::flattening(self : AlbersEqualArea) -> Double {
  self.f
}

///|
/// Project (`lat`, `lon`) (degrees) with central meridian `lon0`. `x` and
/// `y` are relative to (origin_latitude(), lon0). `k` is the scale along
/// the parallel; the scale along the meridian is 1 / k, so areas are kept.
pub fn AlbersEqualArea::forward(
  self : AlbersEqualArea,
  lon0 : Double,
  lat : Double,
  lon : Double,
) -> ConicProjected {
  let ell = self.ell()
  let (dlon, _) = @geomath.ang_diff(lon0, lon)
  let (sphi, cphi0) = @geomath.sincosd(@geomath.lat_fix(lat) * self.sign)
  let cphi = @geomath.fmax(epsx, cphi0)
  let lam = dlon * @geomath.rad_per_deg
  let tphi = sphi / cphi
  let txi = ell.txif(tphi)
  let sxi = txi / hyp(txi)
  let dq = self.qz * dsn(txi, self.txi0, sxi, self.sxi0) * (txi - self.txi0)
  let drho = -self.a *
    dq /
    ((self.m02 - self.n0 * dq).sqrt() + self.nrho0 / self.a)
  let theta = self.k2 * self.n0 * lam
  let stheta = @math.sin(theta)
  let ctheta = @math.cos(theta)
  let t = self.nrho0 + self.n0 * drho
  let x = t *
    (if self.n0 != 0.0 { stheta / self.n0 } else { self.k2 * lam }) /
    self.k0
  let y = (
      self.nrho0 *
      (if self.n0 != 0.0 {
        (if ctheta < 0.0 { 1.0 - ctheta } else { sq(stheta) / (1.0 + ctheta) }) /
        self.n0
      } else {
        0.0
      }) -
      drho * ctheta
    ) /
    self.k0
  let k = self.k0 *
    (if t != 0.0 { t * hyp(self.fm * tphi) / self.a } else { 1.0 })
  { x, y: y * self.sign, gamma: self.sign * theta * @geomath.deg_per_rad, k, }
}

///|
/// Inverse of `forward`: geographic coordinates of (`x`, `y`) with central
/// meridian `lon0`.
pub fn AlbersEqualArea::reverse(
  self : AlbersEqualArea,
  lon0 : Double,
  x : Double,
  y0 : Double,
) -> ConicGeographic {
  let ell = self.ell()
  let y = y0 * self.sign
  let nx = self.k0 * self.n0 * x
  let ny = self.k0 * self.n0 * y
  let y1 = self.nrho0 - ny
  let den = @geomath.hypot(nx, y1) + self.nrho0
  let drho = if den != 0.0 {
    (self.k0 * x * nx - 2.0 * self.k0 * y * self.nrho0 + self.k0 * y * ny) / den
  } else {
    0.0
  }
  let dsxia = -self.scxi0 *
    (2.0 * self.nrho0 + self.n0 * drho) *
    drho /
    (sq(self.a) * self.qz)
  let txi = (self.txi0 + dsxia) /
    @geomath.fmax(1.0 - dsxia * (2.0 * self.txi0 + dsxia), epsx * epsx).sqrt()
  let tphi = ell.tphif(txi)
  let theta = @math.atan2(nx, y1)
  let lam = if self.n0 != 0.0 {
    theta / (self.k2 * self.n0)
  } else {
    x / (y1 * self.k0)
  }
  let lat = @geomath.atan2d(self.sign * tphi, 1.0)
  let lon = @geomath.ang_normalize(
    lam * @geomath.deg_per_rad + @geomath.ang_normalize(lon0),
  )
  let k = self.k0 *
    (if den != 0.0 {
      (self.nrho0 + self.n0 * drho) * hyp(self.fm * tphi) / self.a
    } else {
      1.0
    })
  let gamma = self.sign * theta * @geomath.deg_per_rad
  // Near a pole the meridian scale of an equal-area projection tends to 0:
  // on most conic grids, points within about 1e-6 degrees of the pole fall on
  // the pole's own map coordinates to rounding, so latitude cannot be
  // recovered more finely there. When the input is, to rounding, the image
  // of the pole, return the pole itself; a pole then round-trips exactly.
  if lat.abs() > 89.99 {
    let pole = if lat > 0.0 { 90.0 } else { -90.0 }
    let p = self.forward(lon0, pole, lon)
    let size = @geomath.fmax(
      @geomath.fmax(x.abs(), y0.abs()),
      self.nrho0 / self.k0,
    )
    if @geomath.hypot(p.x - x, p.y - y0) <= 8.0 * eps * size {
      return { lat: pole, lon, gamma, k: p.k, }
    }
  }
  { lat, lon, gamma, k, }
}
