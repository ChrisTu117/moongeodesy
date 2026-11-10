// Lambert conformal conic projection on an ellipsoid, a port of
// GeographicLib::LambertConformalConic (C. F. F. Karney, MIT/X11 License,
// https://geographiclib.sourceforge.io/). Snyder's formulas are rewritten
// with divided differences so that the projection stays accurate when the
// standard parallels are close together, near the equator (where it tends
// to Mercator), and near a pole (polar stereographic).

///|
/// A Lambert conformal conic projection. Coordinates are relative to the
/// origin at latitude `origin_latitude()` and the central meridian passed
/// to `forward` / `reverse`; `ConicGrid` adds a false origin and false
/// easting/northing as in EPSG.
pub struct LambertConformalConic {
  priv a : Double
  priv f : Double
  priv fm : Double
  priv es : Double
  priv sign : Double
  priv n : Double
  priv nc : Double
  priv t0nm1 : Double
  priv scale : Double
  priv lat0 : Double
  priv k0 : Double
  priv scbet0 : Double
  priv tchi0 : Double
  priv scchi0 : Double
  priv psi0 : Double
  priv nrho0 : Double
  priv drhomax : Double
}

///|
/// log of the largest finite double, plus 2 (GeographicLib's ahypover_).
let ahypover : Double = 53.0 * @math.ln(2.0) + 2.0

///|
/// A Lambert conformal conic projection with two standard parallels
/// `stdlat1` and `stdlat2` (degrees), on which the scale is `k1`. `a` is the
/// equatorial radius (m) and `f` the flattening.
pub fn LambertConformalConic::new(
  a : Double,
  f : Double,
  stdlat1 : Double,
  stdlat2 : Double,
  k1~ : Double = 1.0,
) -> LambertConformalConic raise ConicError {
  check_ellipsoid(a, f)
  check_scale(k1)
  check_stdlat(stdlat1)
  check_stdlat(stdlat2)
  let (s1, c1) = @geomath.sincosd(stdlat1)
  let (s2, c2) = @geomath.sincosd(stdlat2)
  if (c1 == 0.0 || c2 == 0.0) && !(c1 == c2 && s1 == s2) {
    raise IncompatibleStandardLatitudes(stdlat1, stdlat2)
  }
  lcc_init(a, f, s1, c1, s2, c2, k1)
}

///|
/// A Lambert conformal conic projection with one standard parallel
/// `stdlat` (degrees), on which the scale is `k0`. This is EPSG's
/// "Lambert Conic Conformal (1SP)" with the natural origin on `stdlat`.
pub fn LambertConformalConic::new_1sp(
  a : Double,
  f : Double,
  stdlat : Double,
  k0~ : Double = 1.0,
) -> LambertConformalConic raise ConicError {
  check_ellipsoid(a, f)
  check_scale(k0)
  check_stdlat(stdlat)
  let (s, c) = @geomath.sincosd(stdlat)
  lcc_init(a, f, s, c, s, c, k0)
}

///|
fn lcc_deatanhe(x : Double, y : Double, e2 : Double, es : Double) -> Double {
  let t = x - y
  let d = 1.0 - e2 * x * y
  if t != 0.0 {
    @geomath.eatanhe(t / d, es) / t
  } else {
    e2 / d
  }
}

///|
fn lcc_init(
  a : Double,
  f : Double,
  sinlat1 : Double,
  coslat1 : Double,
  sinlat2 : Double,
  coslat2 : Double,
  k1 : Double,
) -> LambertConformalConic {
  let fm = 1.0 - f
  let e2 = f * (2.0 - f)
  let es = (if f < 0.0 { -1.0 } else { 1.0 }) * e2.abs().sqrt()
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
  let tbet1 = fm * tphi1
  let scbet1 = hyp(tbet1)
  let tbet2 = fm * tphi2
  let scbet2 = hyp(tbet2)
  let scphi1 = 1.0 / cphi1
  let xi1 = @geomath.eatanhe(sphi1, es)
  let shxi1 = @math.sinh(xi1)
  let chxi1 = hyp(shxi1)
  let tchi1 = chxi1 * tphi1 - shxi1 * scphi1
  let scchi1 = hyp(tchi1)
  let scphi2 = 1.0 / cphi2
  let xi2 = @geomath.eatanhe(sphi2, es)
  let shxi2 = @math.sinh(xi2)
  let chxi2 = hyp(shxi2)
  let tchi2 = chxi2 * tphi2 - shxi2 * scphi2
  let scchi2 = hyp(tchi2)
  let psi1 = @math.asinh(tchi1)
  let mut n = 0.0
  let mut nc = 0.0
  let tphi0 = if tphi2 - tphi1 != 0.0 {
    let num = dlog1p(sq(tbet2) / (1.0 + scbet2), sq(tbet1) / (1.0 + scbet1)) *
      dhyp(tbet2, tbet1, scbet2, scbet1) *
      fm
    let den = dasinh(tphi2, tphi1, scphi2, scphi1) -
      lcc_deatanhe(sphi2, sphi1, e2, es) * dsn(tphi2, tphi1, sphi2, sphi1)
    n = num / den
    if n < 0.25 {
      nc = ((1.0 - n) * (1.0 + n)).sqrt()
    } else {
      // Compute nc = cos(phi0) accurately when n is close to 1.
      let s1 = tphi1 * (2.0 * shxi1 * chxi1 * scphi1 - e2 * tphi1) -
        sq(shxi1) * (1.0 + 2.0 * sq(tphi1))
      let s2 = tphi2 * (2.0 * shxi2 * chxi2 * scphi2 - e2 * tphi2) -
        sq(shxi2) * (1.0 + 2.0 * sq(tphi2))
      let t1 = if tchi1 < 0.0 {
        scbet1 - tchi1
      } else {
        (s1 + 1.0) / (scbet1 + tchi1)
      }
      let t2 = if tchi2 < 0.0 {
        scbet2 - tchi2
      } else {
        (s2 + 1.0) / (scbet2 + tchi2)
      }
      let a2 = -(s2 / (scbet2 + scchi2) + t2) / (2.0 * scbet2)
      let a1 = -(s1 / (scbet1 + scchi1) + t1) / (2.0 * scbet1)
      let mut t = dlog1p(a2, a1) / den
      t = t *
        (
          (
            (if tchi2 >= 0.0 { scchi2 + tchi2 } else { 1.0 / (scchi2 - tchi2) }) +
            (if tchi1 >= 0.0 { scchi1 + tchi1 } else { 1.0 / (scchi1 - tchi1) })
          ) /
          (4.0 * scbet1 * scbet2)
        ) *
        fm
      let tbm = (
          (if tbet1 > 0.0 { 1.0 / (scbet1 + tbet1) } else { scbet1 - tbet1 }) +
          (if tbet2 > 0.0 { 1.0 / (scbet2 + tbet2) } else { scbet2 - tbet2 })
        ) /
        (scbet1 + scbet2)
      let dtchi = den / dasinh(tchi2, tchi1, scchi2, scchi1)
      let dbet = e2 /
        fm *
        (1.0 / (scbet2 + fm * scphi2) + 1.0 / (scbet1 + fm * scphi1))
      let xiz = @geomath.eatanhe(1.0, es)
      let shxiz = @math.sinh(xiz)
      let chxiz = hyp(shxiz)
      let dxiz1 = lcc_deatanhe(1.0, sphi1, e2, es) / (scphi1 * (tphi1 + scphi1))
      let dxiz2 = lcc_deatanhe(1.0, sphi2, e2, es) / (scphi2 * (tphi2 + scphi2))
      let dshxiz1 = dsinh(xiz, xi1, shxiz, shxi1, chxiz, chxi1) * dxiz1
      let dshxiz2 = dsinh(xiz, xi2, shxiz, shxi2, chxiz, chxi2) * dxiz2
      let dchxiz1 = dhyp(shxiz, shxi1, chxiz, chxi1) * dshxiz1
      let dchxiz2 = dhyp(shxiz, shxi2, chxiz, chxi2) * dshxiz2
      let amu12 = -scphi1 * dchxiz1 +
        tphi1 * dshxiz1 -
        scphi2 * dchxiz2 +
        tphi2 * dshxiz2
      let dxi = lcc_deatanhe(sphi1, sphi2, e2, es) *
        dsn(tphi2, tphi1, sphi2, sphi1)
      let dnu12 = (if f * 4.0 * scphi2 * dshxiz2 > f * scphi1 * dshxiz1 {
          (dshxiz1 + dshxiz2) / 2.0 * dhyp(tphi1, tphi2, scphi1, scphi2) -
          (scphi1 + scphi2) /
          2.0 *
          dsinh(xi1, xi2, shxi1, shxi2, chxi1, chxi2) *
          dxi
        } else {
          (scphi2 * dshxiz2 - scphi1 * dshxiz1) / (tphi2 - tphi1)
        }) +
        (tphi1 + tphi2) /
        2.0 *
        dhyp(shxi1, shxi2, chxi1, chxi2) *
        dsinh(xi1, xi2, shxi1, shxi2, chxi1, chxi2) *
        dxi -
        (dchxiz1 + dchxiz2) / 2.0
      let dchia = amu12 - dnu12 * (scphi2 + scphi1)
      let tam = (dchia - dtchi * dbet) / (scchi1 + scchi2)
      t = t * (tbm - tam)
      nc = (@geomath.fmax(0.0, t) * (1.0 + n)).sqrt()
    }
    let r = @geomath.hypot(n, nc)
    n = n / r
    nc = nc / r
    n / nc
  } else {
    nc = 1.0 / hyp(tphi1)
    n = tphi1 * nc
    if polar {
      nc = 0.0
    }
    tphi1
  }
  let scbet0 = hyp(fm * tphi0)
  let shxi0 = @math.sinh(@geomath.eatanhe(n, es))
  let tchi0 = tphi0 * hyp(shxi0) - shxi0 * hyp(tphi0)
  let scchi0 = hyp(tchi0)
  let psi0 = @math.asinh(tchi0)
  let lat0 = @math.atan(sign * tphi0) * @geomath.deg_per_rad
  let t0nm1 = @math.expm1(-n * psi0)
  let ex1 = if tchi1 >= 0.0 { scchi1 + tchi1 } else { 1.0 / (scchi1 - tchi1) }
  let scale = a * k1 / scbet1 * @math.exp(-(sq(nc) / (1.0 + n)) * psi1) * ex1
  let k0 = k1 *
    (scbet0 / scbet1) *
    @math.exp(
      -(sq(nc) / (1.0 + n)) *
      dasinh(tchi1, tchi0, scchi1, scchi0) *
      (tchi1 - tchi0),
    ) *
    ex1 /
    (scchi0 + tchi0)
  let nrho0 = if polar { 0.0 } else { a * k0 / scbet0 }
  // drho at latitude -90, the limit for Reverse.
  let drhomax = {
    let sphi = -1.0
    let cphi = epsx
    let tphi = sphi / cphi
    let scphi = 1.0 / cphi
    let shxi = @math.sinh(@geomath.eatanhe(sphi, es))
    let tchi = hyp(shxi) * tphi - shxi * scphi
    let scchi = hyp(tchi)
    let psi = @math.asinh(tchi)
    let dpsi = dasinh(tchi, tchi0, scchi, scchi0) * (tchi - tchi0)
    -scale *
    (if 2.0 * nc < 1.0 && dpsi != 0.0 {
      (
        @math.exp(sq(nc) / (1.0 + n) * psi) *
        (if tchi > 0.0 { 1.0 / (scchi + tchi) } else { scchi - tchi }) -
        (t0nm1 + 1.0)
      ) /
      -n
    } else {
      dexp(-n * psi, -n * psi0) * dpsi
    })
  }
  {
    a,
    f,
    fm,
    es,
    sign,
    n,
    nc,
    t0nm1,
    scale,
    lat0,
    k0,
    scbet0,
    tchi0,
    scchi0,
    psi0,
    nrho0,
    drhomax,
  }
}

///|
/// Latitude of the origin (degrees): where the scale is smallest (and
/// equal to `central_scale()`), sin(lat0) = n, the cone constant.
pub fn LambertConformalConic::origin_latitude(
  self : LambertConformalConic,
) -> Double {
  self.lat0
}

///|
/// Scale factor at the origin latitude.
pub fn LambertConformalConic::central_scale(
  self : LambertConformalConic,
) -> Double {
  self.k0
}

///|
/// The cone constant n = sin(origin latitude), signed.
pub fn LambertConformalConic::cone_constant(
  self : LambertConformalConic,
) -> Double {
  self.sign * self.n
}

///|
/// Equatorial radius (m).
pub fn LambertConformalConic::equatorial_radius(
  self : LambertConformalConic,
) -> Double {
  self.a
}

///|
/// Flattening.
pub fn LambertConformalConic::flattening(
  self : LambertConformalConic,
) -> Double {
  self.f
}

///|
/// Project (`lat`, `lon`) (degrees) with central meridian `lon0`. `x` and
/// `y` are relative to (origin_latitude(), lon0). A latitude outside
/// [-90, 90] gives NaN; the pole opposite the cone's apex maps to infinity.
pub fn LambertConformalConic::forward(
  self : LambertConformalConic,
  lon0 : Double,
  lat : Double,
  lon : Double,
) -> ConicProjected {
  let (dlon, _) = @geomath.ang_diff(lon0, lon)
  let (sphi, cphi0) = @geomath.sincosd(@geomath.lat_fix(lat) * self.sign)
  let cphi = @geomath.fmax(epsx, cphi0)
  let n = self.n
  let nc = self.nc
  let lam = dlon * @geomath.rad_per_deg
  let tphi = sphi / cphi
  let scbet = hyp(self.fm * tphi)
  let scphi = 1.0 / cphi
  let shxi = @math.sinh(@geomath.eatanhe(sphi, self.es))
  let tchi = hyp(shxi) * tphi - shxi * scphi
  let scchi = hyp(tchi)
  let psi = @math.asinh(tchi)
  let theta = n * lam
  let stheta = @math.sin(theta)
  let ctheta = @math.cos(theta)
  let dpsi = dasinh(tchi, self.tchi0, scchi, self.scchi0) * (tchi - self.tchi0)
  let drho = -self.scale *
    (if 2.0 * nc < 1.0 && dpsi != 0.0 {
      (
        @math.exp(sq(nc) / (1.0 + n) * psi) *
        (if tchi > 0.0 { 1.0 / (scchi + tchi) } else { scchi - tchi }) -
        (self.t0nm1 + 1.0)
      ) /
      -n
    } else {
      dexp(-n * psi, -n * self.psi0) * dpsi
    })
  let x = (self.nrho0 + n * drho) * (if n != 0.0 { stheta / n } else { lam })
  let y = self.nrho0 *
    (if n != 0.0 {
      (if ctheta < 0.0 { 1.0 - ctheta } else { sq(stheta) / (1.0 + ctheta) }) /
      n
    } else {
      0.0
    }) -
    drho * ctheta
  let k = self.k0 *
    (scbet / self.scbet0) /
    (
      @math.exp(-(sq(nc) / (1.0 + n)) * dpsi) *
      (if tchi >= 0.0 { scchi + tchi } else { 1.0 / (scchi - tchi) }) /
      (self.scchi0 + self.tchi0)
    )
  { x, y: y * self.sign, gamma: self.sign * theta * @geomath.deg_per_rad, k, }
}

///|
/// Inverse of `forward`: geographic coordinates of (`x`, `y`) with central
/// meridian `lon0`.
pub fn LambertConformalConic::reverse(
  self : LambertConformalConic,
  lon0 : Double,
  x : Double,
  y0 : Double,
) -> ConicGeographic {
  let y = y0 * self.sign
  let n = self.n
  let nc = self.nc
  let nx = n * x
  let ny = if n != 0.0 { n * y } else { 0.0 }
  let y1 = self.nrho0 - ny
  let den = @geomath.hypot(nx, y1) + self.nrho0
  let mut drho = if den != 0.0 && @geomath.is_finite(den) {
    (x * nx + y * (ny - 2.0 * self.nrho0)) / den
  } else {
    den
  }
  drho = @geomath.fmin(drho, self.drhomax)
  if n == 0.0 {
    drho = @geomath.fmax(drho, -self.drhomax)
  }
  let tnm1 = self.t0nm1 + n * drho / self.scale
  let dpsi = if den == 0.0 {
    0.0
  } else if tnm1 + 1.0 != 0.0 {
    -dlog1p(tnm1, self.t0nm1) * drho / self.scale
  } else {
    ahypover
  }
  let tchi = if 2.0 * n <= 1.0 {
    let psi = self.psi0 + dpsi
    let tchia = @math.sinh(psi)
    let scchi = hyp(tchia)
    let dtchi = dsinh(psi, self.psi0, tchia, self.tchi0, scchi, self.scchi0) *
      dpsi
    self.tchi0 + dtchi
  } else {
    let tn = if tnm1 + 1.0 == 0.0 { epsx } else { tnm1 + 1.0 }
    let sh = @math.sinh(
      -sq(nc) /
      (n * (1.0 + n)) *
      (if 2.0 * tn > 1.0 { @math.ln_1p(tnm1) } else { @math.ln(tn) }),
    )
    sh * (tn + 1.0 / tn) / 2.0 - hyp(sh) * (tnm1 * (tn + 1.0) / tn) / 2.0
  }
  let gamma = @math.atan2(nx, y1)
  let tphi = @geomath.tauf(tchi, self.es)
  let scbet = hyp(self.fm * tphi)
  let scchi = hyp(tchi)
  let lam = if n != 0.0 { gamma / n } else { x / y1 }
  let lat = @geomath.atan2d(self.sign * tphi, 1.0)
  let lon = @geomath.ang_normalize(
    lam * @geomath.deg_per_rad + @geomath.ang_normalize(lon0),
  )
  let k = self.k0 *
    (scbet / self.scbet0) /
    (
      @math.exp(if nc != 0.0 { -(sq(nc) / (1.0 + n)) * dpsi } else { 0.0 }) *
      (if tchi >= 0.0 { scchi + tchi } else { 1.0 / (scchi - tchi) }) /
      (self.scchi0 + self.tchi0)
    )
  { lat, lon, gamma: gamma / self.sign * @geomath.deg_per_rad, k, }
}
