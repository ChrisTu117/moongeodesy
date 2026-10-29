// MoonBit port of GeographicLib::EllipticFunction (Charles Karney, MIT/X11
// License). Carlson symmetric integrals and the incomplete integrals used by
// GeodesicExact. Transcription of EllipticFunction.cpp from GeographicLib 2.5.
// Jacobi amplitude / sncndn are not used by the geodesic and are omitted.

///|
let tol_rf : Double = @math.pow(3.0 * @geomath.dbl_epsilon * 0.01, 1.0 / 8.0)

///|
let tol_rd : Double = @math.pow(0.2 * (@geomath.dbl_epsilon * 0.01), 1.0 / 8.0)

///|
let tol_rg0 : Double = 2.7 * (@geomath.dbl_epsilon * 0.01).sqrt()

///|
let tol_jac : Double = (@geomath.dbl_epsilon * 0.01).sqrt()

///|
fn signbit(x : Double) -> Bool {
  x.reinterpret_as_int64() < 0
}

///|
priv struct Elliptic {
  k2 : Double
  kp2 : Double
  alphap2 : Double
  eps : Double
  ec : Double
  dc : Double
  hc : Double
}

///|
fn elliptic_rc(x : Double, y : Double) -> Double {
  if !(x >= y) {
    @math.atan(((y - x) / x).sqrt()) / (y - x).sqrt()
  } else if x == y {
    1.0 / y.sqrt()
  } else {
    let arg = if y > 0.0 { ((x - y) / y).sqrt() } else { (-x / y).sqrt() }
    @math.asinh(arg) / (x - y).sqrt()
  }
}

///|
fn elliptic_rf3(x : Double, y : Double, z : Double) -> Double {
  let a0 = (x + y + z) / 3.0
  let mut an = a0
  let q = @geomath.fmax(
      @geomath.fmax((a0 - x).abs(), (a0 - y).abs()),
      (a0 - z).abs(),
    ) /
    tol_rf
  let mut x0 = x
  let mut y0 = y
  let mut z0 = z
  let mut mul = 1.0
  while q >= mul * an.abs() {
    let lam = x0.sqrt() * y0.sqrt() +
      y0.sqrt() * z0.sqrt() +
      z0.sqrt() * x0.sqrt()
    an = (an + lam) / 4.0
    x0 = (x0 + lam) / 4.0
    y0 = (y0 + lam) / 4.0
    z0 = (z0 + lam) / 4.0
    mul = mul * 4.0
  }
  let xv = (a0 - x) / (mul * an)
  let yv = (a0 - y) / (mul * an)
  let zv = -(xv + yv)
  let e2 = xv * yv - zv * zv
  let e3 = xv * yv * zv
  (
    e3 * (6930.0 * e3 + e2 * (15015.0 * e2 - 16380.0) + 17160.0) +
    e2 * ((10010.0 - 5775.0 * e2) * e2 - 24024.0) +
    240240.0
  ) /
  (240240.0 * an.sqrt())
}

///|
fn elliptic_rf2(x : Double, y : Double) -> Double {
  let mut xn = x.sqrt()
  let mut yn = y.sqrt()
  if xn < yn {
    let t = xn
    xn = yn
    yn = t
  }
  while (xn - yn).abs() > tol_rg0 * xn {
    let t = (xn + yn) / 2.0
    yn = (xn * yn).sqrt()
    xn = t
  }
  @geomath.pi / (xn + yn)
}

///|
fn elliptic_rg2(x : Double, y : Double) -> Double {
  let x0 = @geomath.fmax(x, y).sqrt()
  let y0 = @geomath.fmin(x, y).sqrt()
  let mut xn = x0
  let mut yn = y0
  let mut s = 0.0
  let mut mul = 0.25
  while (xn - yn).abs() > tol_rg0 * xn {
    let t = (xn + yn) / 2.0
    yn = (xn * yn).sqrt()
    xn = t
    mul = mul * 2.0
    let d = xn - yn
    s = s + mul * d * d
  }
  (@geomath.sq((x0 + y0) / 2.0) - s) * @geomath.pi / (2.0 * (xn + yn))
}

///|
fn elliptic_rd(x : Double, y : Double, z : Double) -> Double {
  let a0 = (x + y + 3.0 * z) / 5.0
  let mut an = a0
  let q = @geomath.fmax(
      @geomath.fmax((a0 - x).abs(), (a0 - y).abs()),
      (a0 - z).abs(),
    ) /
    tol_rd
  let mut x0 = x
  let mut y0 = y
  let mut z0 = z
  let mut mul = 1.0
  let mut s = 0.0
  while q >= mul * an.abs() {
    let lam = x0.sqrt() * y0.sqrt() +
      y0.sqrt() * z0.sqrt() +
      z0.sqrt() * x0.sqrt()
    s = s + 1.0 / (mul * z0.sqrt() * (z0 + lam))
    an = (an + lam) / 4.0
    x0 = (x0 + lam) / 4.0
    y0 = (y0 + lam) / 4.0
    z0 = (z0 + lam) / 4.0
    mul = mul * 4.0
  }
  let xv = (a0 - x) / (mul * an)
  let yv = (a0 - y) / (mul * an)
  let zv = -(xv + yv) / 3.0
  let e2 = xv * yv - 6.0 * zv * zv
  let e3 = (3.0 * xv * yv - 8.0 * zv * zv) * zv
  let e4 = 3.0 * (xv * yv - zv * zv) * zv * zv
  let e5 = xv * yv * zv * zv * zv
  (
    (471240.0 - 540540.0 * e2) * e5 +
    (612612.0 * e2 - 540540.0 * e3 - 556920.0) * e4 +
    e3 * (306306.0 * e3 + e2 * (675675.0 * e2 - 706860.0) + 680680.0) +
    e2 * ((417690.0 - 255255.0 * e2) * e2 - 875160.0) +
    4084080.0
  ) /
  (4084080.0 * mul * an * an.sqrt()) +
  3.0 * s
}

///|
fn elliptic_rj(x : Double, y : Double, z : Double, p : Double) -> Double {
  let a0 = (x + y + z + 2.0 * p) / 5.0
  let mut an = a0
  let delta = (p - x) * (p - y) * (p - z)
  let q = @geomath.fmax(
      @geomath.fmax((a0 - x).abs(), (a0 - y).abs()),
      @geomath.fmax((a0 - z).abs(), (a0 - p).abs()),
    ) /
    tol_rd
  let mut x0 = x
  let mut y0 = y
  let mut z0 = z
  let mut p0 = p
  let mut mul = 1.0
  let mut mul3 = 1.0
  let mut s = 0.0
  while q >= mul * an.abs() {
    let lam = x0.sqrt() * y0.sqrt() +
      y0.sqrt() * z0.sqrt() +
      z0.sqrt() * x0.sqrt()
    let d0 = (p0.sqrt() + x0.sqrt()) *
      (p0.sqrt() + y0.sqrt()) *
      (p0.sqrt() + z0.sqrt())
    let e0 = delta / (mul3 * @geomath.sq(d0))
    s = s + elliptic_rc(1.0, 1.0 + e0) / (mul * d0)
    an = (an + lam) / 4.0
    x0 = (x0 + lam) / 4.0
    y0 = (y0 + lam) / 4.0
    z0 = (z0 + lam) / 4.0
    p0 = (p0 + lam) / 4.0
    mul = mul * 4.0
    mul3 = mul3 * 64.0
  }
  let xv = (a0 - x) / (mul * an)
  let yv = (a0 - y) / (mul * an)
  let zv = (a0 - z) / (mul * an)
  let pv = -(xv + yv + zv) / 2.0
  let e2 = xv * yv + xv * zv + yv * zv - 3.0 * pv * pv
  let e3 = xv * yv * zv + 2.0 * pv * (e2 + 2.0 * pv * pv)
  let e4 = (2.0 * xv * yv * zv + pv * (e2 + 3.0 * pv * pv)) * pv
  let e5 = xv * yv * zv * pv * pv
  (
    (471240.0 - 540540.0 * e2) * e5 +
    (612612.0 * e2 - 540540.0 * e3 - 556920.0) * e4 +
    e3 * (306306.0 * e3 + e2 * (675675.0 * e2 - 706860.0) + 680680.0) +
    e2 * ((417690.0 - 255255.0 * e2) * e2 - 875160.0) +
    4084080.0
  ) /
  (4084080.0 * mul * an * an.sqrt()) +
  6.0 * s
}

///|
fn Elliptic::from_k2(k2 : Double) -> Elliptic {
  Elliptic::from_parts(k2, 0.0, 1.0 - k2, 1.0)
}

///|
fn Elliptic::from_parts(
  k2 : Double,
  alpha2 : Double,
  kp2 : Double,
  alphap2 : Double,
) -> Elliptic {
  let eps = k2 / @geomath.sq(kp2.sqrt() + 1.0)
  let (kc, ec, dc) = if k2 != 0.0 {
    if kp2 != 0.0 {
      (
        elliptic_rf2(kp2, 1.0),
        2.0 * elliptic_rg2(kp2, 1.0),
        elliptic_rd(0.0, kp2, 1.0) / 3.0,
      )
    } else {
      (@double.infinity, 1.0, @double.infinity)
    }
  } else {
    let half_pi = @geomath.pi / 2.0
    (half_pi, half_pi, half_pi / 2.0)
  }
  let hc = if alpha2 != 0.0 {
    let rj = if kp2 != 0.0 && alphap2 != 0.0 {
      elliptic_rj(0.0, kp2, 1.0, alphap2)
    } else {
      @double.infinity
    }
    let rc = if kp2 != 0.0 {
      0.0
    } else if alphap2 != 0.0 {
      elliptic_rc(1.0, alphap2)
    } else {
      @double.infinity
    }
    if kp2 != 0.0 {
      kc - (if alphap2 != 0.0 { alphap2 * rj } else { 0.0 }) / 3.0
    } else {
      rc
    }
  } else if kp2 == 1.0 {
    @geomath.pi / 4.0
  } else if kp2 == 0.0 {
    1.0
  } else {
    kp2 * elliptic_rd(0.0, 1.0, kp2) / 3.0
  }
  { k2, kp2, alphap2, eps, ec, dc, hc, }
}

///|
fn Elliptic::delta_amp(self : Elliptic, sn : Double, cn : Double) -> Double {
  if self.k2 < 0.0 {
    (1.0 - self.k2 * sn * sn).sqrt()
  } else {
    (self.kp2 + self.k2 * cn * cn).sqrt()
  }
}

///|
fn Elliptic::e_inc(
  self : Elliptic,
  sn : Double,
  cn : Double,
  dn : Double,
) -> Double {
  let cn2 = cn * cn
  let dn2 = dn * dn
  let sn2 = sn * sn
  let ei = if cn2 != 0.0 {
    sn.abs() *
    (if self.k2 <= 0.0 {
      elliptic_rf3(cn2, dn2, 1.0) -
      self.k2 * sn2 * elliptic_rd(cn2, dn2, 1.0) / 3.0
    } else if self.kp2 >= 0.0 {
      self.kp2 * elliptic_rf3(cn2, dn2, 1.0) +
      self.k2 * self.kp2 * sn2 * elliptic_rd(cn2, 1.0, dn2) / 3.0 +
      self.k2 * cn.abs() / dn
    } else {
      -self.kp2 * sn2 * elliptic_rd(dn2, 1.0, cn2) / 3.0 + dn / cn.abs()
    })
  } else {
    self.ec
  }
  let ei2 = if signbit(cn) { 2.0 * self.ec - ei } else { ei }
  @geomath.copysign(ei2, sn)
}

///|
fn Elliptic::d_inc(
  self : Elliptic,
  sn : Double,
  cn : Double,
  dn : Double,
) -> Double {
  let cn2 = cn * cn
  let dn2 = dn * dn
  let sn2 = sn * sn
  let di = if cn2 != 0.0 {
    sn.abs() * sn2 * elliptic_rd(cn2, dn2, 1.0) / 3.0
  } else {
    self.dc
  }
  let di2 = if signbit(cn) { 2.0 * self.dc - di } else { di }
  @geomath.copysign(di2, sn)
}

///|
fn Elliptic::h_inc(
  self : Elliptic,
  sn : Double,
  cn : Double,
  dn : Double,
) -> Double {
  let cn2 = cn * cn
  let dn2 = dn * dn
  let sn2 = sn * sn
  let hi = if cn2 != 0.0 {
    sn.abs() *
    (
      elliptic_rf3(cn2, dn2, 1.0) -
      self.alphap2 *
      sn2 *
      elliptic_rj(cn2, dn2, 1.0, cn2 + self.alphap2 * sn2) /
      3.0
    )
  } else {
    self.hc
  }
  let hi2 = if signbit(cn) { 2.0 * self.hc - hi } else { hi }
  @geomath.copysign(hi2, sn)
}

///|
fn Elliptic::delta_e(
  self : Elliptic,
  sn0 : Double,
  cn0 : Double,
  dn : Double,
) -> Double {
  let (sn, cn) = if signbit(cn0) { (-sn0, -cn0) } else { (sn0, cn0) }
  self.e_inc(sn, cn, dn) * (@geomath.pi / 2.0) / self.ec - @math.atan2(sn, cn)
}

///|
fn Elliptic::delta_d(
  self : Elliptic,
  sn0 : Double,
  cn0 : Double,
  dn : Double,
) -> Double {
  let (sn, cn) = if signbit(cn0) { (-sn0, -cn0) } else { (sn0, cn0) }
  self.d_inc(sn, cn, dn) * (@geomath.pi / 2.0) / self.dc - @math.atan2(sn, cn)
}

///|
fn Elliptic::delta_h(
  self : Elliptic,
  sn0 : Double,
  cn0 : Double,
  dn : Double,
) -> Double {
  let (sn, cn) = if signbit(cn0) { (-sn0, -cn0) } else { (sn0, cn0) }
  self.h_inc(sn, cn, dn) * (@geomath.pi / 2.0) / self.hc - @math.atan2(sn, cn)
}

///|
fn Elliptic::einv(self : Elliptic, x0 : Double) -> Double {
  let mut x = x0
  let n = (x / (2.0 * self.ec) + 0.5).floor()
  x = x - 2.0 * self.ec * n
  let mut phi = @geomath.pi * x / (2.0 * self.ec)
  phi = phi - self.eps * @math.sin(2.0 * phi) / 2.0
  let mut ok = false
  for _ in 0..<25 {
    let sn = @math.sin(phi)
    let cn = @math.cos(phi)
    let dn = self.delta_amp(sn, cn)
    let err = (self.e_inc(sn, cn, dn) - x) / dn
    phi = phi - err
    if !(err.abs() > tol_jac) {
      ok = true
      break
    }
  }
  if !ok {
    abort("Convergence failure in EllipticFunction")
  }
  n * @geomath.pi + phi
}

///|
fn Elliptic::delta_einv(
  self : Elliptic,
  stau0 : Double,
  ctau0 : Double,
) -> Double {
  let (stau, ctau) = if signbit(ctau0) {
    (-stau0, -ctau0)
  } else {
    (stau0, ctau0)
  }
  let tau = @math.atan2(stau, ctau)
  self.einv(tau * self.ec / (@geomath.pi / 2.0)) - tau
}
