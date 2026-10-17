// Auxiliary latitudes and divided differences used by the rhumb-line port.
// Transcribed from GeographicLib AuxAngle / AuxLatitude / DAuxLatitude
// (Charles Karney, MIT/X11 License), series path only (order 6).
// C. F. F. Karney, "On auxiliary latitudes", Survey Review 56(395), 165-180
// (2024), https://doi.org/10.1080/00396265.2023.2217604

///|
struct AuxAngle {
  y : Double
  x : Double
}

///|
fn AuxAngle::from_degrees(d : Double) -> AuxAngle {
  let (s, c) = @geomath.sincosd(d)
  { y: s, x: c, }
}

///|
fn AuxAngle::from_radians(r : Double) -> AuxAngle {
  { y: @math.sin(r), x: @math.cos(r), }
}

///|
fn AuxAngle::to_degrees(self : AuxAngle) -> Double {
  @geomath.atan2d(self.y, self.x)
}

///|
fn AuxAngle::to_radians(self : AuxAngle) -> Double {
  @math.atan2(self.y, self.x)
}

///|
fn AuxAngle::tan(self : AuxAngle) -> Double {
  // cos(±90°) is a signed zero. GeographicLib's sincosd adds 0 to force
  // +0; MoonBit folds that add and keeps -0, which flips tan(90°).
  // A zero cosine is an infinite tangent with the sign of the sine.
  if self.x == 0.0 {
    if self.y == 0.0 {
      @double.not_a_number
    } else {
      @geomath.copysign(@double.infinity, self.y)
    }
  } else {
    self.y / self.x
  }
}

///|
fn AuxAngle::lam(self : AuxAngle) -> Double {
  @math.asinh(self.tan())
}

///|
fn AuxAngle::normalized(self : AuxAngle) -> AuxAngle {
  let too_big = self.y.abs() > @double.max_value / 2.0 &&
    self.x.abs() > @double.max_value / 2.0
  if self.tan().is_nan() || too_big {
    return { y: @double.not_a_number, x: @double.not_a_number, }
  }
  let r = @geomath.hypot(self.y, self.x)
  let mut y = self.y / r
  let mut x = self.x / r
  if y.is_nan() {
    y = @geomath.copysign(1.0, self.y)
  }
  if x.is_nan() {
    x = @geomath.copysign(1.0, self.x)
  }
  { y, x, }
}

///|
fn AuxAngle::add(self : AuxAngle, p : AuxAngle) -> AuxAngle {
  if p.tan() != 0.0 {
    { y: self.y * p.x + self.x * p.y, x: self.x * p.x - self.y * p.y, }
  } else {
    self
  }
}

///|
fn aux_sc(tphi : Double) -> Double {
  @geomath.hypot(1.0, tphi)
}

///|
fn aux_sn(tphi : Double) -> Double {
  if tphi.is_inf() {
    @geomath.copysign(1.0, tphi)
  } else {
    tphi / aux_sc(tphi)
  }
}

///|
fn aux_h(tphi : Double) -> Double {
  tphi * aux_sn(tphi) / 2.0
}

///|
fn aux_datan(x : Double, y : Double) -> Double {
  let d = y - x
  let xy = x * y
  if x == y {
    1.0 / (1.0 + xy)
  } else if xy.is_inf() && xy > 0.0 {
    0.0
  } else if 2.0 * xy > -1.0 {
    @math.atan(d / (1.0 + xy)) / d
  } else {
    (@math.atan(y) - @math.atan(x)) / d
  }
}

///|
fn aux_dasinh(x : Double, y : Double) -> Double {
  let d = y - x
  let xy = x * y
  let hx = aux_sc(x)
  let hy = aux_sc(y)
  if x == y {
    1.0 / hx
  } else if d.is_inf() {
    0.0
  } else if xy > 0.0 {
    let inner = if xy < 1.0 {
      d * (x + y) / (x * hy + y * hx)
    } else {
      d * (1.0 / x + 1.0 / y) / (hy / y + hx / x)
    }
    @math.asinh(inner) / d
  } else {
    (@math.asinh(y) - @math.asinh(x)) / d
  }
}

///|
fn aux_dh(x : Double, y : Double) -> Double {
  if (x + y).is_nan() {
    return x + y
  }
  if x.is_inf() {
    return @geomath.copysign(0.5, x)
  }
  if y.is_inf() {
    return @geomath.copysign(0.5, y)
  }
  let sx = aux_sn(x)
  let sy = aux_sn(y)
  let d = sx * x + sy * y
  if d / 2.0 == 0.0 {
    return (x + y) / 2.0
  }
  if x * y <= 0.0 {
    return (aux_h(y) - aux_h(x)) / (y - x)
  }
  let scx = aux_sc(x)
  let scy = aux_sc(y)
  (x + y) /
  (2.0 * d) *
  (@geomath.sq(sx * sy) + @geomath.sq(sy / scx) + @geomath.sq(sx / scy))
}

///|
/// (psi2 - psi1) / (chi2 - chi1), with x = tan(chi1), y = tan(chi2).
fn aux_dlam(x : Double, y : Double) -> Double {
  if x == y {
    aux_sc(x)
  } else if x.is_nan() || y.is_nan() {
    @double.not_a_number
  } else if x.is_inf() || y.is_inf() {
    @double.infinity
  } else {
    aux_dasinh(x, y) / aux_datan(x, y)
  }
}

///|
/// Divided difference of log(sec chi) with respect to isometric latitude.
fn aux_dp0_dpsi(x : Double, y : Double) -> Double {
  if x == y {
    aux_sn(x)
  } else if (x + y).is_nan() {
    x + y
  } else if x.is_inf() {
    @geomath.copysign(1.0, x)
  } else if y.is_inf() {
    @geomath.copysign(1.0, y)
  } else {
    aux_dasinh(aux_h(x), aux_h(y)) * aux_dh(x, y) / aux_dasinh(x, y)
  }
}

///|
fn aux_clenshaw(
  sinp : Bool,
  szeta : Double,
  czeta : Double,
  c : FixedArray[Double],
  off : Int32,
  k0 : Int32,
) -> Double {
  let mut k = k0
  let mut u0 = 0.0
  let mut u1 = 0.0
  let x = 2.0 * (czeta - szeta) * (czeta + szeta)
  while k > 0 {
    k = k - 1
    let t = x * u0 - u1 + c[off + k]
    u1 = u0
    u0 = t
  }
  let f0 = if sinp { 2.0 * szeta * czeta } else { x / 2.0 }
  let fm1 = if sinp { 0.0 } else { 1.0 }
  f0 * u0 - fm1 * u1
}

///|
fn aux_dclenshaw(
  sinp : Bool,
  delta : Double,
  szeta1 : Double,
  czeta1 : Double,
  szeta2 : Double,
  czeta2 : Double,
  c : FixedArray[Double],
  off : Int32,
  k0 : Int32,
) -> Double {
  let d2 = delta * delta
  let czetap = czeta2 * czeta1 - szeta2 * szeta1
  let szetap = szeta2 * czeta1 + czeta2 * szeta1
  let czetam = czeta2 * czeta1 + szeta2 * szeta1
  let szetamd = if delta == 1.0 {
    szeta2 * czeta1 - czeta2 * szeta1
  } else if delta != 0.0 {
    @math.sin(delta) / delta
  } else {
    1.0
  }
  let xa = 2.0 * czetap * czetam
  let xb = -2.0 * szetap * szetamd
  let mut u0a = 0.0
  let mut u0b = 0.0
  let mut u1a = 0.0
  let mut u1b = 0.0
  let mut k = k0 - 1
  while k >= 0 {
    let ta = xa * u0a + d2 * xb * u0b - u1a + c[off + k]
    let tb = xb * u0a + xa * u0b - u1b
    u1a = u0a
    u0a = ta
    u1b = u0b
    u0b = tb
    k = k - 1
  }
  let f0a = (if sinp { szetap } else { czetap }) * czetam
  let f0b = (if sinp { czetap } else { -szetap }) * szetamd
  let fm1a = if sinp { 0.0 } else { 1.0 }
  2.0 * (f0a * u0b + f0b * u0a - fm1a * u1b)
}
