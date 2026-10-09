// Transcription of GeographicLib::Math (Charles Karney, MIT/X11 License).
// https://geographiclib.sourceforge.io/

///|
pub let pi : Double = 3.141592653589793

///|
pub let deg_per_rad : Double = 180.0 / pi

///|
pub let rad_per_deg : Double = pi / 180.0

///|
/// Machine epsilon for Double (2^-52).
pub let dbl_epsilon : Double = 2.220446049250313e-16

///|
/// Smallest positive normal Double.
pub let dbl_min : Double = 2.2250738585072014e-308

///|
pub fn sq(x : Double) -> Double {
  x * x
}

///|
pub fn degrees(x : Double) -> Double {
  x * deg_per_rad
}

///|
pub fn radians(x : Double) -> Double {
  x * rad_per_deg
}

///|
pub fn is_finite(x : Double) -> Bool {
  !x.is_nan() && !x.is_inf()
}

///|
let sign_mask : Int64 = -9223372036854775808L

///|
/// Magnitude of `x` with the sign bit of `y` (handles signed zeros).
pub fn copysign(x : Double, y : Double) -> Double {
  let bits = (x.reinterpret_as_int64() & sign_mask.lnot()) |
    (y.reinterpret_as_int64() & sign_mask)
  bits.reinterpret_as_double()
}

///|
/// Round to nearest integer, ties to even (like C `rint` / Python `round`).
pub fn round_even(x : Double) -> Double {
  let f = x.floor()
  let d = x - f
  if d > 0.5 {
    f + 1.0
  } else if d < 0.5 {
    f
  } else if f % 2.0 == 0.0 {
    f
  } else {
    f + 1.0
  }
}

///|
/// IEEE 754 remainder of x/y, result in [-y/2, y/2] (y > 0).
pub fn remainder(x : Double, y : Double) -> Double {
  if !is_finite(x) {
    return @double.not_a_number
  }
  let mut r = x % y // fmod: exact, sign of x
  let h = y / 2.0
  if r > h {
    r = r - y
  } else if r < -h {
    r = r + y
  } else if r.abs() == h {
    // tie: pick the result whose quotient is even
    let n = (x - r) / y
    if n % 2.0 != 0.0 {
      r = r - copysign(y, r)
    }
  }
  r
}

///|
/// Error-free transformation of a sum: u + v = s + t exactly.
pub fn sum2(u : Double, v : Double) -> (Double, Double) {
  let s = u + v
  let mut up = s - v
  let mut vpp = s - up
  up = up - u
  vpp = vpp - v
  let t = if s == 0.0 { s } else { 0.0 - (up + vpp) }
  (s, t)
}

///|
/// Evaluate polynomial of degree n with coefficients p[s..s+n] by Horner.
pub fn polyval(n : Int, p : FixedArray[Double], s : Int, x : Double) -> Double {
  if n < 0 {
    return 0.0
  }
  let mut y = p[s]
  let mut k = s
  for _ in 0..<n {
    k += 1
    y = y * x + p[k]
  }
  y
}

///|
/// Round an angle so that small values underflow to zero.
pub fn ang_round(x : Double) -> Double {
  let z = 1.0 / 16.0
  let mut y = x.abs()
  if y < z {
    y = z - (z - y)
  }
  copysign(y, x)
}

///|
/// Reduce angle to [-180, 180].
pub fn ang_normalize(x : Double) -> Double {
  let y = remainder(x, 360.0)
  if y.abs() == 180.0 {
    copysign(180.0, x)
  } else {
    y
  }
}

///|
/// Replace latitudes outside [-90, 90] by NaN.
pub fn lat_fix(x : Double) -> Double {
  if x.abs() > 90.0 {
    @double.not_a_number
  } else {
    x
  }
}

///|
/// Compute y - x reduced to [-180, 180] accurately; returns (d, error).
pub fn ang_diff(x : Double, y : Double) -> (Double, Double) {
  let (d0, t0) = sum2(remainder(-x, 360.0), remainder(y, 360.0))
  let (d1, t) = sum2(remainder(d0, 360.0), t0)
  let mut d = d1
  if d == 0.0 || d.abs() == 180.0 {
    d = copysign(d, if t == 0.0 { y - x } else { -t })
  }
  (d, t)
}

///|
fn quadrant(s : Double, c : Double, q : Int) -> (Double, Double) {
  match (q % 4 + 4) % 4 {
    1 => (c, -s)
    2 => (-s, -c)
    3 => (-c, s)
    _ => (s, c)
  }
}

///|
/// Sine and cosine of x in degrees, exact for multiples of 90.
pub fn sincosd(x : Double) -> (Double, Double) {
  let mut r = if is_finite(x) { x % 360.0 } else { @double.not_a_number }
  let q = if r.is_nan() { 0 } else { round_even(r / 90.0).to_int() }
  r = r - 90.0 * q.to_double()
  r = radians(r)
  let (s0, c0) = quadrant(@math.sin(r), @math.cos(r), q)
  let c = c0 + 0.0
  let s = if s0 == 0.0 { copysign(s0, x) } else { s0 }
  (s, c)
}

///|
/// Sine and cosine of (x + t) in degrees with x in [-180, 180].
pub fn sincosde(x : Double, t : Double) -> (Double, Double) {
  let q = if is_finite(x) { round_even(x / 90.0).to_int() } else { 0 }
  let r = radians(ang_round(x - 90.0 * q.to_double() + t))
  let (s0, c0) = quadrant(@math.sin(r), @math.cos(r), q)
  let c = c0 + 0.0
  let s = if s0 == 0.0 { copysign(s0, x) } else { s0 }
  (s, c)
}

///|
/// atan2(y, x) in degrees, exact for the quadrant boundaries.
pub fn atan2d(y0 : Double, x0 : Double) -> Double {
  let mut x = x0
  let mut y = y0
  let mut q = 0
  if y.abs() > x.abs() {
    q = 2
    let t = x
    x = y
    y = t
  }
  if x < 0.0 {
    q += 1
    x = -x
  }
  let ang = degrees(@math.atan2(y, x))
  match q {
    1 => copysign(180.0, y) - ang
    2 => 90.0 - ang
    3 => -90.0 + ang
    _ => ang
  }
}

///|
/// Normalize a two-vector.
pub fn norm2(x : Double, y : Double) -> (Double, Double) {
  let r = hypot(x, y)
  (x / r, y / r)
}

///|
/// max with Python/C++ std::max semantics (returns `a` unless `b > a`).
pub fn fmax(a : Double, b : Double) -> Double {
  if b > a {
    b
  } else {
    a
  }
}

///|
/// min with Python/C++ std::min semantics (returns `a` unless `b < a`).
pub fn fmin(a : Double, b : Double) -> Double {
  if b < a {
    b
  } else {
    a
  }
}

///|
pub fn imin(a : Int, b : Int) -> Int {
  if b < a {
    b
  } else {
    a
  }
}

///|
/// Dekker/Veltkamp error-free product: a * b = p + e exactly.
fn two_prod(a : Double, b : Double) -> (Double, Double) {
  let p = a * b
  let split = 134217729.0 // 2^27 + 1
  let ca = split * a
  let ah = ca - (ca - a)
  let al = a - ah
  let cb = split * b
  let bh = cb - (cb - b)
  let bl = b - bh
  let e = ah * bh - p + ah * bl + al * bh + al * bl
  (p, e)
}

///|
/// Nearly correctly rounded hypot (like Python 3.10+ math.hypot). The
/// geodesic code takes differences of products of normalized vectors, so
/// last-bit errors here are amplified for very short lines.
pub fn hypot(x0 : Double, y0 : Double) -> Double {
  let x = x0.abs()
  let y = y0.abs()
  if x.is_inf() || y.is_inf() {
    return @double.infinity
  }
  if x.is_nan() || y.is_nan() {
    return @double.not_a_number
  }
  let (big, small) = if x >= y { (x, y) } else { (y, x) }
  if big == 0.0 {
    return 0.0
  }
  // keep the scaled values in a safe range
  if big > 1.0e150 || big < 1.0e-150 {
    return @math.hypot(x, y)
  }
  let (p1, e1) = two_prod(big, big)
  let (p2, e2) = two_prod(small, small)
  let (s, t) = sum2(p1, p2)
  let lo = t + e1 + e2
  let h = (s + lo).sqrt()
  // one Newton correction using the exact residual (s + lo) - h^2
  let (hh, he) = two_prod(h, h)
  let r = s - hh + lo - he
  h + r / (2.0 * h)
}
