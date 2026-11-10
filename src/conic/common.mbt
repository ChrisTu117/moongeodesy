// Shared definitions for the conic projections.

///|
/// Errors raised when a conic projection is constructed.
#derive(Debug)
pub suberror ConicError {
  /// Equatorial radius not finite and positive, or flattening not below 1.
  InvalidEllipsoid(Double, Double)
  /// Standard parallel outside [-90, 90] or not finite.
  InvalidStandardLatitude(Double)
  /// The two standard parallels are incompatible (for example opposite
  /// poles, or only one of them at a pole for the conformal projection).
  IncompatibleStandardLatitudes(Double, Double)
  /// Scale factor not finite and positive.
  InvalidScale(Double)
  /// Origin latitude outside [-90, 90] or not finite.
  InvalidOriginLatitude(Double)
}

///|
pub extend ConicError with @debug.Debug::{to_repr}

///|
/// Result of a forward projection.
#derive(Debug)
pub(all) struct ConicProjected {
  /// Easting (m).
  x : Double
  /// Northing (m).
  y : Double
  /// Meridian convergence (degrees): the bearing of grid north measured
  /// clockwise from true north.
  gamma : Double
  /// Point scale factor. For the conformal projection it is the same in
  /// every direction; for the equal-area projection it is the scale along
  /// the parallel (the scale along the meridian is 1 / k).
  k : Double
}

///|
pub extend ConicProjected with @debug.Debug::{to_repr}

///|
/// Result of a reverse projection.
#derive(Debug)
pub(all) struct ConicGeographic {
  /// Latitude (degrees).
  lat : Double
  /// Longitude (degrees), in [-180, 180].
  lon : Double
  /// Meridian convergence (degrees).
  gamma : Double
  /// Point scale factor, as in `ConicProjected`.
  k : Double
}

///|
pub extend ConicGeographic with @debug.Debug::{to_repr}

///|
let eps : Double = @geomath.dbl_epsilon

///|
/// eps^2, used to keep cos(phi) away from 0 at the poles.
let epsx : Double = @geomath.dbl_epsilon * @geomath.dbl_epsilon

///|
fn hyp(x : Double) -> Double {
  @geomath.hypot(1.0, x)
}

///|
fn sq(x : Double) -> Double {
  x * x
}

///|
fn check_ellipsoid(a : Double, f : Double) -> Unit raise ConicError {
  if !(@geomath.is_finite(a) && a > 0.0) || !(@geomath.is_finite(f) && f < 1.0) {
    raise InvalidEllipsoid(a, f)
  }
}

///|
fn check_scale(k : Double) -> Unit raise ConicError {
  if !(@geomath.is_finite(k) && k > 0.0) {
    raise InvalidScale(k)
  }
}

///|
fn check_stdlat(lat : Double) -> Unit raise ConicError {
  if !(lat.abs() <= 90.0) {
    raise InvalidStandardLatitude(lat)
  }
}

// Divided differences, Df(x, y) = (f(x) - f(y)) / (x - y), evaluated
// without cancellation. See W. M. Kahan and R. J. Fateman, "Symbolic
// computation of divided differences", SIGSAM Bull. 33(2), 7-28 (1999), and
// GeographicLib::LambertConformalConic.

///|
/// D hyp, given hx = hyp(x) and hy = hyp(y).
fn dhyp(x : Double, y : Double, hx : Double, hy : Double) -> Double {
  (x + y) / (hx + hy)
}

///|
/// D sn for sn(x) = x / hyp(x), given sx = sn(x) and sy = sn(y).
fn dsn(x : Double, y : Double, sx : Double, sy : Double) -> Double {
  let t = x * y
  if t > 0.0 {
    (x + y) * sq(sx * sy / t) / (sx + sy)
  } else if x - y != 0.0 {
    (sx - sy) / (x - y)
  } else {
    1.0
  }
}

///|
/// D log1p.
fn dlog1p(x : Double, y0 : Double) -> Double {
  let mut t = x - y0
  let mut y = y0
  if t < 0.0 {
    t = -t
    y = x
  }
  if t != 0.0 {
    @math.ln_1p(t / (1.0 + y)) / t
  } else {
    1.0 / (1.0 + x)
  }
}

///|
/// D exp.
fn dexp(x : Double, y : Double) -> Double {
  let t = (x - y) / 2.0
  (if t != 0.0 { @math.sinh(t) / t } else { 1.0 }) * @math.exp((x + y) / 2.0)
}

///|
/// D sinh, given sx = sinh(x), cx = cosh(x), and the same for y.
fn dsinh(
  x : Double,
  y : Double,
  sx : Double,
  sy : Double,
  cx : Double,
  cy : Double,
) -> Double {
  let t = (x - y) / 2.0
  (if t != 0.0 { @math.sinh(t) / t } else { 1.0 }) *
  ((sx * sy + cx * cy + 1.0) / 2.0).sqrt()
}

///|
/// D asinh, given hx = hyp(x) and hy = hyp(y).
fn dasinh(x : Double, y : Double, hx : Double, hy : Double) -> Double {
  let t = x - y
  if t != 0.0 {
    @math.asinh(
      if x * y > 0.0 {
        t * (x + y) / (x * hy + y * hx)
      } else {
        x * hy - y * hx
      },
    ) /
    t
  } else {
    1.0 / hx
  }
}
