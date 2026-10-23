// MoonBit port of GeographicLib's Geocentric (Charles Karney, MIT/X11 License,
// https://geographiclib.sourceforge.io/). Converts geodetic (lat, lon, h) to
// geocentric ECEF (X, Y, Z) and back. This is a close transcription of
// Geocentric.cpp (GeographicLib 2.x).

///|
/// Error raised for an invalid ellipsoid.
pub suberror GeocentricError {
  InvalidEllipsoid(String)
} derive(Debug)

///|
pub extend GeocentricError with @debug.Debug::{to_repr}

///|
/// An ellipsoid of revolution for geocentric (ECEF) conversions.
pub struct Geocentric {
  /// Equatorial radius in metres.
  a : Double
  /// Flattening.
  f : Double
  e2 : Double
  e2m : Double
  e2a : Double
  e4a : Double
  maxrad : Double
}

///|
/// Geocentric ECEF coordinates in metres.
pub(all) struct GeocentricPoint {
  /// X (metres), through the equator at lon = 0.
  x : Double
  /// Y (metres), through the equator at lon = 90°.
  y : Double
  /// Z (metres), through the North Pole.
  z : Double
} derive(Debug)

///|
pub extend GeocentricPoint with @debug.Debug::{to_repr}

///|
/// Geodetic position: latitude and longitude in degrees, ellipsoidal height
/// in metres.
pub(all) struct GeodeticPosition {
  /// Latitude (degrees).
  lat : Double
  /// Longitude (degrees).
  lon : Double
  /// Height above the ellipsoid (metres).
  h : Double
} derive(Debug)

///|
pub extend GeodeticPosition with @debug.Debug::{to_repr}

///|
/// Construct an ellipsoid with equatorial radius `a` (metres) and flattening
/// `f`. Raises `InvalidEllipsoid` when `a` is not finite and positive, or when
/// the polar semi-axis is not positive (`f >= 1`).
pub fn Geocentric::new(
  a : Double,
  f : Double,
) -> Geocentric raise GeocentricError {
  let e2 = f * (2.0 - f)
  let e2m = @geomath.sq(1.0 - f)
  if !(@geomath.is_finite(a) && a > 0.0) {
    raise InvalidEllipsoid("Equatorial radius is not positive")
  }
  if !(@geomath.is_finite(f) && f < 1.0) {
    raise InvalidEllipsoid("Polar semi-axis is not positive")
  }
  {
    a,
    f,
    e2,
    e2m,
    e2a: e2.abs(),
    e4a: @geomath.sq(e2),
    maxrad: 2.0 * a / @geomath.dbl_epsilon,
  }
}

///|
let wgs84_geocentric : Geocentric = Geocentric::new(
  6378137.0,
  1.0 / 298.257223563,
) catch {
  _ => abort("unreachable")
}

///|
/// The WGS-84 ellipsoid.
pub fn Geocentric::wgs84() -> Geocentric {
  wgs84_geocentric
}

///|
let grs80_geocentric : Geocentric = Geocentric::new(
  6378137.0,
  1.0 / 298.257222101,
) catch {
  _ => abort("unreachable")
}

///|
/// The GRS80 ellipsoid (also used by CGCS2000).
pub fn Geocentric::grs80() -> Geocentric {
  grs80_geocentric
}

///|
/// The CGCS2000 ellipsoid (same `a` and `f` as GRS80).
pub fn Geocentric::cgcs2000() -> Geocentric {
  grs80_geocentric
}

///|
/// Equatorial radius `a` in metres.
pub fn Geocentric::equatorial_radius(self : Geocentric) -> Double {
  self.a
}

///|
/// Flattening `f` (dimensionless).
pub fn Geocentric::flattening(self : Geocentric) -> Double {
  self.f
}

///|
/// Build a Geocentric ellipsoid matching `geod`.
pub fn Geocentric::from_geodesic(geod : Geodesic) -> Geocentric {
  Geocentric::new(geod.a, geod.f) catch {
    _ => abort("Geodesic already validated its ellipsoid")
  }
}

///|
/// Local ENU → geocentric rotation matrix (row-major, 9 elements). Columns are
/// the east, north, and up unit vectors in ECEF.
fn geocentric_rotation(
  sphi : Double,
  cphi : Double,
  slam : Double,
  clam : Double,
  m : FixedArray[Double],
) -> Unit {
  // Local X axis (east) in geocentric coords
  m[0] = -slam
  m[3] = clam
  m[6] = 0.0
  // Local Y axis (north) in geocentric coords
  m[1] = -clam * sphi
  m[4] = -slam * sphi
  m[7] = cphi
  // Local Z axis (up) in geocentric coords
  m[2] = clam * cphi
  m[5] = slam * cphi
  m[8] = sphi
}

///|
/// Convert geodetic (lat, lon, h) to geocentric ECEF. `lat` should be in
/// [-90°, 90°]; values outside that range yield NaN. `h` is height above the
/// ellipsoid in metres.
pub fn Geocentric::forward(
  self : Geocentric,
  lat : Double,
  lon : Double,
  h : Double,
) -> GeocentricPoint {
  let (sphi, cphi) = @geomath.sincosd(@geomath.lat_fix(lat))
  let (slam, clam) = @geomath.sincosd(lon)
  let n = self.a / (1.0 - self.e2 * @geomath.sq(sphi)).sqrt()
  let z = (self.e2m * n + h) * sphi
  let mut x = (n + h) * cphi
  let y = x * slam
  x = x * clam
  { x, y, z, }
}

///|
/// Convert geocentric ECEF (X, Y, Z) to geodetic (lat, lon, h). The solution
/// that minimises `|h|` is returned (the closest point on the ellipsoid).
/// Longitude is in [-180°, 180°].
pub fn Geocentric::reverse(
  self : Geocentric,
  x : Double,
  y : Double,
  z : Double,
) -> GeodeticPosition {
  let mut r = @geomath.hypot(x, y)
  let mut slam = if r != 0.0 { y / r } else { 0.0 }
  let mut clam = if r != 0.0 { x / r } else { 1.0 }
  let mut h = @geomath.hypot(r, z)
  let mut sphi = 0.0
  let mut cphi = 1.0
  if h > self.maxrad {
    // Far away: treat the Earth as a point.
    r = @geomath.hypot(x / 2.0, y / 2.0)
    slam = if r != 0.0 { y / 2.0 / r } else { 0.0 }
    clam = if r != 0.0 { x / 2.0 / r } else { 1.0 }
    let hh = @geomath.hypot(z / 2.0, r)
    sphi = z / 2.0 / hh
    cphi = r / hh
  } else if self.e4a == 0.0 {
    let hh = @geomath.hypot(if h == 0.0 { 1.0 } else { z }, r)
    sphi = (if h == 0.0 { 1.0 } else { z }) / hh
    cphi = r / hh
    h = h - self.a
  } else {
    let mut p = @geomath.sq(r / self.a)
    let mut q = self.e2m * @geomath.sq(z / self.a)
    let rr = (p + q - self.e4a) / 6.0
    if self.f < 0.0 {
      let tmp = p
      p = q
      q = tmp
    }
    if !(self.e4a * q == 0.0 && rr <= 0.0) {
      let s = self.e4a * p * q / 4.0
      let r2 = @geomath.sq(rr)
      let r3 = rr * r2
      let disc = s * (2.0 * r3 + s)
      let mut u = rr
      if disc >= 0.0 {
        let mut t3 = s + r3
        t3 = t3 + (if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() })
        let t = @math.cbrt(t3)
        u = u + t + (if t != 0.0 { r2 / t } else { 0.0 })
      } else {
        let ang = @math.atan2((-disc).sqrt(), -(s + r3))
        u = u + 2.0 * rr * @math.cos(ang / 3.0)
      }
      let v = (@geomath.sq(u) + self.e4a * q).sqrt()
      let uv = if u < 0.0 { self.e4a * q / (v - u) } else { u + v }
      let w = @geomath.fmax(0.0, self.e2a * (uv - q) / (2.0 * v))
      let k = uv / ((uv + @geomath.sq(w)).sqrt() + w)
      let k1 = if self.f >= 0.0 { k } else { k - self.e2 }
      let k2 = if self.f >= 0.0 { k + self.e2 } else { k }
      let d = k1 * r / k2
      let hh = @geomath.hypot(z / k1, r / k2)
      sphi = z / k1 / hh
      cphi = r / k2 / hh
      h = (1.0 - self.e2m / k1) * @geomath.hypot(d, z)
    } else {
      let zz = ((if self.f >= 0.0 { self.e4a - p } else { p }) / self.e2m).sqrt()
      let xx = (if self.f < 0.0 { self.e4a - p } else { p }).sqrt()
      let hh = @geomath.hypot(zz, xx)
      sphi = zz / hh
      cphi = xx / hh
      if z < 0.0 {
        sphi = -sphi
      }
      h = -self.a * (if self.f >= 0.0 { self.e2m } else { 1.0 }) * hh / self.e2a
    }
  }
  { lat: @geomath.atan2d(sphi, cphi), lon: @geomath.atan2d(slam, clam), h, }
}
