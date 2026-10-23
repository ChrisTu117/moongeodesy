// MoonBit port of GeographicLib's LocalCartesian (Charles Karney, MIT/X11
// License, https://geographiclib.sourceforge.io/). Converts geodetic
// (lat, lon, h) to a local East-North-Up cartesian frame relative to an origin
// (lat0, lon0, h0), via geocentric ECEF. Transcription of LocalCartesian.cpp.

///|
/// Local East-North-Up cartesian coordinates relative to an origin on an
/// ellipsoid. The z axis is normal to the ellipsoid; y points north; x points
/// east. The plane z = -h0 is tangent to the ellipsoid at the origin.
pub struct LocalCartesian {
  earth : Geocentric
  mut lat0 : Double
  mut lon0 : Double
  mut h0 : Double
  mut x0 : Double
  mut y0 : Double
  mut z0 : Double
  /// Row-major 3×3: columns are local E, N, U expressed in ECEF.
  r : FixedArray[Double]
}

///|
/// Local cartesian ENU point in metres.
pub(all) struct LocalCartesianPoint {
  /// East (metres).
  x : Double
  /// North (metres).
  y : Double
  /// Up, normal to the ellipsoid (metres).
  z : Double
} derive(Debug)

///|
pub extend LocalCartesianPoint with @debug.Debug::{to_repr}

///|
fn local_cartesian_init(
  earth : Geocentric,
  lat0 : Double,
  lon0 : Double,
  h0 : Double,
) -> LocalCartesian {
  let lc : LocalCartesian = {
    earth,
    lat0: 0.0,
    lon0: 0.0,
    h0: 0.0,
    x0: 0.0,
    y0: 0.0,
    z0: 0.0,
    r: FixedArray::make(9, 0.0),
  }
  lc.reset(lat0, lon0, h0~)
  lc
}

///|
/// Local cartesian frame on ellipsoid `(a, f)` with origin `(lat0, lon0, h0)`.
/// `h0` defaults to 0. Raises `GeocentricError` when the ellipsoid is invalid.
pub fn LocalCartesian::new(
  a : Double,
  f : Double,
  lat0 : Double,
  lon0 : Double,
  h0~ : Double = 0.0,
) -> LocalCartesian raise GeocentricError {
  local_cartesian_init(Geocentric::new(a, f), lat0, lon0, h0)
}

///|
/// Local cartesian frame on the same ellipsoid as `geod`.
pub fn LocalCartesian::from_geodesic(
  geod : Geodesic,
  lat0 : Double,
  lon0 : Double,
  h0~ : Double = 0.0,
) -> LocalCartesian {
  local_cartesian_init(Geocentric::from_geodesic(geod), lat0, lon0, h0)
}

///|
/// Local cartesian frame on WGS-84 with origin `(lat0, lon0, h0)`.
pub fn LocalCartesian::wgs84(
  lat0 : Double,
  lon0 : Double,
  h0~ : Double = 0.0,
) -> LocalCartesian {
  local_cartesian_init(Geocentric::wgs84(), lat0, lon0, h0)
}

///|
/// Local cartesian frame on GRS80 / CGCS2000.
pub fn LocalCartesian::grs80(
  lat0 : Double,
  lon0 : Double,
  h0~ : Double = 0.0,
) -> LocalCartesian {
  local_cartesian_init(Geocentric::grs80(), lat0, lon0, h0)
}

///|
/// Move the origin to `(lat0, lon0, h0)`. `lat0` must be in [-90°, 90°].
pub fn LocalCartesian::reset(
  self : LocalCartesian,
  lat0 : Double,
  lon0 : Double,
  h0~ : Double = 0.0,
) -> Unit {
  self.lat0 = @geomath.lat_fix(lat0)
  self.lon0 = @geomath.ang_normalize(lon0)
  self.h0 = h0
  let p = self.earth.forward(self.lat0, self.lon0, self.h0)
  self.x0 = p.x
  self.y0 = p.y
  self.z0 = p.z
  let (sphi, cphi) = @geomath.sincosd(self.lat0)
  let (slam, clam) = @geomath.sincosd(self.lon0)
  geocentric_rotation(sphi, cphi, slam, clam, self.r)
}

///|
/// Origin latitude (degrees).
pub fn LocalCartesian::latitude_origin(self : LocalCartesian) -> Double {
  self.lat0
}

///|
/// Origin longitude (degrees).
pub fn LocalCartesian::longitude_origin(self : LocalCartesian) -> Double {
  self.lon0
}

///|
/// Origin height above the ellipsoid (metres).
pub fn LocalCartesian::height_origin(self : LocalCartesian) -> Double {
  self.h0
}

///|
/// Equatorial radius of the underlying ellipsoid (metres).
pub fn LocalCartesian::equatorial_radius(self : LocalCartesian) -> Double {
  self.earth.a
}

///|
/// Flattening of the underlying ellipsoid.
pub fn LocalCartesian::flattening(self : LocalCartesian) -> Double {
  self.earth.f
}

///|
/// Convert geodetic `(lat, lon, h)` to local ENU `(x, y, z)`. `h` is height
/// above the ellipsoid in metres.
pub fn LocalCartesian::forward(
  self : LocalCartesian,
  lat : Double,
  lon : Double,
  h : Double,
) -> LocalCartesianPoint {
  let p = self.earth.forward(lat, lon, h)
  let xc = p.x - self.x0
  let yc = p.y - self.y0
  let zc = p.z - self.z0
  let r = self.r
  {
    x: r[0] * xc + r[3] * yc + r[6] * zc,
    y: r[1] * xc + r[4] * yc + r[7] * zc,
    z: r[2] * xc + r[5] * yc + r[8] * zc,
  }
}

///|
/// Convert local ENU `(x, y, z)` to geodetic `(lat, lon, h)`.
pub fn LocalCartesian::reverse(
  self : LocalCartesian,
  x : Double,
  y : Double,
  z : Double,
) -> GeodeticPosition {
  let r = self.r
  let xc = self.x0 + r[0] * x + r[1] * y + r[2] * z
  let yc = self.y0 + r[3] * x + r[4] * y + r[5] * z
  let zc = self.z0 + r[6] * x + r[7] * y + r[8] * z
  self.earth.reverse(xc, yc, zc)
}
