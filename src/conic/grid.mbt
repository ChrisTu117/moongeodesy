// EPSG-style map grids on top of the conic projections: a false origin
// (latitude and longitude) and a false easting and northing, as in EPSG
// methods 9801 (Lambert Conic Conformal 1SP), 9802 (Lambert Conic
// Conformal 2SP), and 9822 (Albers Equal Area).

///|
priv enum ConicKind {
  Lambert(LambertConformalConic)
  Albers(AlbersEqualArea)
}

///|
/// A projected coordinate system: a conic projection with a false origin
/// and false easting/northing. Easting and northing are in metres.
pub struct ConicGrid {
  priv kind : ConicKind
  /// Latitude of the false origin (degrees).
  lat_origin : Double
  /// Longitude of the false origin, the central meridian (degrees).
  lon_origin : Double
  /// Easting at the false origin (m).
  false_easting : Double
  /// Northing at the false origin (m).
  false_northing : Double
  priv y_origin : Double
}

///|
/// Grid coordinates of a point.
#derive(Debug)
pub(all) struct GridPoint {
  /// Easting (m).
  easting : Double
  /// Northing (m).
  northing : Double
  /// Meridian convergence (degrees), clockwise from true north to grid north.
  gamma : Double
  /// Point scale factor (see `ConicProjected`).
  k : Double
}

///|
pub extend GridPoint with @debug.Debug::{to_repr}

///|
fn check_origin(lat : Double, lon : Double) -> Unit raise ConicError {
  if !(lat.abs() <= 90.0) || !@geomath.is_finite(lon) {
    raise InvalidOriginLatitude(lat)
  }
}

///|
fn ConicKind::forward(
  self : ConicKind,
  lon0 : Double,
  lat : Double,
  lon : Double,
) -> ConicProjected {
  match self {
    Lambert(p) => p.forward(lon0, lat, lon)
    Albers(p) => p.forward(lon0, lat, lon)
  }
}

///|
fn ConicKind::reverse(
  self : ConicKind,
  lon0 : Double,
  x : Double,
  y : Double,
) -> ConicGeographic {
  match self {
    Lambert(p) => p.reverse(lon0, x, y)
    Albers(p) => p.reverse(lon0, x, y)
  }
}

///|
fn make_grid(
  kind : ConicKind,
  lat_origin : Double,
  lon_origin : Double,
  false_easting : Double,
  false_northing : Double,
) -> ConicGrid {
  let y_origin = kind.forward(lon_origin, lat_origin, lon_origin).y
  { kind, lat_origin, lon_origin, false_easting, false_northing, y_origin, }
}

///|
/// EPSG method 9802, Lambert Conic Conformal (2SP): standard parallels
/// `stdlat1` and `stdlat2`, false origin (`lat_origin`, `lon_origin`), and
/// false easting/northing at that origin.
pub fn ConicGrid::lambert_2sp(
  a : Double,
  f : Double,
  lat_origin~ : Double,
  lon_origin~ : Double,
  stdlat1~ : Double,
  stdlat2~ : Double,
  false_easting~ : Double = 0.0,
  false_northing~ : Double = 0.0,
) -> ConicGrid raise ConicError {
  check_origin(lat_origin, lon_origin)
  let p = LambertConformalConic::new(a, f, stdlat1, stdlat2)
  make_grid(Lambert(p), lat_origin, lon_origin, false_easting, false_northing)
}

///|
/// EPSG method 9801, Lambert Conic Conformal (1SP): natural origin
/// (`lat_origin`, `lon_origin`) with scale `k0` there, and false
/// easting/northing at that origin.
pub fn ConicGrid::lambert_1sp(
  a : Double,
  f : Double,
  lat_origin~ : Double,
  lon_origin~ : Double,
  k0~ : Double = 1.0,
  false_easting~ : Double = 0.0,
  false_northing~ : Double = 0.0,
) -> ConicGrid raise ConicError {
  check_origin(lat_origin, lon_origin)
  let p = LambertConformalConic::new_1sp(a, f, lat_origin, k0~)
  make_grid(Lambert(p), lat_origin, lon_origin, false_easting, false_northing)
}

///|
/// EPSG method 9822, Albers Equal Area: standard parallels `stdlat1` and
/// `stdlat2`, false origin (`lat_origin`, `lon_origin`), and false
/// easting/northing at that origin.
pub fn ConicGrid::albers(
  a : Double,
  f : Double,
  lat_origin~ : Double,
  lon_origin~ : Double,
  stdlat1~ : Double,
  stdlat2~ : Double,
  false_easting~ : Double = 0.0,
  false_northing~ : Double = 0.0,
) -> ConicGrid raise ConicError {
  check_origin(lat_origin, lon_origin)
  let p = AlbersEqualArea::new(a, f, stdlat1, stdlat2)
  make_grid(Albers(p), lat_origin, lon_origin, false_easting, false_northing)
}

///|
/// Grid coordinates of (`lat`, `lon`) in degrees.
pub fn ConicGrid::forward(
  self : ConicGrid,
  lat : Double,
  lon : Double,
) -> GridPoint {
  let p = self.kind.forward(self.lon_origin, lat, lon)
  {
    easting: self.false_easting + p.x,
    northing: self.false_northing + (p.y - self.y_origin),
    gamma: p.gamma,
    k: p.k,
  }
}

///|
/// Geographic coordinates of (`easting`, `northing`).
pub fn ConicGrid::reverse(
  self : ConicGrid,
  easting : Double,
  northing : Double,
) -> ConicGeographic {
  self.kind.reverse(
    self.lon_origin,
    easting - self.false_easting,
    northing - self.false_northing + self.y_origin,
  )
}
