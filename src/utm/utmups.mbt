// UTM / UPS coordinates (WGS-84), following GeographicLib::UTMUPS conventions.

///|
/// Errors raised by UTM/UPS/MGRS conversions.
pub suberror UtmError {
  /// Latitude outside [-90, 90] or not finite.
  InvalidLatitude(Double)
  /// Zone outside 0..=60 (0 denotes UPS) or inconsistent with the latitude.
  InvalidZone(Int)
  /// Easting/northing outside the allowed range for the zone.
  OutOfRange(String)
  /// Malformed MGRS string.
  InvalidMgrs(String)
} derive(Debug)

///|
pub extend UtmError with @debug.Debug::{to_repr}

///|
/// A UTM or UPS coordinate. `zone` is 1..=60 for UTM and 0 for UPS.
pub(all) struct UtmCoord {
  /// UTM zone 1..=60, or 0 for UPS.
  zone : Int
  /// Northern hemisphere (UTM) or north polar (UPS).
  north : Bool
  /// Easting (m), including the false easting.
  easting : Double
  /// Northing (m), including the false northing.
  northing : Double
  /// Meridian convergence (degrees).
  gamma : Double
  /// Point scale factor (dimensionless).
  k : Double
} derive(Debug)

///|
pub extend UtmCoord with @debug.Debug::{to_repr}

///|
let utm_false_easting : Double = 500000.0

///|
let utm_false_northing_south : Double = 10000000.0

///|
let ups_false_en : Double = 2000000.0

///|
/// Latitude band index in [-10, 9] (band C is -10, X is 9).
fn latitude_band(lat : Double) -> Int {
  let ilat = lat.floor().to_int()
  let b = if ilat >= 0 { ilat / 8 } else { -((-ilat + 7) / 8) }
  if b < -10 {
    -10
  } else if b > 9 {
    9
  } else {
    b
  }
}

///|
/// The standard UTM zone for (lat, lon), including the Norway and Svalbard
/// exceptions; returns 0 (UPS) for lat >= 84 or lat < -80.
pub fn standard_zone(lat : Double, lon : Double) -> Int {
  if lat >= 84.0 || lat < -80.0 {
    return 0
  }
  let mut ilon = @geomath.ang_normalize(lon).floor().to_int()
  if ilon == 180 {
    ilon = -180
  }
  let mut zone = (ilon + 186) / 6
  let band = latitude_band(lat)
  if band == 7 && zone == 31 && ilon >= 3 {
    zone = 32
  } else if band == 9 && ilon >= 0 && ilon < 42 {
    zone = 2 * ((ilon + 183) / 12) + 1
  }
  zone
}

///|
/// Central meridian of a UTM zone (degrees).
pub fn central_meridian(zone : Int) -> Double {
  (6 * zone - 183).to_double()
}

///|
/// Convert geographic coordinates (degrees) to UTM or UPS metres.
/// A latitude outside [-90, 90], including NaN, raises `InvalidLatitude`.
/// The default zone is `standard_zone`: UPS when `lat >= 84` or `lat < -80`,
/// otherwise UTM, including the Norway and Svalbard exceptions.
/// Pass `zone` 1..=60 or 0 (UPS) to force a zone.
/// Any other zone raises `InvalidZone`.
/// Easting and northing include the standard false offsets
/// (UTM easting 500 km; southern UTM northing 10_000 km; UPS 2_000 km).
/// Latitude exactly ±90 is accepted but the UPS easting and northing are NaN.
/// A point a fraction of a degree off the pole is finite.
pub fn from_latlon(
  lat : Double,
  lon : Double,
  zone? : Int,
) -> UtmCoord raise UtmError {
  if !(lat.abs() <= 90.0) {
    raise InvalidLatitude(lat)
  }
  let z = match zone {
    Some(z) => z
    None => standard_zone(lat, lon)
  }
  if z < 0 || z > 60 {
    raise InvalidZone(z)
  }
  let north = lat > 0.0 || (lat == 0.0 && @geomath.copysign(1.0, lat) > 0.0)
  if z == 0 {
    let p = PolarStereographic::ups().forward(north, lat, lon)
    {
      zone: 0,
      north,
      easting: p.x + ups_false_en,
      northing: p.y + ups_false_en,
      gamma: p.gamma,
      k: p.k,
    }
  } else {
    let p = TransverseMercator::utm().forward(central_meridian(z), lat, lon)
    {
      zone: z,
      north,
      easting: p.x + utm_false_easting,
      northing: p.y + (if north { 0.0 } else { utm_false_northing_south }),
      gamma: p.gamma,
      k: p.k,
    }
  }
}

///|
/// Convert UTM or UPS metres back to geographic degrees.
/// `zone` 0 is UPS and 1..=60 is UTM. Any other zone raises `InvalidZone`.
/// A non-finite easting or northing, or one with absolute value
/// >= 1e8 m, raises `OutOfRange`.
/// Longitude ±180 falls in UTM zone 1.
pub fn to_latlon(
  zone : Int,
  north : Bool,
  easting : Double,
  northing : Double,
) -> Geographic raise UtmError {
  if zone < 0 || zone > 60 {
    raise InvalidZone(zone)
  }
  if !(easting.abs() < 1.0e8 && northing.abs() < 1.0e8) {
    raise OutOfRange("easting/northing not finite or too large")
  }
  if zone == 0 {
    PolarStereographic::ups().reverse(
      north,
      easting - ups_false_en,
      northing - ups_false_en,
    )
  } else {
    TransverseMercator::utm().reverse(
      central_meridian(zone),
      easting - utm_false_easting,
      northing - (if north { 0.0 } else { utm_false_northing_south }),
    )
  }
}

///|
/// Human-readable form, e.g. "50N 448252.0 4419382.1" or "UPS N 2000000.0 ...".
pub fn UtmCoord::to_string(self : UtmCoord) -> String {
  let h = if self.north { "N" } else { "S" }
  let e = (self.easting * 10.0).round() / 10.0
  let n = (self.northing * 10.0).round() / 10.0
  if self.zone == 0 {
    "UPS \{h} \{e} \{n}"
  } else {
    "\{self.zone}\{h} \{e} \{n}"
  }
}
