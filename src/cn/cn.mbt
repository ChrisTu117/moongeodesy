///|
/// Conversions between the coordinate systems used by maps in mainland China.
///
/// - WGS-84: GPS / international standard.
/// - GCJ-02: "Mars coordinates", mandated for maps published in China
///   (AMap/Gaode, Tencent, Google China, Apple Maps China).
/// - BD-09: Baidu Maps, an additional offset applied on top of GCJ-02.
///
/// The forward GCJ-02 transform is the widely published reverse-engineered
/// algorithm (compatible with `coordtransform`, `gcoord`, `coord-convert`); outside mainland China the
/// transform is the identity, as with the official encryption.

///|
/// A geographic position in decimal degrees.
pub(all) struct LatLon {
  lat : Double
  lon : Double
} derive(Eq, Debug)

///|
pub impl Show for LatLon with fn output(self, logger) {
  logger.write_string("LatLon(\{self.lat}, \{self.lon})")
}

///|
pub extend LatLon with Show::{to_string, output}

///|
pub extend LatLon with @debug.Debug::{to_repr}

///|
pub extend LatLon with Eq::{equal, not_equal}

///|
pub fn LatLon::new(lat : Double, lon : Double) -> LatLon {
  { lat, lon, }
}

///|
let pi : Double = 3.14159265358979324

///|
/// Krasovsky 1940 semi-major axis used by GCJ-02.
let gcj_a : Double = 6378245.0

///|
/// Krasovsky 1940 first eccentricity squared used by GCJ-02.
let gcj_ee : Double = 0.00669342162296594323

///|
let x_pi : Double = 3.14159265358979324 * 3000.0 / 180.0

///|
/// Returns true when the point lies outside the rough bounding box of
/// mainland China, where GCJ-02 equals WGS-84.
pub fn out_of_china(lat : Double, lon : Double) -> Bool {
  lon < 72.004 || lon > 137.8347 || lat < 0.8293 || lat > 55.8271
}

///|
fn transform_lat(x : Double, y : Double) -> Double {
  let mut ret = -100.0 +
    2.0 * x +
    3.0 * y +
    0.2 * y * y +
    0.1 * x * y +
    0.2 * x.abs().sqrt()
  ret += (20.0 * @math.sin(6.0 * x * pi) + 20.0 * @math.sin(2.0 * x * pi)) *
    2.0 /
    3.0
  ret += (20.0 * @math.sin(y * pi) + 40.0 * @math.sin(y / 3.0 * pi)) * 2.0 / 3.0
  ret += (160.0 * @math.sin(y / 12.0 * pi) + 320.0 * @math.sin(y * pi / 30.0)) *
    2.0 /
    3.0
  ret
}

///|
fn transform_lon(x : Double, y : Double) -> Double {
  let mut ret = 300.0 +
    x +
    2.0 * y +
    0.1 * x * x +
    0.1 * x * y +
    0.1 * x.abs().sqrt()
  ret += (20.0 * @math.sin(6.0 * x * pi) + 20.0 * @math.sin(2.0 * x * pi)) *
    2.0 /
    3.0
  ret += (20.0 * @math.sin(x * pi) + 40.0 * @math.sin(x / 3.0 * pi)) * 2.0 / 3.0
  ret += (150.0 * @math.sin(x / 12.0 * pi) + 300.0 * @math.sin(x / 30.0 * pi)) *
    2.0 /
    3.0
  ret
}

///|
/// GCJ-02 offset (dlat, dlon) in degrees for a WGS-84 point.
fn gcj_delta(lat : Double, lon : Double) -> (Double, Double) {
  let mut dlat = transform_lat(lon - 105.0, lat - 35.0)
  let mut dlon = transform_lon(lon - 105.0, lat - 35.0)
  let rad_lat = lat / 180.0 * pi
  let s = @math.sin(rad_lat)
  let magic = 1.0 - gcj_ee * s * s
  let sqrt_magic = magic.sqrt()
  dlat = dlat * 180.0 / (gcj_a * (1.0 - gcj_ee) / (magic * sqrt_magic) * pi)
  dlon = dlon * 180.0 / (gcj_a / sqrt_magic * @math.cos(rad_lat) * pi)
  (dlat, dlon)
}

///|
/// WGS-84 → GCJ-02.
pub fn wgs84_to_gcj02(p : LatLon) -> LatLon {
  if out_of_china(p.lat, p.lon) {
    return p
  }
  let (dlat, dlon) = gcj_delta(p.lat, p.lon)
  { lat: p.lat + dlat, lon: p.lon + dlon, }
}

///|
/// GCJ-02 → WGS-84, single-step approximation (error up to ~1–2 m).
pub fn gcj02_to_wgs84_approx(p : LatLon) -> LatLon {
  if out_of_china(p.lat, p.lon) {
    return p
  }
  let (dlat, dlon) = gcj_delta(p.lat, p.lon)
  { lat: p.lat - dlat, lon: p.lon - dlon, }
}

///|
/// Solve f(w) = target for w by fixed-point iteration, where f is a small
/// smooth offset map. Converges to ~1e-12 degrees (sub-millimetre).
fn invert(f : (LatLon) -> LatLon, target : LatLon) -> LatLon {
  let mut lat = target.lat
  let mut lon = target.lon
  for _ in 0..<50 {
    let g = f({ lat, lon, })
    let elat = target.lat - g.lat
    let elon = target.lon - g.lon
    lat = lat + elat
    lon = lon + elon
    if elat.abs() < 1.0e-13 && elon.abs() < 1.0e-13 {
      break
    }
  }
  { lat, lon, }
}

///|
/// GCJ-02 → WGS-84, exact inverse (round-trips to < 1e-9 degrees).
pub fn gcj02_to_wgs84(p : LatLon) -> LatLon {
  if out_of_china(p.lat, p.lon) {
    return p
  }
  invert(wgs84_to_gcj02, p)
}

///|
/// GCJ-02 → BD-09 (Baidu).
pub fn gcj02_to_bd09(p : LatLon) -> LatLon {
  let x = p.lon
  let y = p.lat
  let z = (x * x + y * y).sqrt() + 0.00002 * @math.sin(y * x_pi)
  let theta = @math.atan2(y, x) + 0.000003 * @math.cos(x * x_pi)
  { lat: z * @math.sin(theta) + 0.006, lon: z * @math.cos(theta) + 0.0065, }
}

///|
/// BD-09 → GCJ-02 using Baidu's published approximate inverse formula
/// (error up to ~1e-6 degrees). This matches most existing libraries.
pub fn bd09_to_gcj02_approx(p : LatLon) -> LatLon {
  let x = p.lon - 0.0065
  let y = p.lat - 0.006
  let z = (x * x + y * y).sqrt() - 0.00002 * @math.sin(y * x_pi)
  let theta = @math.atan2(y, x) - 0.000003 * @math.cos(x * x_pi)
  { lat: z * @math.sin(theta), lon: z * @math.cos(theta), }
}

///|
/// BD-09 → GCJ-02, exact inverse (round-trips to < 1e-9 degrees).
pub fn bd09_to_gcj02(p : LatLon) -> LatLon {
  invert(gcj02_to_bd09, p)
}

///|
/// WGS-84 → BD-09.
pub fn wgs84_to_bd09(p : LatLon) -> LatLon {
  gcj02_to_bd09(wgs84_to_gcj02(p))
}

///|
/// BD-09 → WGS-84 (exact).
pub fn bd09_to_wgs84(p : LatLon) -> LatLon {
  gcj02_to_wgs84(bd09_to_gcj02(p))
}
