///|
/// Conversions between the coordinate systems used by maps in mainland China.
///
/// - WGS-84: GPS / international standard.
/// - GCJ-02: "Mars coordinates", mandated for maps published in China
///   (AMap/Gaode, Tencent, Google China, Apple Maps China).
/// - BD-09: Baidu Maps, an additional offset applied on top of GCJ-02.
///
/// The forward GCJ-02 transform is the widely published reverse-engineered
/// algorithm (compatible with `coordtransform`, `gcoord`, `coord-convert`
/// inside China). GCJ-02 is the identity wherever `out_of_china` is true:
/// outside the usual bounding box, and on the land of other countries inside
/// it (Seoul, Hanoi, Ulaanbaatar, ...). BD-09 is still applied everywhere.
/// Angles are decimal degrees. These functions do not raise.

///|
/// A geographic position in decimal degrees.
#derive(Eq, Debug)
pub(all) struct LatLon {
  /// Latitude (decimal degrees). Not range-checked.
  lat : Double
  /// Longitude (decimal degrees). Not range-checked.
  lon : Double
}

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
/// Latitude and longitude in decimal degrees. Values are not range-checked.
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
/// Returns true outside the usual GCJ-02 bounding box: longitude
/// 72.004°..137.8347°, latitude 0.8293°..55.8271°. This is the box-only
/// test used by `coordtransform`, `gcoord` and `coord-convert`. The box also
/// covers neighbouring countries; `out_of_china` excludes them.
pub fn out_of_china_bbox(lat : Double, lon : Double) -> Bool {
  lon < 72.004 || lon > 137.8347 || lat < 0.8293 || lat > 55.8271
}

///|
/// Returns true where GCJ-02 equals WGS-84: outside the bounding box (see
/// `out_of_china_bbox`), or on the land of another country inside it.
/// China, including Hong Kong, Macau and Taiwan, and the sea inside the box
/// return false, so GCJ-02 is applied there.
/// Country borders come from Natural Earth 1:10m admin-0, China point of view
/// (public domain), simplified to about 1 km near China. Within about 2 km
/// of a land border the answer may fall on either side.
pub fn out_of_china(lat : Double, lon : Double) -> Bool {
  out_of_china_bbox(lat, lon) || on_foreign_land(lat, lon)
}

///|
/// True when (lat, lon) lies inside one of the foreign land rings in
/// `boundary_data.mbt` (even-odd rule).
fn on_foreign_land(lat : Double, lon : Double) -> Bool {
  let rings = foreign_ring_start.length() - 1
  for r in 0..<rings {
    let b = 4 * r
    if lon < foreign_ring_bbox[b] ||
      lat < foreign_ring_bbox[b + 1] ||
      lon > foreign_ring_bbox[b + 2] ||
      lat > foreign_ring_bbox[b + 3] {
      continue
    }
    let start = foreign_ring_start[r]
    let end = foreign_ring_start[r + 1]
    let mut inside = false
    let mut j = end - 1
    for i in start..<end {
      let xi = foreign_xy[2 * i]
      let yi = foreign_xy[2 * i + 1]
      let xj = foreign_xy[2 * j]
      let yj = foreign_xy[2 * j + 1]
      if (yi > lat) != (yj > lat) &&
        lon < xj + (lat - yj) * (xi - xj) / (yi - yj) {
        inside = !inside
      }
      j = i
    }
    if inside {
      return true
    }
  }
  false
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
/// The GCJ-02 offset applied without the `out_of_china` test.
fn gcj_shift(p : LatLon) -> LatLon {
  let (dlat, dlon) = gcj_delta(p.lat, p.lon)
  { lat: p.lat + dlat, lon: p.lon + dlon, }
}

///|
/// GCJ-02 offsets inside the box are at most about 0.011° (730 m), so a point
/// this far outside the box has no preimage inside it.
fn far_outside_box(p : LatLon) -> Bool {
  let m = 0.05
  p.lon < 72.004 - m ||
  p.lon > 137.8347 + m ||
  p.lat < 0.8293 - m ||
  p.lat > 55.8271 + m
}

///|
/// WGS-84 → GCJ-02 (degrees).
/// Where `out_of_china` is true the result equals the input.
/// Elsewhere this is the published Krasovsky algorithm.
/// Agreement with the official cipher is about 1 m.
pub fn wgs84_to_gcj02(p : LatLon) -> LatLon {
  if out_of_china(p.lat, p.lon) {
    return p
  }
  gcj_shift(p)
}

///|
/// GCJ-02 → WGS-84, one-step approximation (error up to about 1–2 m).
/// The result is the input unless the corrected point is inside China (see
/// `out_of_china`). Prefer `gcj02_to_wgs84` when the round trip must stay
/// under 1e-9°.
pub fn gcj02_to_wgs84_approx(p : LatLon) -> LatLon {
  if far_outside_box(p) {
    return p
  }
  let (dlat, dlon) = gcj_delta(p.lat, p.lon)
  match china_source({ lat: p.lat - dlat, lon: p.lon - dlon, }) {
    Some(w) => w
    None => p
  }
}

///|
/// `w` if it is inside China. A source lying exactly on a boundary can be
/// computed a rounding error outside it; such a point is moved by 1e-10°
/// back to the China side. Otherwise `None`.
fn china_source(w : LatLon) -> LatLon? {
  if !out_of_china(w.lat, w.lon) {
    return Some(w)
  }
  let e = 1.0e-10
  for d in [(e, 0.0), (-e, 0.0), (0.0, e), (0.0, -e)] {
    let q = { lat: w.lat + d.0, lon: w.lon + d.1, }
    if !out_of_china(q.lat, q.lon) {
      return Some(q)
    }
  }
  None
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
/// GCJ-02 → WGS-84 by fixed-point iteration.
/// For every WGS-84 point inside China (see `out_of_china`), a round trip
/// through `wgs84_to_gcj02` stays under 1e-9 degrees, including points next
/// to the box edges and to land borders.
///
/// The offset is solved without the region test, and the solution is used
/// when it lies inside China; otherwise the input is returned unchanged.
/// Next to a boundary the forward map is not one-to-one: inside points move
/// by up to about 730 m, so some GCJ-02 values just outside a boundary have
/// two WGS-84 sources. The inside source is returned for those.
pub fn gcj02_to_wgs84(p : LatLon) -> LatLon {
  if far_outside_box(p) {
    return p
  }
  match china_source(invert(gcj_shift, p)) {
    Some(w) => w
    None => p
  }
}

///|
/// GCJ-02 → BD-09 (Baidu), in degrees.
/// The offset is applied everywhere, including points outside China.
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
/// WGS-84 → BD-09. This is `gcj02_to_bd09(wgs84_to_gcj02(p))`.
/// Outside China only the BD-09 offset remains.
pub fn wgs84_to_bd09(p : LatLon) -> LatLon {
  gcj02_to_bd09(wgs84_to_gcj02(p))
}

///|
/// BD-09 → WGS-84 by the two exact inverses.
/// A round trip through `wgs84_to_bd09` stays under 1e-9 degrees.
pub fn bd09_to_wgs84(p : LatLon) -> LatLon {
  gcj02_to_wgs84(bd09_to_gcj02(p))
}
