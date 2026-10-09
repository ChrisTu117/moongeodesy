// Port of GeographicLib's PolygonArea and Accumulator (Charles Karney,
// MIT/X11 License). Computes the area and perimeter of geodesic polygons.

///|
/// Double-double accumulator for accurate summation.
struct Accumulator {
  mut s : Double
  mut t : Double
}

///|
fn Accumulator::add(self : Accumulator, y0 : Double) -> Unit {
  let (y, u0) = sum2(y0, self.t)
  let (s, t) = sum2(y, self.s)
  self.s = s
  self.t = t
  if self.s == 0.0 {
    self.s = u0
  } else {
    self.t = self.t + u0
  }
}

///|
fn Accumulator::sum(self : Accumulator, y : Double) -> Double {
  if y == 0.0 {
    return self.s
  }
  let b = { s: self.s, t: self.t, }
  b.add(y)
  b.s
}

///|
fn Accumulator::negate(self : Accumulator) -> Unit {
  self.s = -self.s
  self.t = -self.t
}

///|
fn Accumulator::remainder(self : Accumulator, y : Double) -> Unit {
  self.s = remainder(self.s, y)
  self.add(0.0)
}

///|
/// Result of a polygon computation.
pub(all) struct PolygonResult {
  /// Number of vertices.
  num : Int
  /// Perimeter (m); for a polyline, its length.
  perimeter : Double
  /// Area (m²), NaN for a polyline. Counter-clockwise traversal gives a
  /// positive area unless `reverse=true`.
  area : Double
} derive(Debug)

///|
pub extend PolygonResult with @debug.Debug::{to_repr}

///|
/// Incrementally built geodesic polygon (or polyline).
pub struct PolygonArea {
  geod : Geodesic
  polyline : Bool
  area0 : Double
  areasum : Accumulator
  perimetersum : Accumulator
  mut num : Int
  mut lat1 : Double
  mut lon1 : Double
  mut crossings : Int
  mut lat0 : Double
  mut lon0 : Double
}

///|
/// Start a new polygon (or polyline when `polyline=true`).
pub fn Geodesic::polygon(
  self : Geodesic,
  polyline? : Bool = false,
) -> PolygonArea {
  let nan = @double.not_a_number
  {
    geod: self,
    polyline,
    area0: 4.0 * pi * self.c2,
    areasum: { s: 0.0, t: 0.0, },
    perimetersum: { s: 0.0, t: 0.0, },
    num: 0,
    lat1: nan,
    lon1: nan,
    crossings: 0,
    lat0: nan,
    lon0: nan,
  }
}

///|
fn transit(lon1_in : Double, lon2_in : Double) -> Int {
  let (lon12, _) = ang_diff(lon1_in, lon2_in)
  let lon1 = ang_normalize(lon1_in)
  let lon2 = ang_normalize(lon2_in)
  if lon12 > 0.0 && ((lon1 < 0.0 && lon2 >= 0.0) || (lon1 > 0.0 && lon2 == 0.0)) {
    1
  } else if lon12 < 0.0 && lon2 < 0.0 && lon1 >= 0.0 {
    -1
  } else {
    0
  }
}

///|
fn transitdirect(lon1_in : Double, lon2_in : Double) -> Int {
  let lon1 = remainder(lon1_in, 720.0)
  let lon2 = remainder(lon2_in, 720.0)
  (if lon2 >= 0.0 && lon2 < 360.0 { 0 } else { 1 }) -
  (if lon1 >= 0.0 && lon1 < 360.0 { 0 } else { 1 })
}

///|
fn areareduce_a(
  area : Accumulator,
  area0 : Double,
  crossings : Int,
  reverse : Bool,
  sign : Bool,
) -> Double {
  area.remainder(area0)
  if (crossings & 1) != 0 {
    area.add((if area.sum(0.0) < 0.0 { 1.0 } else { -1.0 }) * area0 / 2.0)
  }
  if !reverse {
    area.negate()
  }
  if sign {
    if area.sum(0.0) > area0 / 2.0 {
      area.add(-area0)
    } else if area.sum(0.0) <= -area0 / 2.0 {
      area.add(area0)
    }
  } else if area.sum(0.0) >= area0 {
    area.add(-area0)
  } else if area.sum(0.0) < 0.0 {
    area.add(area0)
  }
  0.0 + area.sum(0.0)
}

///|
fn areareduce_b(
  area_in : Double,
  area0 : Double,
  crossings : Int,
  reverse : Bool,
  sign : Bool,
) -> Double {
  let mut area = remainder(area_in, area0)
  if (crossings & 1) != 0 {
    area = area + (if area < 0.0 { 1.0 } else { -1.0 }) * area0 / 2.0
  }
  if !reverse {
    area = -area
  }
  if sign {
    if area > area0 / 2.0 {
      area = area - area0
    } else if area <= -area0 / 2.0 {
      area = area + area0
    }
  } else if area >= area0 {
    area = area - area0
  } else if area < 0.0 {
    area = area + area0
  }
  0.0 + area
}

///|
/// Remove all points.
pub fn PolygonArea::clear(self : PolygonArea) -> Unit {
  let nan = @double.not_a_number
  self.num = 0
  self.crossings = 0
  self.areasum.s = 0.0
  self.areasum.t = 0.0
  self.perimetersum.s = 0.0
  self.perimetersum.t = 0.0
  self.lat0 = nan
  self.lon0 = nan
  self.lat1 = nan
  self.lon1 = nan
}

///|
/// Add a vertex (degrees).
pub fn PolygonArea::add_point(
  self : PolygonArea,
  lat : Double,
  lon : Double,
) -> Unit {
  if self.num == 0 {
    self.lat0 = lat
    self.lat1 = lat
    self.lon0 = lon
    self.lon1 = lon
  } else {
    let r = self.geod.gen_inverse(self.lat1, self.lon1, lat, lon)
    self.perimetersum.add(r.s12)
    if !self.polyline {
      self.areasum.add(r.ss12)
      self.crossings += transit(self.lon1, lon)
    }
    self.lat1 = lat
    self.lon1 = lon
  }
  self.num += 1
}

///|
/// Add an edge given by azimuth `azi` (degrees) and length `s` (m)
/// from the current last vertex.
pub fn PolygonArea::add_edge(
  self : PolygonArea,
  azi : Double,
  s : Double,
) -> Unit {
  if self.num != 0 {
    let r = self.geod.direct(self.lat1, self.lon1, azi, s, long_unroll=true)
    self.perimetersum.add(s)
    if !self.polyline {
      self.areasum.add(r.area12)
      self.crossings += transitdirect(self.lon1, r.lon2)
    }
    self.lat1 = r.lat2
    self.lon1 = r.lon2
    self.num += 1
  }
}

///|
/// Compute perimeter and area of the polygon so far. With `sign=true`
/// (default) the area is signed (counter-clockwise positive); with
/// `reverse=true`, clockwise traversal counts as positive.
pub fn PolygonArea::compute(
  self : PolygonArea,
  reverse? : Bool = false,
  sign? : Bool = true,
) -> PolygonResult {
  if self.num < 2 {
    return {
      num: self.num,
      perimeter: 0.0,
      area: if self.polyline {
        @double.not_a_number
      } else {
        0.0
      },
    }
  }
  if self.polyline {
    return {
      num: self.num,
      perimeter: self.perimetersum.sum(0.0),
      area: @double.not_a_number,
    }
  }
  let r = self.geod.gen_inverse(self.lat1, self.lon1, self.lat0, self.lon0)
  let perimeter = self.perimetersum.sum(r.s12)
  let tempsum = { s: self.areasum.s, t: self.areasum.t, }
  tempsum.add(r.ss12)
  let crossings = self.crossings + transit(self.lon1, self.lon0)
  let area = areareduce_a(tempsum, self.area0, crossings, reverse, sign)
  { num: self.num, perimeter, area, }
}

///|
/// Result if vertex (lat, lon) were added, without modifying the polygon.
pub fn PolygonArea::test_point(
  self : PolygonArea,
  lat : Double,
  lon : Double,
  reverse? : Bool = false,
  sign? : Bool = true,
) -> PolygonResult {
  if self.num == 0 {
    return {
      num: 1,
      perimeter: 0.0,
      area: if self.polyline {
        @double.not_a_number
      } else {
        0.0
      },
    }
  }
  let mut perimeter = self.perimetersum.sum(0.0)
  let mut tempsum = if self.polyline { 0.0 } else { self.areasum.sum(0.0) }
  let mut crossings = self.crossings
  let num = self.num + 1
  let n = if self.polyline { 1 } else { 2 }
  for i in 0..<n {
    let (la1, lo1, la2, lo2) = if i == 0 {
      (self.lat1, self.lon1, lat, lon)
    } else {
      (lat, lon, self.lat0, self.lon0)
    }
    let r = self.geod.gen_inverse(la1, lo1, la2, lo2)
    perimeter = perimeter + r.s12
    if !self.polyline {
      tempsum = tempsum + r.ss12
      crossings += transit(lo1, lo2)
    }
  }
  if self.polyline {
    return { num, perimeter, area: @double.not_a_number, }
  }
  let area = areareduce_b(tempsum, self.area0, crossings, reverse, sign)
  { num, perimeter, area, }
}

///|
/// Result if an edge (azi, s) were added, without modifying the polygon.
pub fn PolygonArea::test_edge(
  self : PolygonArea,
  azi : Double,
  s : Double,
  reverse? : Bool = false,
  sign? : Bool = true,
) -> PolygonResult {
  let nan = @double.not_a_number
  if self.num == 0 {
    return { num: 0, perimeter: nan, area: nan, }
  }
  let num = self.num + 1
  let mut perimeter = self.perimetersum.sum(0.0) + s
  if self.polyline {
    return { num, perimeter, area: nan, }
  }
  let mut tempsum = self.areasum.sum(0.0)
  let mut crossings = self.crossings
  let d = self.geod.direct(self.lat1, self.lon1, azi, s, long_unroll=true)
  tempsum = tempsum + d.area12
  crossings += transitdirect(self.lon1, d.lon2)
  let r = self.geod.gen_inverse(d.lat2, d.lon2, self.lat0, self.lon0)
  perimeter = perimeter + r.s12
  tempsum = tempsum + r.ss12
  crossings += transit(d.lon2, self.lon0)
  let area = areareduce_b(tempsum, self.area0, crossings, reverse, sign)
  { num, perimeter, area, }
}

///|
/// Convenience: area (m², counter-clockwise positive) and perimeter (m) of a
/// closed polygon given as (lat, lon) vertices in degrees.
pub fn Geodesic::polygon_area(
  self : Geodesic,
  points : ArrayView[(Double, Double)],
) -> PolygonResult {
  let p = self.polygon()
  for pt in points {
    p.add_point(pt.0, pt.1)
  }
  p.compute()
}
