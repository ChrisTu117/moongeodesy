// MoonBit port of GeographicLib's Intersect class (Charles Karney,
// MIT/X11 License, https://geographiclib.sourceforge.io/).
// C. F. F. Karney, "Geodesic intersections", J. Surveying Eng. 150(3),
// 04024005 (2024), https://doi.org/10.1061/JSUED2.SUENG-1483
// preprint https://arxiv.org/abs/2308.00495 (2023).
// Transcription of Intersect.cpp from GeographicLib 2.5 (Closest, Next,
// and Segment). Not a binding. `All` is not included.

///|
/// Error raised when the ellipsoid is too eccentric, or when the
/// intersection iteration does not converge.
pub suberror IntersectError {
  TooEccentric(String)
  NoConvergence(String)
} derive(Debug)

///|
pub extend IntersectError with @debug.Debug::{to_repr}

///|
/// Signed distances from the origins of geodesics X and Y to an
/// intersection, plus a coincidence indicator.
pub(all) struct IntersectPoint {
  /// Signed distance along geodesic X (m).
  x : Double
  /// Signed distance along geodesic Y (m).
  y : Double
  /// 0 when the geodesics merely cross. +1 when they are coincident and
  /// parallel. -1 when they are coincident and antiparallel.
  c : Int32
} derive(Debug)

///|
pub extend IntersectPoint with @debug.Debug::{to_repr}

///|
/// Intersection of two geodesic segments.
pub(all) struct SegmentIntersect {
  /// Signed distance from the first endpoint of segment X (m).
  x : Double
  /// Signed distance from the first endpoint of segment Y (m).
  y : Double
  /// Coincidence indicator. See `IntersectPoint::c`.
  c : Int32
  /// 0 when both segments contain the intersection.
  /// Otherwise `3 * kx + ky`, where `kx` is -1 if `x < 0`, 0 if
  /// `0 <= x <=` the length of X, and 1 if `x` is past the end of X.
  /// `ky` is the same flag for Y.
  segmode : Int32
} derive(Debug)

///|
pub extend SegmentIntersect with @debug.Debug::{to_repr}

///|
/// Geodesic intersections on one ellipsoid.
///
/// The geodesics are the series solution (`Geodesic`), accurate for
/// modest flattening (`|f|` up to about 1/50). The constructor rejects
/// an ellipsoid that is too eccentric for the intersection search.
pub struct Intersect {
  geod : Geodesic
  a : Double
  f : Double
  /// Authalic radius (m).
  r_r : Double
  /// Half the authalic circumference, `pi * r_r`.
  d : Double
  eps : Double
  tol : Double
  delta : Double
  /// Minimum distance between distinct intersections (m).
  mut t1 : Double
  /// Tile spacing for `closest` (m).
  mut d1 : Double
  /// Tile spacing for `next` (m).
  mut d2 : Double
}

///|
priv struct XPoint {
  x : Double
  y : Double
  c : Int32
}

///|
fn x_nan() -> XPoint {
  { x: @double.not_a_number, y: @double.not_a_number, c: 0, }
}

///|
fn xp(x : Double, y : Double) -> XPoint {
  { x, y, c: 0, }
}

///|
fn xpc(x : Double, y : Double, c : Int32) -> XPoint {
  { x, y, c, }
}

///|
fn XPoint::add(self : XPoint, p : XPoint) -> XPoint {
  { x: self.x + p.x, y: self.y + p.y, c: if p.c != 0 { p.c } else { self.c }, }
}

///|
/// L1 distance from the origin, `|x| + |y|`.
fn XPoint::l1(self : XPoint) -> Double {
  self.x.abs() + self.y.abs()
}

///|
/// L1 distance between two points.
fn XPoint::l1_to(self : XPoint, p : XPoint) -> Double {
  (self.x - p.x).abs() + (self.y - p.y).abs()
}

///|
fn Intersect::same_point(self : Intersect, p : XPoint, q : XPoint) -> Bool {
  (p.x - q.x).abs() + (p.y - q.y).abs() <= self.delta
}

///|
/// L1 distance `|x - x0| + |y - y0|` between two intersection parameters.
pub fn Intersect::dist(
  x : Double,
  y : Double,
  x0~ : Double = 0.0,
  y0~ : Double = 0.0,
) -> Double {
  (x - x0).abs() + (y - y0).abs()
}

///|
fn Intersect::conjugate_dist(
  self : Intersect,
  line : GeodesicLine,
  s3 : Double,
  semi : Bool,
  m12~ : Double = 0.0,
  mm12~ : Double = 1.0,
  mm21~ : Double = 1.0,
) -> Double {
  // semi = false: solve m23 = 0. semi = true: solve M23 = 0.
  // Karney, "Algorithms for geodesics", eqs. 31-33.
  let mut s = s3
  for _ in 0..<100 {
    let pos = line.position(s)
    let m13 = pos.m12
    let mm13 = pos.scale12
    let mm31 = pos.scale21
    let m23 = m13 * mm12 - m12 * mm13
    let mm23 = mm13 * mm21 +
      (if m12 == 0.0 { 0.0 } else { (1.0 - mm12 * mm21) * m13 / m12 })
    let mm32 = mm31 * mm12 +
      (if m13 == 0.0 { 0.0 } else { (1.0 - mm13 * mm31) * m12 / m13 })
    let ds = if semi { m23 * mm23 / (1.0 - mm23 * mm32) } else { -m23 / mm32 }
    s = s + ds
    if !(ds.abs() > self.tol) {
      break
    }
  }
  s
}

///|
fn Intersect::dist_polar(self : Intersect, lat1 : Double) -> Double {
  let line = self.geod.line(lat1, 0.0, 0.0)
  self.conjugate_dist(
    line,
    (1.0 + self.f / 2.0) * self.a * @geomath.pi / 2.0,
    true,
  )
}

///|
fn Intersect::polar_b(self : Intersect) -> Double {
  if self.f == 0.0 {
    return self.d
  }
  let mut lat0 = 63.0
  let mut s0 = self.dist_polar(lat0)
  let mut lat1 = 65.0
  let mut s1 = self.dist_polar(lat1)
  let mut lat2 = 64.0
  let mut s2 = self.dist_polar(lat2)
  let mut sx = s2
  for _ in 0..<10 {
    let den = (lat1 - lat0) * s2 + (lat0 - lat2) * s1 + (lat2 - lat1) * s0
    if !(den < 0.0 || den > 0.0) {
      break
    }
    let latn = (
        (lat1 - lat0) * (lat1 + lat0) * s2 +
        (lat0 - lat2) * (lat0 + lat2) * s1 +
        (lat2 - lat1) * (lat2 + lat1) * s0
      ) /
      (2.0 * den)
    lat0 = lat1
    s0 = s1
    lat1 = lat2
    s1 = s2
    lat2 = latn
    s2 = self.dist_polar(lat2)
    let better = if self.f < 0.0 { s2 < sx } else { s2 > sx }
    if better {
      sx = s2
    }
  }
  2.0 * sx
}

///|
fn Intersect::conj_dist(
  self : Intersect,
  azi : Double,
) -> (Double, Double) raise IntersectError {
  let line = self.geod.line(0.0, 0.0, azi)
  let s = self.conjugate_dist(line, self.d, false)
  let p = self.basic(line, line, xp(s / 2.0, -3.0 * s / 2.0))
  (s, p.l1() - 2.0 * s)
}

///|
fn Intersect::dist_oblique(self : Intersect) -> Double raise IntersectError {
  if self.f == 0.0 {
    return self.d
  }
  let (_, ds0_0) = self.conj_dist(46.0)
  let (s1_0, ds1_0) = self.conj_dist(44.0)
  let mut ds0 = ds0_0
  let mut s1 = s1_0
  let mut ds1 = ds1_0
  let mut azi0 = 46.0
  let mut azi1 = 44.0
  let mut sx = s1
  let mut dsx = ds1.abs()
  for _ in 0..<10 {
    if ds1 == ds0 {
      break
    }
    let azin = (azi0 * ds1 - azi1 * ds0) / (ds1 - ds0)
    azi0 = azi1
    ds0 = ds1
    azi1 = azin
    let (s_new, ds_new) = self.conj_dist(azi1)
    s1 = s_new
    ds1 = ds_new
    if ds1.abs() < dsx {
      sx = s1
      dsx = ds1.abs()
      if ds1 == 0.0 {
        break
      }
    }
  }
  sx
}

///|
/// Build an intersection solver for `geod`.
/// Raises `TooEccentric` when the flattening is outside the range the
/// search can handle. WGS-84, GRS80, and CGCS2000 are inside that range.
/// Raises `NoConvergence` when a constructor distance fails to converge.
pub fn Intersect::new(geod : Geodesic) -> Intersect raise IntersectError {
  let a = geod.equatorial_radius()
  let f = geod.flattening()
  let r_r = (geod.ellipsoid_area() / (4.0 * @geomath.pi)).sqrt()
  let d = r_r * @geomath.pi
  let eps = 3.0 * @geomath.dbl_epsilon
  let tol = d * @math.pow(@geomath.dbl_epsilon, 0.75)
  let delta = d * @math.pow(@geomath.dbl_epsilon, 0.2)
  let ix : Intersect = {
    geod,
    a,
    f,
    r_r,
    d,
    eps,
    tol,
    delta,
    t1: 0.0,
    d1: 0.0,
    d2: 0.0,
  }
  // GeographicLib Intersect::Intersect. Tile sizes depend on conjugate
  // distances, which themselves call the intersection iteration.
  let mut t1 = a * (1.0 - f) * @geomath.pi
  let mut t2 = 2.0 * ix.dist_polar(90.0)
  let t5 = 2.0 * geod.inverse(0.0, 0.0, 90.0, 0.0).s12
  let mut t4 = t1
  let t3 = if f > 0.0 {
    ix.dist_oblique()
  } else {
    t4 = ix.polar_b()
    let swap = t1
    t1 = t2
    t2 = swap
    t5
  }
  let d1 = t2 / 2.0
  let d2 = 2.0 * t3 / 3.0
  let d3 = t4 - delta
  if !(d1 < d3 && d2 < d3 && d2 < 2.0 * t1) {
    raise TooEccentric("Ellipsoid too eccentric for Closest")
  }
  ix.t1 = t1
  ix.d1 = d1
  ix.d2 = d2
  ix
}

///|
let wgs84_intersect : Intersect = Intersect::new(
  Geodesic::new(wgs84_a, wgs84_f) catch {
    _ => abort("unreachable")
  },
) catch {
  _ => abort("unreachable")
}

///|
/// Intersections on the WGS-84 ellipsoid.
pub fn Intersect::wgs84() -> Intersect {
  wgs84_intersect
}

///|
fn Intersect::spherical(
  self : Intersect,
  line_x : GeodesicLine,
  line_y : GeodesicLine,
  p : XPoint,
) -> XPoint {
  // Spherical triangle at the current guesses. About 4.3 nm on WGS-84
  // is the coincidence threshold `_eps * _rR`.
  let px = line_x.position(p.x)
  let py = line_y.position(p.y)
  let inv = self.geod.inverse(px.lat2, px.lon2, py.lat2, py.lon2)
  let z = inv.s12
  let (ang_x, dx) = @geomath.ang_diff(px.azi2, inv.azi1)
  let (ang_y, dy) = @geomath.ang_diff(py.azi2, inv.azi2)
  let (ang_xy, dxy) = @geomath.ang_diff(ang_x, ang_y)
  // Sign of the oriented triangle. The error terms keep the sign stable.
  let s = @geomath.copysign(1.0, ang_xy + (dxy + dy - dx))
  let (sin_x, cos_x) = @geomath.sincosde(s * ang_x, s * dx)
  let (sin_y, cos_y) = @geomath.sincosde(s * ang_y, s * dy)
  if z <= self.eps * self.r_r {
    let c = if (sin_x - sin_y).abs() <= self.eps &&
      (cos_x - cos_y).abs() <= self.eps {
      1
    } else if (sin_x + sin_y).abs() <= self.eps &&
      (cos_x + cos_y).abs() <= self.eps {
      -1
    } else {
      0
    }
    return xpc(0.0, 0.0, c)
  }
  if sin_x.abs() <= self.eps && sin_y.abs() <= self.eps {
    let c = if cos_x * cos_y > 0.0 { 1 } else { -1 }
    // Coincident geodesics: put the intersection at the midpoint.
    return xpc(cos_x * z / 2.0, -cos_y * z / 2.0, c)
  }
  let sx = self.r_r *
    @math.atan2(
      sin_y * @math.sin(z / self.r_r),
      sin_y * cos_x * @math.cos(z / self.r_r) - cos_y * sin_x,
    )
  let sy = self.r_r *
    @math.atan2(
      sin_x * @math.sin(z / self.r_r),
      -sin_x * cos_y * @math.cos(z / self.r_r) + cos_x * sin_y,
    )
  xpc(sx, sy, 0)
}

///|
fn Intersect::basic(
  self : Intersect,
  line_x : GeodesicLine,
  line_y : GeodesicLine,
  p0 : XPoint,
) -> XPoint raise IntersectError {
  let mut q = p0
  for _ in 0..<100 {
    let dq = self.spherical(line_x, line_y, q)
    q = q.add(dq)
    // Stop on a coincidence, on a short step, or on NaN.
    if q.c != 0 || !(dq.l1() > self.tol) {
      return q
    }
  }
  raise NoConvergence("Convergence failure in Intersect")
}

///|
fn fix_coincident(p0 : XPoint, p : XPoint) -> XPoint {
  if p.c == 0 {
    return p
  }
  // Slide along the common geodesic so the point is centred on p0.
  let c = p.c.to_double()
  let s = (p0.x + c * p0.y - (p.x + c * p.y)) / 2.0
  p.add(xp(s, c * s))
}

///|
fn segment_mode(sx : Double, sy : Double, p : XPoint) -> Int32 {
  let kx = if p.x < 0.0 { -1 } else if p.x <= sx { 0 } else { 1 }
  let ky = if p.y < 0.0 { -1 } else if p.y <= sy { 0 } else { 1 }
  kx * 3 + ky
}

///|
fn fix_segment(sx : Double, sy : Double, p : XPoint) -> XPoint {
  if p.c == 0 {
    return p
  }
  // Place a coincident intersection inside the segment rectangle when
  // the common geodesic cuts it; otherwise in the smallest gap.
  let c = p.c.to_double()
  let pya = p.y - c * p.x
  let sa = -p.x
  let pyb = p.y - c * (p.x - sx)
  let sb = sx - p.x
  let pxc = p.x - c * p.y
  let sc = c * -p.y
  let pxd = p.x - c * (p.y - sy)
  let sd = c * (sy - p.y)
  let ga = 0.0 <= pya && pya <= sy
  let gb = 0.0 <= pyb && pyb <= sy
  let gc = 0.0 <= pxc && pxc <= sx
  let gd = 0.0 <= pxd && pxd <= sx
  let s = if ga && gb {
    (sa + sb) / 2.0
  } else if gc && gd {
    (sc + sd) / 2.0
  } else if ga && gc {
    (sa + sc) / 2.0
  } else if ga && gd {
    (sa + sd) / 2.0
  } else if gb && gc {
    (sb + sc) / 2.0
  } else if gb && gd {
    (sb + sd) / 2.0
  } else if p.c > 0 {
    if (p.x - p.y + sy).abs() < (p.x - p.y - sx).abs() {
      (sy - (p.x + p.y)) / 2.0
    } else {
      (sx - (p.x + p.y)) / 2.0
    }
  } else if (p.x + p.y).abs() < (p.x + p.y - (sx + sy)).abs() {
    (0.0 - (p.x - p.y)) / 2.0
  } else {
    (sx - sy - (p.x - p.y)) / 2.0
  }
  p.add(xp(s, c * s))
}

///|
fn Intersect::closest_int(
  self : Intersect,
  line_x : GeodesicLine,
  line_y : GeodesicLine,
  p0 : XPoint,
) -> XPoint raise IntersectError {
  let ix = [0, 1, -1, 0, 0]
  let iy = [0, 0, 0, 1, -1]
  let skip = [false, false, false, false, false]
  let mut q = x_nan()
  for n in 0..<5 {
    if skip[n] {
      continue
    }
    let mut qx = self.basic(
      line_x,
      line_y,
      p0.add(xp(ix[n].to_double() * self.d1, iy[n].to_double() * self.d1)),
    )
    qx = fix_coincident(p0, qx)
    if self.same_point(q, qx) {
      continue
    }
    if qx.l1_to(p0) < self.t1 {
      q = qx
      break
    }
    if n == 0 || qx.l1_to(p0) < q.l1_to(p0) {
      q = qx
    }
    let mut m = n + 1
    while m < 5 {
      let tile = xp(
        p0.x + ix[m].to_double() * self.d1,
        p0.y + iy[m].to_double() * self.d1,
      )
      skip[m] = skip[m] || qx.l1_to(tile) < 2.0 * self.t1 - self.d1 - self.delta
      m = m + 1
    }
  }
  q
}

///|
fn Intersect::next_int(
  self : Intersect,
  line_x : GeodesicLine,
  line_y : GeodesicLine,
) -> XPoint raise IntersectError {
  let ix = [-1, -1, 1, 1, -2, 0, 2, 0]
  let iy = [-1, 1, -1, 1, 0, 2, 0, -2]
  let skip = [false, false, false, false, false, false, false, false]
  let z = xp(0.0, 0.0)
  let mut q = xp(@double.infinity, 0.0)
  for n in 0..<8 {
    if skip[n] {
      continue
    }
    let mut qx = self.basic(
      line_x,
      line_y,
      xp(ix[n].to_double() * self.d2, iy[n].to_double() * self.d2),
    )
    qx = fix_coincident(z, qx)
    let zerop = self.same_point(z, qx)
    if qx.c == 0 && zerop {
      continue
    }
    if qx.c != 0 && zerop {
      for sgn in [-1, 1] {
        let s = self.conjugate_dist(line_x, sgn.to_double() * self.d, false)
        let qa = xpc(s, qx.c.to_double() * s, qx.c)
        if qa.l1() < q.l1() {
          q = qa
        }
      }
    } else if qx.l1() < q.l1() {
      q = qx
    }
    for sgn in [-1, 0, 1] {
      if (qx.c == 0 && sgn != 0) || (zerop && sgn == 0) {
        continue
      }
      let qy = if qx.c != 0 {
        qx.add(
          xp(
            sgn.to_double() * self.d2,
            qx.c.to_double() * sgn.to_double() * self.d2,
          ),
        )
      } else {
        qx
      }
      let mut m = n + 1
      while m < 8 {
        let tile = xp(ix[m].to_double() * self.d2, iy[m].to_double() * self.d2)
        skip[m] = skip[m] ||
          qy.l1_to(tile) < 2.0 * self.t1 - self.d2 - self.delta
        m = m + 1
      }
    }
  }
  q
}

///|
fn Intersect::segment_int(
  self : Intersect,
  line_x : GeodesicLine,
  line_y : GeodesicLine,
) -> (XPoint, Int32) raise IntersectError {
  let sx = line_x.distance()
  let sy = line_y.distance()
  let p0 = xp(sx / 2.0, sy / 2.0)
  let mut q = self.closest_int(line_x, line_y, p0)
  q = fix_segment(sx, sy, q)
  let mut segmode = segment_mode(sx, sy, q)
  // GeographicLib also checks the four corners. The conjecture is that
  // this never overrides the midpoint result. The check stays, matching
  // Intersect.cpp (conjectureproved = false).
  if segmode != 0 && p0.l1() >= p0.l1_to(q) {
    let mut segmodex = 1
    let mut qx = x_nan()
    let mut ix = 0
    while ix < 2 && segmodex != 0 {
      let mut iy = 0
      while iy < 2 && segmodex != 0 {
        let t = xp(ix.to_double() * sx, iy.to_double() * sy)
        if q.l1_to(t) >= 2.0 * self.t1 {
          qx = self.basic(line_x, line_y, t)
          qx = fix_coincident(t, qx)
          segmodex = segment_mode(sx, sy, qx)
        }
        iy = iy + 1
      }
      ix = ix + 1
    }
    if segmodex == 0 {
      segmode = 0
      q = qx
    }
  }
  (q, segmode)
}

///|
fn point_of(p : XPoint) -> IntersectPoint {
  { x: p.x, y: p.y, c: p.c, }
}

///|
/// Closest intersection of two geodesics.
///
/// Each geodesic is a point `(lat, lon)` in degrees and an azimuth in
/// degrees, clockwise from north. `x` and `y` are signed metres from
/// those points. The result minimises `|x - x0| + |y - y0|`.
/// `x0` and `y0` default to 0.
///
/// The two geodesics meet: the point at distance `x` on the first
/// geodesic is the point at distance `y` on the second.
pub fn Intersect::closest(
  self : Intersect,
  lat_x : Double,
  lon_x : Double,
  azi_x : Double,
  lat_y : Double,
  lon_y : Double,
  azi_y : Double,
  x0~ : Double = 0.0,
  y0~ : Double = 0.0,
) -> IntersectPoint raise IntersectError {
  let line_x = self.geod.line(lat_x, lon_x, azi_x)
  let line_y = self.geod.line(lat_y, lon_y, azi_y)
  point_of(self.closest_int(line_x, line_y, xp(x0, y0)))
}

///|
/// Intersection of two shortest geodesic segments.
///
/// Each segment is two endpoints in degrees. `x` and `y` are signed
/// metres from the first endpoint of each segment. `segmode` is 0 when
/// the intersection lies on both segments.
///
/// The result is well defined when each segment has one shortest
/// geodesic between its endpoints. If the segments do not intersect,
/// the point is the intersection closest to the two midpoints.
pub fn Intersect::segment(
  self : Intersect,
  lat_x1 : Double,
  lon_x1 : Double,
  lat_x2 : Double,
  lon_x2 : Double,
  lat_y1 : Double,
  lon_y1 : Double,
  lat_y2 : Double,
  lon_y2 : Double,
) -> SegmentIntersect raise IntersectError {
  let line_x = self.geod.inverse_line(lat_x1, lon_x1, lat_x2, lon_x2)
  let line_y = self.geod.inverse_line(lat_y1, lon_y1, lat_y2, lon_y2)
  let (q, segmode) = self.segment_int(line_x, line_y)
  { x: q.x, y: q.y, c: q.c, segmode, }
}

///|
/// Next intersection of two geodesics that already meet at one point.
///
/// Both geodesics start at `(lat, lon)`. `azi_x` and `azi_y` are their
/// azimuths in degrees. The result minimises `|x| + |y|` among
/// intersections other than the starting point.
///
/// Equal closest intersections do occur. This method returns one of them.
pub fn Intersect::next(
  self : Intersect,
  lat : Double,
  lon : Double,
  azi_x : Double,
  azi_y : Double,
) -> IntersectPoint raise IntersectError {
  let line_x = self.geod.line(lat, lon, azi_x)
  let line_y = self.geod.line(lat, lon, azi_y)
  point_of(self.next_int(line_x, line_y))
}
