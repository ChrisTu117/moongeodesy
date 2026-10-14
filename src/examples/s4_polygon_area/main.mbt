///|
/// Scenario 4: land and region areas on the ellipsoid.
/// A surveyor traces a field on Gaode Maps (GCJ-02); we convert the vertices to
/// WGS-84 and compute the geodesic area and perimeter, then compare a large
/// region with Python geographiclib 2.1.
/// Run: `moon run src/examples/s4_polygon_area`
fn main {
  let g = @geodesic.Geodesic::wgs84()

  // 1. A field traced in GCJ-02 (counter-clockwise), converted to WGS-84.
  let field_gcj = [
    @cn.LatLon::new(40.00030, 116.38624),
    @cn.LatLon::new(40.00030, 116.40924),
    @cn.LatLon::new(40.01830, 116.40924),
    @cn.LatLon::new(40.01830, 116.38624),
  ]
  let field = field_gcj.map(fn(p) {
    let w = @cn.gcj02_to_wgs84(p)
    (w.lat, w.lon)
  })
  let r = g.polygon_area(field[:])
  println("field: \{r.num} vertices")
  println("  perimeter : \{r.perimeter} m")
  println(
    "  area      : \{r.area} m2 = \{r.area / 1.0e4} ha = \{r.area / (10000.0 / 15.0)} mu",
  )
  // Using the raw GCJ-02 numbers instead shifts the field but barely changes
  // its area; the offset matters for position, not for size.
  let raw = g.polygon_area(field_gcj.map(fn(p) { (p.lat, p.lon) })[:])
  println(
    "  area from raw GCJ-02 coordinates differs by \{(raw.area - r.area).abs()} m2",
  )

  // 2. Incremental building: what if the surveyor adds a corner bulging west?
  let p = g.polygon()
  for v in field {
    p.add_point(v.0, v.1)
  }
  let more = p.test_point(40.00900, 116.37000)
  println("  with a fifth corner the area would be \{more.area} m2")

  // 3. A continental triangle (Mohe, Kashgar region, Sanya) vs geographiclib.
  let big = [(53.56, 123.26), (36.0, 73.5), (18.16, 109.5)]
  let b = g.polygon_area(big[:])
  let ref_area = 7842200048910.174 // geographiclib 2.1, m2
  let ref_perim = 12434486.692466922 // m
  println("triangle Mohe-Kashgar-Sanya:")
  println(
    "  area      : \{b.area / 1.0e6} km2 (geographiclib \{ref_area / 1.0e6} km2)",
  )
  println(
    "  relative area diff      : \{((b.area - ref_area) / ref_area).abs()}",
  )
  println(
    "  relative perimeter diff : \{((b.perimeter - ref_perim) / ref_perim).abs()}",
  )

  // 4. The same quad in WGS-84 checked against geographiclib directly.
  let quad = [
    (39.9990, 116.3800),
    (39.9990, 116.4030),
    (40.0170, 116.4030),
    (40.0170, 116.3800),
  ]
  let q = g.polygon_area(quad[:])
  let q_ref = 3924960.9460105896
  println(
    "WGS-84 quad area \{q.area} m2, |diff| vs geographiclib \{(q.area - q_ref).abs()} m2",
  )
  if ((b.area - ref_area) / ref_area).abs() > 1.0e-12 ||
    (q.area - q_ref).abs() > 1.0e-6 {
    abort("polygon area disagrees with geographiclib")
  }
}
