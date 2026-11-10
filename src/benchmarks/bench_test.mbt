///|
/// Benchmarks for the hot paths. Run: `moon bench -p ChrisTu117/moongeodesy/benchmarks`
let g : @geodesic.Geodesic = @geodesic.Geodesic::wgs84()

///|
test "geodesic inverse Beijing-Shanghai" (b : @bench.T) {
  b.bench(fn() { b.keep(g.inverse(39.9087, 116.3975, 31.2304, 121.4737).s12) })
}

///|
test "geodesic inverse nearly antipodal" (b : @bench.T) {
  b.bench(fn() { b.keep(g.inverse(0.0, 0.0, 0.5, 179.5).s12) })
}

///|
test "geodesic direct 1000 km" (b : @bench.T) {
  b.bench(fn() { b.keep(g.direct(22.5431, 114.0579, 90.0, 1.0e6).lat2) })
}

///|
test "polygon area 4 vertices" (b : @bench.T) {
  let pts = [(39.0, 116.0), (39.0, 117.0), (40.0, 117.0), (40.0, 116.0)]
  b.bench(fn() { b.keep(g.polygon_area(pts[:]).area) })
}

///|
test "cn wgs84_to_gcj02" (b : @bench.T) {
  let p = @cn.LatLon::new(39.9087, 116.3975)
  b.bench(fn() { b.keep(@cn.wgs84_to_gcj02(p).lat) })
}

///|
test "cn gcj02_to_wgs84 (iterative)" (b : @bench.T) {
  let p = @cn.LatLon::new(39.9101, 116.40374)
  b.bench(fn() { b.keep(@cn.gcj02_to_wgs84(p).lat) })
}

///|
test "utm from_latlon" (b : @bench.T) {
  b.bench(fn() {
    let u = @utm.from_latlon(39.9087, 116.3975) catch { _ => panic() }
    b.keep(u.zone)
  })
}

///|
test "mgrs encode precision 5" (b : @bench.T) {
  b.bench(fn() {
    let s = @utm.to_mgrs(39.9087, 116.3975) catch { _ => panic() }
    b.keep(s.length())
  })
}

///|
let wmm : @magnetic.MagneticModel = @magnetic.MagneticModel::wmm2025()

///|
test "wmm2025 field (all elements and rates)" (b : @bench.T) {
  b.bench(fn() {
    let e = wmm.field(39.9087, 116.3975, 0.0, 2026.5) catch { _ => panic() }
    b.keep(e.decl)
  })
}

///|
let albers : @conic.AlbersEqualArea = @conic.AlbersEqualArea::new(
  6378137.0,
  1.0 / 298.257222101,
  25.0,
  47.0,
) catch {
  _ => panic()
}

///|
test "albers forward" (b : @bench.T) {
  b.bench(fn() { b.keep(albers.forward(105.0, 39.9087, 116.3975).x) })
}

///|
let lcc : @conic.LambertConformalConic = @conic.LambertConformalConic::new(
  6378137.0,
  1.0 / 298.257222101,
  25.0,
  47.0,
) catch {
  _ => panic()
}

///|
test "lcc forward" (b : @bench.T) {
  b.bench(fn() { b.keep(lcc.forward(105.0, 39.9087, 116.3975).x) })
}
