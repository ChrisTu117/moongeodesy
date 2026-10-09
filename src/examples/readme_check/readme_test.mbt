// Compiles and runs the code snippets from README.md (keep in sync).

///|
test "readme" {
  let g = @geodesic.Geodesic::wgs84()

  // Inverse problem: distance and bearings between two points
  let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737) // Beijing → Shanghai
  // r.s12 = 1066676.4 m, r.azi1 = 152.94°, r.azi2, r.m12, r.area12 ...

  // Direct problem: go 1000 km due east from Shenzhen
  let d = g.direct(22.5431, 114.0579, 90.0, 1000000.0) // d.lat2, d.lon2, d.azi2

  // Geodesic line: waypoints at most 100 km apart, e.g. for drawing a route
  let line = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
  let pts = line.densify(100000.0)

  // Polygon area (m², counter-clockwise positive) and perimeter (m)
  let a = g.polygon_area([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])

  // Other ellipsoids: CGCS2000 / GRS80 built in, or any (a, f)
  let cgcs = @geodesic.Geodesic::cgcs2000()
  let sphere = @geodesic.Geodesic::new(6371000.0, 0.0)

  // Chinese map coordinate systems
  let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → Gaode
  let bd = @cn.gcj02_to_bd09(gcj) // Gaode → Baidu
  let wgs = @cn.bd09_to_wgs84(bd) // exact inverse (round trip < 1e-9°)

  // UTM / UPS / MGRS (functions raise @utm.UtmError on bad input)
  let u = @utm.from_latlon(39.9087, 116.3975) // zone 50N, easting/northing, convergence, scale
  let m = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
  let p = @utm.from_mgrs("50S MK 48502 17797") // centre of the 1 m square
  ignore((r, d, pts, a, cgcs, sphere, wgs, u, m, p))
}
