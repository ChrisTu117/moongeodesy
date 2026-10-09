///|
/// Tour of moongeodesy. Run with `moon run src/examples/demo`.
fn main {
  let g = @geodesic.Geodesic::wgs84()

  // 1. Shortest path between Beijing and Shanghai (WGS-84).
  let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737)
  println(
    "Beijing → Shanghai: \{r.s12 / 1000.0} km, initial bearing \{r.azi1}°",
  )

  // 2. Travel 1000 km due east from Shenzhen.
  let d = g.direct(22.5431, 114.0579, 90.0, 1000000.0)
  println("1000 km east of Shenzhen: \{d.lat2}, \{d.lon2}")

  // 3. Points at most 2000 km apart along the Shanghai → San Francisco great ellipse.
  let line = g.inverse_line(31.2304, 121.4737, 37.7749, -122.4194)
  println("Shanghai → San Francisco: \{line.distance() / 1000.0} km")
  for p in line.densify(2000000.0) {
    println("  \{p.0}, \{p.1}")
  }

  // 4. Area of a polygon (counter-clockwise), here a 1°×1° cell at the equator.
  let area = g.polygon_area([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])
  println(
    "1°×1° cell: \{area.area / 1.0e6} km², perimeter \{area.perimeter / 1000.0} km",
  )

  // 5. GPS (WGS-84) → Gaode/Tencent (GCJ-02) → Baidu (BD-09) and back.
  let wgs = @cn.LatLon::new(39.9087, 116.3975)
  let gcj = @cn.wgs84_to_gcj02(wgs)
  let bd = @cn.gcj02_to_bd09(gcj)
  let back = @cn.bd09_to_wgs84(bd)
  println("WGS-84 \{wgs} → GCJ-02 \{gcj} → BD-09 \{bd} → WGS-84 \{back}")

  // 6. UTM and MGRS.
  let u = @utm.from_latlon(39.9087, 116.3975) catch { _ => panic() }
  let m = @utm.to_mgrs(39.9087, 116.3975) catch { _ => panic() }
  println("UTM \{u.to_string()}, MGRS \{m}")
  let p = @utm.from_mgrs("50SMK4850217797") catch { _ => panic() }
  println("MGRS 50SMK4850217797 → \{p.lat}, \{p.lon}")
}
