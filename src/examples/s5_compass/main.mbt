///|
/// Scenario 5: a survey team in the Western Hills of Beijing has its own
/// position from GPS (WGS-84) and the target from Gaode (GCJ-02). They walk
/// with a magnetic compass and read a UTM grid map. Which compass bearing
/// should they follow, and how far off would they end up without the
/// declination? Run: `moon run src/examples/s5_compass`
fn main {
  let g = @geodesic.Geodesic::wgs84()
  let wmm = @magnetic.MagneticModel::wmm2025()
  let year = @magnetic.decimal_year(2026, 11, 10) catch {
    e => abort("\{e.to_repr()}")
  }
  let start = @cn.LatLon::new(39.9917, 116.1883) // GPS, WGS-84
  let target_gcj = @cn.LatLon::new(40.0091, 116.2050) // Gaode, GCJ-02
  let target = @cn.gcj02_to_wgs84(target_gcj)
  let leg = g.inverse(start.lat, start.lon, target.lat, target.lon)
  let e = wmm.field(start.lat, start.lon, 0.0, year) catch {
    err => abort("\{err.to_repr()}")
  }
  let compass = @magnetic.magnetic_bearing(leg.azi1, e.decl)
  println("start (WGS-84)       \{start}")
  println("target GCJ-02 \{target_gcj} -> WGS-84 \{target}")
  println("distance             \{leg.s12} m")
  println("true bearing         \{leg.azi1} deg")
  println(
    "declination (WMM2025, 2026-11-10) \{e.decl} deg, changing \{e.decl_dot} deg/yr",
  )
  println("compass bearing      \{compass} deg")
  // Walking the true bearing on the compass instead puts the team off
  // to the side by s * sin(D) at the end of the leg.
  let rad = 3.141592653589793 / 180.0
  let miss = (leg.s12 * @math.sin(e.decl * rad)).abs()
  println("miss if declination is ignored: \{miss} m")

  // The same bearing from a UTM map: grid bearing + convergence = true
  // bearing, and the grivation (grid -> magnetic) is D - convergence.
  let u = @utm.from_latlon(start.lat, start.lon) catch {
    err => abort("\{err.to_repr()}")
  }
  let grid_bearing = leg.azi1 - u.gamma
  let griv = @magnetic.grivation(e.decl, u.gamma)
  let compass_from_grid = @magnetic.magnetic_bearing(grid_bearing, griv)
  println(
    "UTM zone \{u.zone}: convergence \{u.gamma} deg, grid bearing \{grid_bearing} deg, grivation \{griv} deg",
  )
  println("compass bearing from the grid map: \{compass_from_grid} deg")
  if (compass_from_grid - compass).abs() > 1.0e-9 {
    abort("grid and true-north routes disagree")
  }

  // Declination across China on the same date.
  println("")
  println("declination across China, 2026-11-10 (WMM2025):")
  let cities = [
    ("Mohe", 52.97, 122.54),
    ("Harbin", 45.80, 126.53),
    ("Beijing", 39.9087, 116.3975),
    ("Shanghai", 31.2304, 121.4737),
    ("Shenzhen", 22.5431, 114.0579),
    ("Lhasa", 29.65, 91.13),
    ("Urumqi", 43.83, 87.62),
    ("Kashgar", 39.47, 75.99),
  ]
  for c in cities {
    let d = wmm.declination(c.1, c.2, 0.0, year) catch {
      err => abort("\{err.to_repr()}")
    }
    println("  \{c.0}: \{d} deg")
  }
}
