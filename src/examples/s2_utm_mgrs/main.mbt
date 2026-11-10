///|
/// Scenario 2: field survey points (incl. Norway/Svalbard exceptions and a
/// polar UPS point) go WGS-84 -> UTM/UPS -> MGRS -> back, with residuals.
/// Run: `moon run src/examples/s2_utm_mgrs`
fn main {
  let g = @geodesic.Geodesic::wgs84()
  // Expected zone (0 = UPS) and MGRS grid zone designator.
  let pts : Array[(String, Double, Double, Int32, String)] = [
    ("Beijing", 39.9087, 116.3975, 50, "50S"),
    ("Sydney", -33.8688, 151.2093, 56, "56H"),
    ("Bergen (zone 32V exception)", 60.3913, 5.3221, 32, "32V"),
    ("Svalbard (zone 33X exception)", 78.2232, 15.6267, 33, "33X"),
    ("South Pole station (UPS)", -89.9978, 139.2729, 0, "B"),
  ]
  for p in pts {
    let (name, lat, lon, zone, gzd) = p
    show(g, name, lat, lon, zone, gzd) catch {
      e => abort("\{name}: conversion failed: \{e}")
    }
  }
}

///|
fn show(
  g : @geodesic.Geodesic,
  name : String,
  lat : Double,
  lon : Double,
  zone : Int32,
  gzd : String,
) -> Unit raise {
  let u = @utm.from_latlon(lat, lon)
  let m = u.to_mgrs(lat, precision=5)
  if u.zone != zone || !m.has_prefix(gzd) {
    abort("\{name}: expected zone \{zone} and \{gzd}, got \{u.zone} and \{m}")
  }
  let d = @utm.from_mgrs(m, centered=false)
  let back = @utm.to_latlon(u.zone, u.north, u.easting, u.northing)
  let utm_rt = g.inverse(lat, lon, back.lat, back.lon).s12
  // 1 m MGRS truncation: the SW corner must lie within sqrt(2) m.
  let mgrs_off = g.inverse(lat, lon, d.lat, d.lon).s12
  println(name)
  println("  UTM/UPS \{u.to_string()}  (k=\{u.k})")
  println("  MGRS    \{m}")
  println(
    "  UTM round trip \{utm_rt} m; MGRS 1 m cell corner offset \{mgrs_off} m",
  )
  if utm_rt > 1.0e-6 || mgrs_off > 1.4143 {
    abort("\{name}: residual too large")
  }
}
