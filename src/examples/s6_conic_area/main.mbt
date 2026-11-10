///|
/// Scenario 6: national statistics measure land area on a CGCS2000 Albers
/// grid (standard parallels 25°N and 47°N, central meridian 105°E) and
/// print maps on a Lambert grid with the same parallels. How good are the
/// areas, and how large is the Lambert scale error across China? Many web
/// GIS tools measure area on Web Mercator instead; compare that too.
/// Run: `moon run src/examples/s6_conic_area`
fn main {
  let a = 6378137.0
  let f = 1.0 / 298.257222101 // CGCS2000
  let g = @geodesic.Geodesic::cgcs2000()
  let albers = @conic.ConicGrid::albers(
    a,
    f,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  let lambert = @conic.ConicGrid::lambert_2sp(
    a,
    f,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  println("Area of a 0.5° x 0.5° cell (km^2), edges sampled every 0.005°:")
  println(
    "  place      ellipsoid      Albers  (rel. err)       Lambert  (rel. err)   Web Mercator  (rel. err)",
  )
  let places = [
    ("Sanya", 18.0, 109.25),
    ("Shenzhen", 22.5, 114.0),
    ("Shanghai", 31.0, 121.25),
    ("Beijing", 39.75, 116.25),
    ("Urumqi", 43.75, 87.5),
    ("Harbin", 45.75, 126.5),
    ("Mohe", 52.75, 122.25),
  ]
  for p in places {
    let ring = cell(p.1, p.2, 0.5, 100)
    let exact = g.polygon_area(ring).area.abs()
    let aa = planar_area(ring, (lat, lon) => {
      let q = albers.forward(lat, lon)
      (q.easting, q.northing)
    })
    let la = planar_area(ring, (lat, lon) => {
      let q = lambert.forward(lat, lon)
      (q.easting, q.northing)
    })
    let wm = planar_area(ring, web_mercator)
    println(
      "  \{p.0}: \{exact / 1.0e6}  \{aa / 1.0e6} (\{(aa - exact) / exact})  \{la / 1.0e6} (\{(la - exact) / exact})  \{wm / 1.0e6} (\{(wm - exact) / exact})",
    )
    if ((aa - exact) / exact).abs() > 1.0e-6 {
      abort("Albers area is off")
    }
  }
  println("")
  println("Lambert 25°N/47°N scale factor and meridian convergence:")
  for lat in [18.0, 25.0, 30.0, 36.0, 42.0, 47.0, 53.0] {
    let at105 = lambert.forward(lat, 105.0)
    let at75 = lambert.forward(lat, 75.0)
    let at135 = lambert.forward(lat, 135.0)
    println(
      "  \{lat}°N: k = \{at105.k} (\{(at105.k - 1.0) * 1000.0} m per km); convergence at 75°E \{at75.gamma}°, at 135°E \{at135.gamma}°",
    )
    if (at105.k - at75.k).abs() > 1.0e-15 {
      abort("Lambert scale should not depend on longitude")
    }
  }
  // Grid round trip.
  let p = albers.forward(39.9087, 116.3975)
  let back = albers.reverse(p.easting, p.northing)
  println("")
  println(
    "Beijing on the Albers grid: E \{p.easting} m, N \{p.northing} m; back to \{back.lat}, \{back.lon}",
  )
  if (back.lat - 39.9087).abs() > 1.0e-12 ||
    (back.lon - 116.3975).abs() > 1.0e-12 {
    abort("round trip failed")
  }
}

///|
/// A closed lat/lon cell with `n` points per side, counter-clockwise.
fn cell(
  lat0 : Double,
  lon0 : Double,
  size : Double,
  n : Int32,
) -> Array[(Double, Double)] {
  let ring = []
  let step = size / n.to_double()
  for i in 0..<n {
    ring.push((lat0, lon0 + step * i.to_double()))
  }
  for i in 0..<n {
    ring.push((lat0 + step * i.to_double(), lon0 + size))
  }
  for i in 0..<n {
    ring.push((lat0 + size, lon0 + size - step * i.to_double()))
  }
  for i in 0..<n {
    ring.push((lat0 + size - step * i.to_double(), lon0))
  }
  ring
}

///|
/// Shoelace area of the projected ring.
fn planar_area(
  ring : Array[(Double, Double)],
  proj : (Double, Double) -> (Double, Double),
) -> Double {
  let m = ring.length()
  let mut s = 0.0
  for i in 0..<m {
    let (x1, y1) = proj(ring[i].0, ring[i].1)
    let (x2, y2) = proj(ring[(i + 1) % m].0, ring[(i + 1) % m].1)
    s = s + (x1 * y2 - x2 * y1)
  }
  s.abs() / 2.0
}

///|
/// Web Mercator (EPSG:3857) easting and northing.
fn web_mercator(lat : Double, lon : Double) -> (Double, Double) {
  let r = 6378137.0
  let rad = 3.141592653589793 / 180.0
  (
    r * lon * rad,
    r * @math.ln(@math.tan(3.141592653589793 / 4.0 + lat * rad / 2.0)),
  )
}
