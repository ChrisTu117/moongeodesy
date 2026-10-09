// Polar stereographic projection, following GeographicLib::PolarStereographic
// (Charles Karney, MIT/X11 License).

///|
/// Polar stereographic projection on an ellipsoid.
pub struct PolarStereographic {
  a : Double
  f : Double
  k0 : Double
  e2 : Double
  es : Double
  e2m : Double
  c : Double
}

///|
pub fn PolarStereographic::new(
  a : Double,
  f : Double,
  k0 : Double,
) -> PolarStereographic {
  let e2 = f * (2.0 - f)
  let es = (if f < 0.0 { -1.0 } else { 1.0 }) * e2.abs().sqrt()
  let e2m = 1.0 - e2
  let c = (1.0 - f) * @math.exp(eatanhe(1.0, es))
  { a, f, k0, e2, es, e2m, c, }
}

///|
let ups_ps : PolarStereographic = PolarStereographic::new(
  6378137.0,
  1.0 / 298.257223563,
  0.994,
)

///|
/// The UPS projection (WGS-84, k0 = 0.994).
pub fn PolarStereographic::ups() -> PolarStereographic {
  ups_ps
}

///|
/// Project (lat, lon) about the north (`northp=true`) or south pole.
pub fn PolarStereographic::forward(
  self : PolarStereographic,
  northp : Bool,
  lat_in : Double,
  lon : Double,
) -> Projected {
  let lat = @geomath.lat_fix(lat_in) * (if northp { 1.0 } else { -1.0 })
  let (sphi, cphi) = @geomath.sincosd(lat)
  let tau = sphi / cphi
  let secphi = @geomath.hypot(1.0, tau)
  let taup = taupf(tau, self.es)
  let mut rho = @geomath.hypot(1.0, taup) + taup.abs()
  rho = if taup >= 0.0 {
    if lat != 90.0 {
      1.0 / rho
    } else {
      0.0
    }
  } else {
    rho
  }
  rho = rho * (2.0 * self.k0 * self.a / self.c)
  let k = if lat != 90.0 {
    rho / self.a * secphi * (self.e2m + self.e2 / (secphi * secphi)).sqrt()
  } else {
    self.k0
  }
  let (slon, clon) = @geomath.sincosd(lon)
  let x = rho * slon
  let y = (if northp { -rho } else { rho }) * clon
  let gamma = @geomath.ang_normalize(if northp { lon } else { -lon })
  { x, y, gamma, k, }
}

///|
/// Inverse polar stereographic projection.
pub fn PolarStereographic::reverse(
  self : PolarStereographic,
  northp : Bool,
  x : Double,
  y : Double,
) -> Geographic {
  let rho = @geomath.hypot(x, y)
  let t = if rho != 0.0 {
    rho / (2.0 * self.k0 * self.a / self.c)
  } else {
    @geomath.dbl_epsilon * @geomath.dbl_epsilon
  }
  let taup = (1.0 / t - t) / 2.0
  let tau = tauf(taup, self.es)
  let secphi = @geomath.hypot(1.0, tau)
  let k = if rho != 0.0 {
    rho / self.a * secphi * (self.e2m + self.e2 / (secphi * secphi)).sqrt()
  } else {
    self.k0
  }
  let lat = (if northp { 1.0 } else { -1.0 }) * @geomath.atan2d(tau, 1.0)
  let lon = @geomath.atan2d(x, if northp { -y } else { y })
  let gamma = @geomath.ang_normalize(if northp { lon } else { -lon })
  { lat, lon, gamma, k, }
}
