// World Magnetic Model (WMM2025): the Earth's main magnetic field as a
// degree-12 spherical-harmonic expansion with a linear secular variation.
// Equations follow Section 1.2 of Chulliat et al. (2025), "The US/UK World
// Magnetic Model for 2025-2030: Technical Report", NOAA NCEI,
// https://doi.org/10.25923/prbc-s316. The code is written from those
// equations; it is not a port of NOAA's C library or of GeographicLib.
//
// The Schmidt semi-normalized associated Legendre functions are evaluated as
// P(n, m)(sin p) = cos(p)^m * T(n, m)(sin p), where T is a polynomial. The
// easterly component needs P(n, m) / cos(p) = cos(p)^(m - 1) * T(n, m), which
// stays finite at the geographic poles, so no special pole branch is needed.

///|
/// Errors raised by `MagneticModel::field` and `decimal_year`.
#derive(Debug)
pub suberror MagneticError {
  /// Latitude outside [-90, 90] or not finite.
  InvalidLatitude(Double)
  /// Longitude not finite.
  InvalidLongitude(Double)
  /// Height (m above the WGS-84 ellipsoid) outside the model's range.
  HeightOutOfRange(Double)
  /// Decimal year outside the model's validity period.
  DateOutOfRange(Double)
  /// Calendar date that does not exist (year, month, day).
  InvalidDate(Int32, Int32, Int32)
}

///|
pub extend MagneticError with @debug.Debug::{to_repr}

///|
/// A spherical-harmonic main-field model with linear secular variation.
pub struct MagneticModel {
  /// Model name, e.g. "WMM2025".
  name : String
  /// Reference epoch t0 (decimal year).
  epoch : Double
  /// First decimal year of validity.
  valid_from : Double
  /// Last decimal year of validity.
  valid_to : Double
  /// Lowest valid height above the WGS-84 ellipsoid (m).
  min_height : Double
  /// Highest valid height above the WGS-84 ellipsoid (m).
  max_height : Double
  /// Maximum degree N.
  nmax : Int32
  /// Geomagnetic reference radius a (m).
  radius : Double
  priv g : FixedArray[Double]
  priv h : FixedArray[Double]
  priv gdot : FixedArray[Double]
  priv hdot : FixedArray[Double]
}

///|
/// The seven magnetic elements, their rates of change, and the WMM grid
/// variation at one place and time. Angles are degrees, intensities nT,
/// rates per year. Components are in the ellipsoidal frame: `x` north,
/// `y` east, `z` down (perpendicular to the WGS-84 ellipsoid).
#derive(Debug)
pub(all) struct MagneticElements {
  /// Northerly intensity X (nT).
  x : Double
  /// Easterly intensity Y (nT).
  y : Double
  /// Vertical intensity Z, positive down (nT).
  z : Double
  /// Horizontal intensity H (nT).
  h : Double
  /// Total intensity F (nT).
  f : Double
  /// Inclination (dip) I, positive down (degrees).
  incl : Double
  /// Declination (magnetic variation) D, clockwise from true north (degrees).
  decl : Double
  /// Grid variation GV for the polar UPS grids (degrees): D - lon at or
  /// north of 55°N, D + lon at or south of 55°S, both reduced to
  /// [-180, 180); NaN between 55°S and 55°N.
  gv : Double
  /// dX/dt (nT/yr).
  xdot : Double
  /// dY/dt (nT/yr).
  ydot : Double
  /// dZ/dt (nT/yr).
  zdot : Double
  /// dH/dt (nT/yr).
  hdot : Double
  /// dF/dt (nT/yr).
  fdot : Double
  /// dI/dt (degrees/yr).
  incl_dot : Double
  /// dD/dt (degrees/yr).
  decl_dot : Double
}

///|
pub extend MagneticElements with @debug.Debug::{to_repr}

///|
/// How far a compass can be trusted, from the horizontal intensity H.
#derive(Debug, Eq)
pub(all) enum CompassZone {
  /// H >= 6000 nT.
  Normal
  /// 2000 nT <= H < 6000 nT: the WMM2025 "Caution Zone".
  Caution
  /// H < 2000 nT: the WMM2025 "Blackout Zone"; the declination is unreliable.
  Blackout
}

///|
pub extend CompassZone with @debug.Debug::{to_repr}

///|
pub extend CompassZone with Eq::{equal, not_equal}

///|
/// WGS-84 semi-major axis (m), used to convert geodetic to geocentric
/// coordinates.
let wgs84_a : Double = 6378137.0

///|
/// WGS-84 flattening.
let wgs84_f : Double = 1.0 / 298.257223563

///|
/// The World Magnetic Model 2025 (NOAA NCEI and the British Geological
/// Survey), valid from 2025.0 to 2030.0 and from 1 km below to 850 km above
/// the WGS-84 ellipsoid.
pub fn MagneticModel::wmm2025() -> MagneticModel {
  {
    name: "WMM2025",
    epoch: wmm2025_epoch,
    valid_from: 2025.0,
    valid_to: 2030.0,
    min_height: -1000.0,
    max_height: 850000.0,
    nmax: 12,
    radius: 6371200.0,
    g: wmm2025_g,
    h: wmm2025_h,
    gdot: wmm2025_gdot,
    hdot: wmm2025_hdot,
  }
}

///|
/// Geocentric spherical coordinates (latitude in radians, radius in m) of a
/// point at geodetic latitude `lat` (degrees) and height `height` (m) on
/// WGS-84.
fn geocentric(lat : Double, height : Double) -> (Double, Double) {
  let (sphi, cphi) = @geomath.sincosd(lat)
  let e2 = wgs84_f * (2.0 - wgs84_f)
  let rc = wgs84_a / (1.0 - e2 * sphi * sphi).sqrt()
  let p = (rc + height) * cphi
  let z = (rc * (1.0 - e2) + height) * sphi
  (@math.atan2(z, p), @geomath.hypot(p, z))
}

///|
/// Gauss coefficient arrays at decimal year `t`: g(t) = g + (t - t0) gdot.
fn MagneticModel::coefficients_at(
  self : MagneticModel,
  t : Double,
) -> (FixedArray[Double], FixedArray[Double]) {
  let dt = t - self.epoch
  let n = self.g.length()
  let gt = FixedArray::make(n, 0.0)
  let ht = FixedArray::make(n, 0.0)
  for k in 0..<n {
    gt[k] = self.g[k] + dt * self.gdot[k]
    ht[k] = self.h[k] + dt * self.hdot[k]
  }
  (gt, ht)
}

///|
/// Polynomial parts T(n, m)(x) and dT/dx of the Schmidt semi-normalized
/// associated Legendre functions, P(n, m)(x) = (1 - x^2)^(m/2) T(n, m)(x),
/// indexed by n * (n + 1) / 2 + m.
fn legendre_poly(
  nmax : Int32,
  x : Double,
) -> (FixedArray[Double], FixedArray[Double]) {
  let size = (nmax + 1) * (nmax + 2) / 2
  let t = FixedArray::make(size, 0.0)
  let dt = FixedArray::make(size, 0.0)
  for m in 0..<=nmax {
    let kmm = m * (m + 1) / 2 + m
    // Sectoral term: T(0,0) = T(1,1) = 1; T(m,m) = sqrt((2m-1)/(2m)) T(m-1,m-1).
    t[kmm] = if m <= 1 {
      1.0
    } else {
      let km1 = (m - 1) * m / 2 + (m - 1)
      ((2 * m - 1).to_double() / (2 * m).to_double()).sqrt() * t[km1]
    }
    // Recurrence in n, starting from T(m-1, m) = 0.
    for n in (m + 1)..<=nmax {
      let k = n * (n + 1) / 2 + m
      let k1 = (n - 1) * n / 2 + m
      let nn = n.to_double()
      let mm = m.to_double()
      let a = (2 * n - 1).to_double() / (nn * nn - mm * mm).sqrt()
      let b = if n - 2 >= m {
        (((nn - 1.0) * (nn - 1.0) - mm * mm) / (nn * nn - mm * mm)).sqrt()
      } else {
        0.0
      }
      let (t2, dt2) = if n - 2 >= m {
        let k2 = (n - 2) * (n - 1) / 2 + m
        (t[k2], dt[k2])
      } else {
        (0.0, 0.0)
      }
      t[k] = a * x * t[k1] - b * t2
      dt[k] = a * (t[k1] + x * dt[k1]) - b * dt2
    }
  }
  (t, dt)
}

///|
/// Field components in the geocentric spherical frame (north, east, down),
/// from Gauss coefficients `gc` and `hc`, at geocentric latitude `phip`
/// (radians), radius `r` (m), and longitude `lon` (degrees).
fn MagneticModel::sum_field(
  self : MagneticModel,
  gc : FixedArray[Double],
  hc : FixedArray[Double],
  t : FixedArray[Double],
  dt : FixedArray[Double],
  sp : Double,
  cp : Double,
  ratio : Double,
  lon : Double,
) -> (Double, Double, Double) {
  let nmax = self.nmax
  // cpow[k] = cos(phi')^k, with cpow[0] = 1 even at the poles.
  let cpow = FixedArray::make(nmax + 2, 1.0)
  for k in 1..<(nmax + 2) {
    cpow[k] = cpow[k - 1] * cp
  }
  let cosm = FixedArray::make(nmax + 1, 1.0)
  let sinm = FixedArray::make(nmax + 1, 0.0)
  for m in 1..<=nmax {
    let (s, c) = @geomath.sincosd(m.to_double() * lon)
    sinm[m] = s
    cosm[m] = c
  }
  let mut bx = 0.0
  let mut by = 0.0
  let mut bz = 0.0
  let mut rn = ratio * ratio // (a/r)^(n+2) for n = 0
  for n in 1..<=nmax {
    rn = rn * ratio
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    for m in 0..<=n {
      let k = n * (n + 1) / 2 + m
      let cs = gc[k] * cosm[m] + hc[k] * sinm[m]
      // P(n,m) and dP(n,m)/dphi'
      let pnm = cpow[m] * t[k]
      let dpnm = if m == 0 {
        cpow[1] * dt[k]
      } else {
        -m.to_double() * sp * cpow[m - 1] * t[k] + cpow[m + 1] * dt[k]
      }
      sx = sx + cs * dpnm
      sz = sz + cs * pnm
      if m > 0 {
        // m * P(n,m) / cos(phi') = m * cos(phi')^(m-1) * T(n,m)
        sy = sy +
          m.to_double() *
          (gc[k] * sinm[m] - hc[k] * cosm[m]) *
          cpow[m - 1] *
          t[k]
      }
    }
    bx = bx - rn * sx
    by = by + rn * sy
    bz = bz - (n + 1).to_double() * rn * sz
  }
  (bx, by, bz)
}

///|
/// Reduce an angle in degrees to [-180, 180).
fn reduce180(x : Double) -> Double {
  let y = @geomath.remainder(x, 360.0)
  if y == 180.0 {
    -180.0
  } else {
    y
  }
}

///|
/// The magnetic elements at geodetic latitude `lat` and longitude `lon`
/// (degrees), `height` metres above the WGS-84 ellipsoid, and decimal year
/// `year` (see `decimal_year`).
///
/// Raises `InvalidLatitude`, `InvalidLongitude`, `HeightOutOfRange`, or
/// `DateOutOfRange` outside the model's domain.
pub fn MagneticModel::field(
  self : MagneticModel,
  lat : Double,
  lon : Double,
  height : Double,
  year : Double,
) -> MagneticElements raise MagneticError {
  if !(lat >= -90.0 && lat <= 90.0) {
    raise InvalidLatitude(lat)
  }
  if !@geomath.is_finite(lon) {
    raise InvalidLongitude(lon)
  }
  if !(height >= self.min_height && height <= self.max_height) {
    raise HeightOutOfRange(height)
  }
  if !(year >= self.valid_from && year <= self.valid_to) {
    raise DateOutOfRange(year)
  }
  let (phip, r) = geocentric(lat, height)
  let sp = @math.sin(phip)
  let cp = @math.cos(phip)
  let (t, dt) = legendre_poly(self.nmax, sp)
  let ratio = self.radius / r
  let (gc, hc) = self.coefficients_at(year)
  let (xp, yp, zp) = self.sum_field(gc, hc, t, dt, sp, cp, ratio, lon)
  let (xpd, ypd, zpd) = self.sum_field(
    self.gdot,
    self.hdot,
    t,
    dt,
    sp,
    cp,
    ratio,
    lon,
  )
  // Rotate from the geocentric to the ellipsoidal frame.
  let psi = phip - lat * @geomath.rad_per_deg
  let spsi = @math.sin(psi)
  let cpsi = @math.cos(psi)
  let x = xp * cpsi - zp * spsi
  let z = xp * spsi + zp * cpsi
  let y = yp
  let xdot = xpd * cpsi - zpd * spsi
  let zdot = xpd * spsi + zpd * cpsi
  let ydot = ypd
  let h = @geomath.hypot(x, y)
  let f = @geomath.hypot(h, z)
  let incl = @geomath.atan2d(z, h)
  let decl = @geomath.atan2d(y, x)
  let hdot = (x * xdot + y * ydot) / h
  let fdot = (x * xdot + y * ydot + z * zdot) / f
  let incl_dot = (h * zdot - z * hdot) / (f * f) * @geomath.deg_per_rad
  let decl_dot = (x * ydot - y * xdot) / (h * h) * @geomath.deg_per_rad
  let gv = if lat >= 55.0 {
    reduce180(decl - lon)
  } else if lat <= -55.0 {
    reduce180(decl + lon)
  } else {
    @double.not_a_number
  }
  {
    x,
    y,
    z,
    h,
    f,
    incl,
    decl,
    gv,
    xdot,
    ydot,
    zdot,
    hdot,
    fdot,
    incl_dot,
    decl_dot,
  }
}

///|
/// Declination D (degrees) at a place and time; shorthand for
/// `field(...).decl`.
pub fn MagneticModel::declination(
  self : MagneticModel,
  lat : Double,
  lon : Double,
  height : Double,
  year : Double,
) -> Double raise MagneticError {
  self.field(lat, lon, height, year).decl
}

///|
/// Whether a compass is usable here, by the WMM2025 Blackout Zone (H below
/// 2000 nT) and Caution Zone (H below 6000 nT) thresholds.
pub fn MagneticElements::zone(self : MagneticElements) -> CompassZone {
  if self.h < 2000.0 {
    Blackout
  } else if self.h < 6000.0 {
    Caution
  } else {
    Normal
  }
}

///|
/// Grivation (grid magnetic angle), the clockwise angle from grid north to
/// magnetic north, in degrees reduced to [-180, 180): `decl - gamma`, where
/// `gamma` is the meridian convergence of the grid, the clockwise angle from
/// true north to grid north (for UTM, `@utm.from_latlon(...).gamma`).
pub fn grivation(decl : Double, gamma : Double) -> Double {
  reduce180(decl - gamma)
}

///|
/// Magnetic bearing (degrees, [0, 360)) shown by a compass for a true
/// bearing `true_bearing`, given the declination `decl`.
pub fn magnetic_bearing(true_bearing : Double, decl : Double) -> Double {
  let b = @geomath.remainder(true_bearing - decl, 360.0)
  if b < 0.0 {
    b + 360.0
  } else {
    b
  }
}

///|
/// True bearing (degrees, [0, 360)) for a compass reading `magnetic_bearing`,
/// given the declination `decl`.
pub fn true_bearing(magnetic_bearing : Double, decl : Double) -> Double {
  let b = @geomath.remainder(magnetic_bearing + decl, 360.0)
  if b < 0.0 {
    b + 360.0
  } else {
    b
  }
}

///|
fn is_leap(year : Int32) -> Bool {
  (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

///|
/// Decimal year of a calendar date at 00:00, as in the WMM software:
/// year + (day of year - 1) / (days in year). Raises `InvalidDate` for a
/// date that does not exist.
pub fn decimal_year(
  year : Int32,
  month : Int32,
  day : Int32,
) -> Double raise MagneticError {
  let days = [
    31,
    if is_leap(year) {
      29
    } else {
      28
    },
    31,
    30,
    31,
    30,
    31,
    31,
    30,
    31,
    30,
    31,
  ]
  if month < 1 || month > 12 || day < 1 || day > days[month - 1] {
    raise InvalidDate(year, month, day)
  }
  let mut doy = day
  for i in 0..<(month - 1) {
    doy = doy + days[i]
  }
  let len = if is_leap(year) { 366.0 } else { 365.0 }
  year.to_double() + (doy - 1).to_double() / len
}
