// MoonBit port of GeographicLib's Rhumb and RhumbLine (Charles Karney,
// MIT/X11 License, https://geographiclib.sourceforge.io/).
// Series formulation (order 6), the default used by RhumbSolve:
// C. F. F. Karney, "The area of rhumb polygons", Stud. Geophys. Geod. 68(3-4),
// 99-120 (2024), https://doi.org/10.1007/s11200-024-0709-z
// The series is accurate for |f| <= 1/150. WGS-84, GRS80, and CGCS2000
// are inside that range. This is a transcription, not a binding.

///|
let aux_phi : Int32 = 0

///|
let aux_beta : Int32 = 1

///|
let aux_mu : Int32 = 3

///|
let aux_chi : Int32 = 4

///|
let rm_coeff : FixedArray[Double] = [1.0 / 256.0, 1.0 / 64.0, 1.0 / 4.0, 1.0]

///|
let c2_coeff : FixedArray[Double] = [
  4.0 / 1287.0,
  4.0 / 693.0,
  4.0 / 315.0,
  4.0 / 105.0,
  4.0 / 15.0,
  -1.0 / 3.0,
  1.0,
]

///|
// Fourier coefficients of the rhumb-area integral, order 6.
// GeographicLib Rhumb.cpp, GEOGRAPHICLIB_RHUMBAREA_ORDER == 6.
let area_coeff : FixedArray[Double] = [
  138734126.0 / 638512875.0,
  -102614.0 / 467775.0,
  596.0 / 2025.0,
  -398.0 / 945.0,
  22.0 / 45.0,
  -1.0 / 3.0,
  17749373.0 / 425675250.0,
  -24562.0 / 155925.0,
  1543.0 / 4725.0,
  -118.0 / 315.0,
  1.0 / 5.0,
  1882432.0 / 8513505.0,
  -38068.0 / 155925.0,
  152.0 / 945.0,
  -17.0 / 315.0,
  268864.0 / 2027025.0,
  -752.0 / 10395.0,
  5.0 / 252.0,
  62464.0 / 2027025.0,
  -101.0 / 17325.0,
  11537.0 / 4054050.0,
]

///|
/// Error raised for an invalid ellipsoid.
pub suberror RhumbError {
  InvalidEllipsoid(String)
} derive(Debug)

///|
pub extend RhumbError with @debug.Debug::{to_repr}

///|
/// An ellipsoid of revolution on which rhumb-line problems are solved.
pub struct Rhumb {
  /// Equatorial radius in metres.
  a : Double
  /// Flattening.
  f : Double
  n : Double
  /// Rectifying radius (m).
  rm : Double
  /// Authalic radius squared, times 1 degree in radians. Used for area.
  c2 : Double
  /// Fourier coefficients for auxiliary-latitude conversions.
  c : FixedArray[Double]
  /// Fourier coefficients for the area integral.
  p_p : FixedArray[Double]
}

///|
/// A rhumb line of constant azimuth. Build one with `Rhumb::line`, then
/// call `position` for points along it.
pub struct RhumbLine {
  rh : Rhumb
  lat1 : Double
  lon1 : Double
  azi12 : Double
  salp : Double
  calp : Double
  /// Rectifying latitude of point 1, in degrees.
  mu1 : Double
  phi1 : AuxAngle
  chi1 : AuxAngle
}

///|
/// Result of the inverse rhumb problem.
pub(all) struct RhumbInverseResult {
  /// Latitude of point 1 (degrees). NaN when the input was outside [-90, 90].
  lat1 : Double
  /// Longitude of point 1 (degrees), reduced to [-180, 180].
  lon1 : Double
  /// Latitude of point 2 (degrees).
  lat2 : Double
  /// Longitude of point 2 (degrees), reduced to [-180, 180].
  lon2 : Double
  /// Rhumb distance (m).
  s12 : Double
  /// Constant azimuth (degrees, clockwise from north).
  azi12 : Double
  /// Area between the rhumb line and the equator, S12 (m²).
  area12 : Double
} derive(Debug)

///|
pub extend RhumbInverseResult with @debug.Debug::{to_repr}

///|
/// Result of the direct rhumb problem.
pub(all) struct RhumbDirectResult {
  /// Latitude of point 1 (degrees). NaN when the input was outside [-90, 90].
  lat1 : Double
  /// Longitude of point 1 (degrees).
  lon1 : Double
  /// Constant azimuth (degrees, clockwise from north).
  azi12 : Double
  /// Rhumb distance (m).
  s12 : Double
  /// Latitude of point 2 (degrees).
  lat2 : Double
  /// Longitude of point 2 (degrees). NaN when the line crosses a pole.
  lon2 : Double
  /// Area between the rhumb line and the equator, S12 (m²). NaN when the
  /// line crosses a pole.
  area12 : Double
} derive(Debug)

///|
pub extend RhumbDirectResult with @debug.Debug::{to_repr}

///|
fn Rhumb::fill_coeffs(self : Rhumb) -> Unit {
  let n2 = self.n * self.n
  for k in 0..<(aux_number * aux_number) {
    let auxout = k / aux_number
    let auxin = k % aux_number
    if auxout != auxin {
      let mut o = aux_ptrs[k]
      let mut d = self.n
      if auxin <= aux_rectifying && auxout <= aux_rectifying {
        for l in 0..<aux_lmax {
          let m = (aux_lmax - l - 1) / 2
          self.c[aux_lmax * k + l] = d * @geomath.polyval(m, aux_coeffs, o, n2)
          o = o + m + 1
          d = d * self.n
        }
      } else {
        for l in 0..<aux_lmax {
          let m = aux_lmax - l - 1
          self.c[aux_lmax * k + l] = d *
            @geomath.polyval(m, aux_coeffs, o, self.n)
          o = o + m + 1
          d = d * self.n
        }
      }
    }
  }
}

///|
fn Rhumb::fill_area(self : Rhumb) -> Unit {
  let mut d = 1.0
  let mut o = 0
  for l in 0..<aux_lmax {
    let m = aux_lmax - l - 1
    d = d * self.n
    self.p_p[l] = d * @geomath.polyval(m, area_coeff, o, self.n)
    o = o + m + 1
  }
}

///|
/// Construct an ellipsoid with equatorial radius `a` (metres) and
/// flattening `f`. Raises `InvalidEllipsoid` when `a` or the polar
/// semi-axis `a * (1 - f)` is not finite and positive.
/// The series assume `|f|` up to about 1/150.
/// WGS-84, GRS80, and CGCS2000 are inside that range.
pub fn Rhumb::new(a : Double, f : Double) -> Rhumb raise RhumbError {
  let b = a * (1.0 - f)
  if !(@geomath.is_finite(a) && a > 0.0) {
    raise InvalidEllipsoid("Equatorial radius is not positive")
  }
  if !(@geomath.is_finite(b) && b > 0.0) {
    raise InvalidEllipsoid("Polar semi-axis is not positive")
  }
  let n = f / (2.0 - f)
  let rm = (a + b) / 2.0 * @geomath.polyval(aux_lmax / 2, rm_coeff, 0, n * n)
  let c2 = a *
    (a + b) /
    2.0 *
    @geomath.polyval(aux_lmax, c2_coeff, 0, n) *
    @geomath.rad_per_deg
  let rh : Rhumb = {
    a,
    f,
    n,
    rm,
    c2,
    c: FixedArray::make(aux_lmax * aux_number * aux_number, 0.0),
    p_p: FixedArray::make(aux_lmax, 0.0),
  }
  rh.fill_coeffs()
  rh.fill_area()
  rh
}

///|
let wgs84_rhumb : Rhumb = Rhumb::new(wgs84_a, wgs84_f) catch {
  _ => abort("unreachable")
}

///|
/// The WGS-84 ellipsoid.
pub fn Rhumb::wgs84() -> Rhumb {
  wgs84_rhumb
}

///|
let grs80_rhumb : Rhumb = Rhumb::new(6378137.0, 1.0 / 298.257222101) catch {
  _ => abort("unreachable")
}

///|
/// The GRS80 ellipsoid (also used by CGCS2000).
pub fn Rhumb::grs80() -> Rhumb {
  grs80_rhumb
}

///|
/// The CGCS2000 ellipsoid (a = 6378137 m, f = 1/298.257222101).
pub fn Rhumb::cgcs2000() -> Rhumb {
  grs80_rhumb
}

///|
/// Equatorial radius `a` in metres.
pub fn Rhumb::equatorial_radius(self : Rhumb) -> Double {
  self.a
}

///|
/// Flattening `f` (dimensionless).
pub fn Rhumb::flattening(self : Rhumb) -> Double {
  self.f
}

///|
/// Total area of the ellipsoid in m².
pub fn Rhumb::ellipsoid_area(self : Rhumb) -> Double {
  // c2 includes a factor of one degree in radians, and 720 * c2 = 4 * pi * R².
  720.0 * self.c2
}

///|
fn Rhumb::convert(
  self : Rhumb,
  auxin : Int32,
  auxout : Int32,
  zeta : AuxAngle,
) -> AuxAngle {
  if auxin == auxout {
    return zeta
  }
  let k = aux_number * auxout + auxin
  let z = zeta.normalized()
  let d = aux_clenshaw(true, z.y, z.x, self.c, aux_lmax * k, aux_lmax)
  z.add(AuxAngle::from_radians(d))
}

///|
fn Rhumb::dconvert(
  self : Rhumb,
  auxin : Int32,
  auxout : Int32,
  zeta1 : AuxAngle,
  zeta2 : AuxAngle,
) -> Double {
  if auxin == auxout {
    return 1.0
  }
  let k = aux_number * auxout + auxin
  let z1 = zeta1.normalized()
  let z2 = zeta2.normalized()
  1.0 +
  aux_dclenshaw(
    true,
    z2.to_radians() - z1.to_radians(),
    z1.y,
    z1.x,
    z2.y,
    z2.x,
    self.c,
    aux_lmax * k,
    aux_lmax,
  )
}

///|
fn Rhumb::mean_sin_xi(self : Rhumb, chix : AuxAngle, chiy : AuxAngle) -> Double {
  let phix = self.convert(aux_chi, aux_phi, chix)
  let phiy = self.convert(aux_chi, aux_phi, chiy)
  let betax = self.convert(aux_phi, aux_beta, phix).normalized()
  let betay = self.convert(aux_phi, aux_beta, phiy).normalized()
  let dp = aux_dclenshaw(
    false,
    betay.to_radians() - betax.to_radians(),
    betax.y,
    betax.x,
    betay.y,
    betay.x,
    self.p_p,
    0,
    aux_lmax,
  )
  let tx = chix.tan()
  let ty = chiy.tan()
  let dbeta_dpsi = self.dconvert(aux_chi, aux_beta, chix, chiy) /
    aux_dlam(tx, ty)
  aux_dp0_dpsi(tx, ty) + dp * dbeta_dpsi
}

///|
fn nan_inverse(
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> RhumbInverseResult {
  {
    lat1,
    lon1,
    lat2,
    lon2,
    s12: @double.not_a_number,
    azi12: @double.not_a_number,
    area12: @double.not_a_number,
  }
}

///|
/// Solve the inverse rhumb problem: distance and constant azimuth between
/// two points, plus the area between that rhumb line and the equator.
///
/// Latitudes and longitudes are degrees. `s12` is metres. `area12` is
/// square metres. Azimuth is clockwise from north, in [-180, 180].
/// The shortest rhumb line is returned (it wraps at most halfway around
/// the earth). On opposite meridians the east-going line is chosen.
/// A latitude outside [-90, 90] comes back as NaN. No error is raised.
/// Latitude ±90 is a valid pole.
pub fn Rhumb::inverse(
  self : Rhumb,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> RhumbInverseResult {
  let lat1f = @geomath.lat_fix(lat1)
  let lat2f = @geomath.lat_fix(lat2)
  let lon1n = @geomath.ang_normalize(lon1)
  let lon2n = @geomath.ang_normalize(lon2)
  if lat1f.is_nan() || lat2f.is_nan() {
    return nan_inverse(lat1f, lon1n, lat2f, lon2n)
  }
  let phi1 = AuxAngle::from_degrees(lat1f)
  let phi2 = AuxAngle::from_degrees(lat2f)
  let chi1 = self.convert(aux_phi, aux_chi, phi1)
  let chi2 = self.convert(aux_phi, aux_chi, phi2)
  let (lon12, _) = @geomath.ang_diff(lon1, lon2)
  let lam12 = lon12 * @geomath.rad_per_deg
  let psi1 = chi1.lam()
  let psi2 = chi2.lam()
  let psi12 = psi2 - psi1
  let azi12 = @geomath.atan2d(lam12, psi12)
  let s12 = if psi1.is_inf() || psi2.is_inf() {
    let mu1 = self.convert(aux_phi, aux_mu, phi1)
    let mu2 = self.convert(aux_phi, aux_mu, phi2)
    (mu2.to_radians() - mu1.to_radians()).abs() * self.rm
  } else {
    let h = @geomath.hypot(lam12, psi12)
    let dmudpsi = self.dconvert(aux_chi, aux_mu, chi1, chi2) /
      aux_dlam(chi1.tan(), chi2.tan())
    h * dmudpsi * self.rm
  }
  let area12 = self.c2 * lon12 * self.mean_sin_xi(chi1, chi2)
  { lat1: lat1f, lon1: lon1n, lat2: lat2f, lon2: lon2n, s12, azi12, area12, }
}

///|
/// Rhumb distance between two points, in metres.
pub fn Rhumb::distance(
  self : Rhumb,
  lat1 : Double,
  lon1 : Double,
  lat2 : Double,
  lon2 : Double,
) -> Double {
  self.inverse(lat1, lon1, lat2, lon2).s12
}

///|
/// Start a rhumb line at (lat1, lon1) with constant azimuth `azi12`
/// (degrees, clockwise from north).
pub fn Rhumb::line(
  self : Rhumb,
  lat1 : Double,
  lon1 : Double,
  azi12 : Double,
) -> RhumbLine {
  let lat1f = @geomath.lat_fix(lat1)
  let azi = @geomath.ang_normalize(azi12)
  let (salp, calp) = @geomath.sincosd(azi)
  let phi1 = AuxAngle::from_degrees(lat1f)
  let mu1 = self.convert(aux_phi, aux_mu, phi1).to_degrees()
  let chi1 = self.convert(aux_phi, aux_chi, phi1)
  { rh: self, lat1: lat1f, lon1, azi12: azi, salp, calp, mu1, phi1, chi1, }
}

///|
fn nan_direct(
  lat1 : Double,
  lon1 : Double,
  azi12 : Double,
  s12 : Double,
) -> RhumbDirectResult {
  {
    lat1,
    lon1,
    azi12,
    s12,
    lat2: @double.not_a_number,
    lon2: @double.not_a_number,
    area12: @double.not_a_number,
  }
}

///|
/// Point a distance `s12` metres along this rhumb line.
/// `s12` may be negative (travel the other way).
/// Longitude is reduced to [-180, 180] unless `long_unroll` is true.
/// If the line crosses a pole, `lon2` and `area12` are NaN.
pub fn RhumbLine::position(
  self : RhumbLine,
  s12 : Double,
  long_unroll~ : Bool = false,
) -> RhumbDirectResult {
  let lon1o = if long_unroll {
    self.lon1
  } else {
    @geomath.ang_normalize(self.lon1)
  }
  if self.lat1.is_nan() {
    return nan_direct(self.lat1, lon1o, self.azi12, s12)
  }
  let r12 = s12 / (self.rh.rm * @geomath.rad_per_deg)
  let mu2 = self.mu1 + r12 * self.calp
  let (lat2, lon2, area12) = if mu2.abs() <= 90.0 {
    let phi2 = self.rh.convert(aux_mu, aux_phi, AuxAngle::from_degrees(mu2))
    let chi2 = self.rh.convert(aux_phi, aux_chi, phi2)
    let dmudpsi = self.rh.dconvert(aux_chi, aux_mu, self.chi1, chi2) /
      aux_dlam(self.chi1.tan(), chi2.tan())
    let lon2x = r12 * self.salp / dmudpsi
    let area = self.rh.c2 * lon2x * self.rh.mean_sin_xi(self.chi1, chi2)
    let lon2 = if long_unroll {
      self.lon1 + lon2x
    } else {
      @geomath.ang_normalize(@geomath.ang_normalize(self.lon1) + lon2x)
    }
    (phi2.to_degrees(), lon2, area)
  } else {
    let mut mu2n = @geomath.ang_normalize(mu2)
    if mu2n.abs() > 90.0 {
      mu2n = @geomath.ang_normalize(180.0 - mu2n)
    }
    let lat2 = self.rh
      .convert(aux_mu, aux_phi, AuxAngle::from_degrees(mu2n))
      .to_degrees()
    (lat2, @double.not_a_number, @double.not_a_number)
  }
  { lat1: self.lat1, lon1: lon1o, azi12: self.azi12, s12, lat2, lon2, area12, }
}

///|
/// Solve the direct rhumb problem: travel `s12` metres from (lat1, lon1)
/// at constant azimuth `azi12`. See `RhumbLine::position` for poles and
/// `long_unroll`.
pub fn Rhumb::direct(
  self : Rhumb,
  lat1 : Double,
  lon1 : Double,
  azi12 : Double,
  s12 : Double,
  long_unroll~ : Bool = false,
) -> RhumbDirectResult {
  self.line(lat1, lon1, azi12).position(s12, long_unroll~)
}

///|
/// Latitude of the start point (degrees).
pub fn RhumbLine::latitude(self : RhumbLine) -> Double {
  self.lat1
}

///|
/// Longitude of the start point (degrees), as given to `Rhumb::line`.
pub fn RhumbLine::longitude(self : RhumbLine) -> Double {
  self.lon1
}

///|
/// Constant azimuth (degrees, clockwise from north), reduced to [-180, 180].
pub fn RhumbLine::azimuth(self : RhumbLine) -> Double {
  self.azi12
}
