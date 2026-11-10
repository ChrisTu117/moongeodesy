# magnetic

地磁模型 WMM2025：任意地点和日期的磁偏角、磁倾角、磁场强度及年变化，罗盘方位与真方位、格网方位互换。

Package `ChrisTu117/moongeodesy/magnetic`. Angles are degrees, heights are metres above the WGS-84 ellipsoid, intensities are nanotesla (nT), and rates are per year. The command-line front end is `src/examples/cli` (`mag`); `src/examples/s5_compass` is a worked field example.

The World Magnetic Model 2025 is produced by NOAA NCEI and the British Geological Survey for the US NGA and the UK DGC. It is a degree-12 spherical-harmonic model of the Earth's main field with a linear secular variation. The coefficients come from NOAA's `WMM.COF` (U.S. Government material, not subject to copyright protection; NOAA NCEI Geomagnetic Modeling Team and British Geological Survey, 2024, [doi:10.25921/aqfd-sd83](https://doi.org/10.25921/aqfd-sd83)). The code follows the equations in Section 1.2 of the technical report (Chulliat et al. 2025, [doi:10.25923/prbc-s316](https://doi.org/10.25923/prbc-s316)); it is not a port of NOAA's C library or of GeographicLib. `tools/gen_wmm.py` writes the coefficients and checks the SHA-256 of the official file.

## The model and its domain

模型与适用范围：2025.0 至 2030.0，椭球面下 1 km 至上 850 km。范围外抛出 `MagneticError`。

```moonbit nocheck
pub fn MagneticModel::wmm2025() -> MagneticModel
pub fn MagneticModel::field(Self, lat : Double, lon : Double, height : Double, year : Double) -> MagneticElements raise MagneticError
pub fn MagneticModel::declination(Self, lat : Double, lon : Double, height : Double, year : Double) -> Double raise MagneticError
pub fn decimal_year(year : Int, month : Int, day : Int) -> Double raise MagneticError
```

`field` raises `InvalidLatitude` for a latitude outside [-90, 90] or NaN, `InvalidLongitude` for a non-finite longitude, `HeightOutOfRange` outside [-1000, 850000] m, and `DateOutOfRange` outside [2025.0, 2030.0]. These are the limits in the model specification (MIL-PRF-89500B). The height is above the ellipsoid, not above mean sea level. The report (Section 1.2) puts the effect of that difference, the geoid height, at about 1 nT or less.

`decimal_year` uses the convention of NOAA's software: the year plus (day of year − 1) / (days in the year), so 1 January is `.0` and leap years have 366 days. It raises `InvalidDate` for a date that does not exist.

```mbt check
///|
test "declination in Beijing" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let t = @magnetic.decimal_year(2026, 11, 10)
  let e = wmm.field(39.9087, 116.3975, 0.0, t)
  // About 7.6 degrees west, drifting further west by about 0.04 deg/yr.
  assert_true(e.decl > -7.6 && e.decl < -7.5)
  assert_true(e.decl_dot < 0.0)
  assert_eq(e.zone(), @magnetic.Normal)
  let late = try wmm.field(39.9087, 116.3975, 0.0, 2031.0) catch {
    @magnetic.DateOutOfRange(_) => true
    _ => false
  } noraise {
    _ => false
  }
  assert_true(late)
}
```

## Magnetic elements

磁场七要素及其年变化。X 北、Y 东、Z 向下，垂直于 WGS-84 椭球面。

```moonbit nocheck
pub(all) struct MagneticElements {
  x : Double; y : Double; z : Double   // north, east, down (nT)
  h : Double; f : Double               // horizontal and total intensity (nT)
  incl : Double; decl : Double         // inclination (down positive), declination (east positive)
  gv : Double                          // polar grid variation, NaN between 55°S and 55°N
  xdot : Double; ydot : Double; zdot : Double; hdot : Double; fdot : Double  // nT/yr
  incl_dot : Double; decl_dot : Double // deg/yr
}
pub fn MagneticElements::zone(Self) -> CompassZone
pub(all) enum CompassZone { Normal; Caution; Blackout }
```

`decl` is the clockwise angle from true north to the horizontal field, in (−180, 180]. `incl` is positive when the field points down. `gv` is the WMM grid variation for the polar UPS grids: `decl - lon` at or north of 55°N and `decl + lon` at or south of 55°S, reduced to [−180, 180), and NaN in between, as in NOAA's test values. At an exact dip pole `h` is 0, the declination is meaningless, and the angle rates divide by zero.

`zone` applies the WMM2025 thresholds: `Blackout` where the horizontal intensity is below 2000 nT (the compass cannot be trusted), `Caution` below 6000 nT, otherwise `Normal`.

At the geographic poles the field is finite. North and east there refer to the meridian `lon`, so `x`, `y`, and `decl` rotate with `lon` while `h`, `z`, and `f` do not.

```mbt check
///|
test "near the north magnetic pole" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  // First official test point: 89°N 121°W, 28 km, 2025.0, H = 1504.3 nT.
  let e = wmm.field(89.0, -121.0, 28000.0, 2025.0)
  assert_true((e.h - 1504.298146).abs() < 1.0e-4)
  assert_eq(e.zone(), @magnetic.Blackout)
  assert_true(e.gv == e.decl + 121.0 || e.gv == e.decl + 121.0 - 360.0)
}
```

## Bearings and grids

罗盘方位、真方位与格网方位换算。格网磁偏角 = 磁偏角 − 子午线收敛角。

```moonbit nocheck
pub fn magnetic_bearing(true_bearing : Double, decl : Double) -> Double
pub fn true_bearing(magnetic_bearing : Double, decl : Double) -> Double
pub fn grivation(decl : Double, gamma : Double) -> Double
```

`magnetic_bearing` is what a compass shows for a true bearing, `true_bearing` is the inverse; both return [0, 360). `grivation` is the angle from grid north to magnetic north, `decl - gamma` reduced to [−180, 180), where `gamma` is the grid's meridian convergence, the clockwise angle from true north to grid north. For UTM it is `@utm.from_latlon(...).gamma`; for the north UPS grid it is `lon` and for the south UPS grid `-lon`. A grid bearing `b` read off a map is then `magnetic_bearing(b, grivation(decl, gamma))` on the compass.

```mbt check
///|
test "grid map to compass" {
  let decl = -7.5
  let gamma = -0.52 // UTM zone 50 convergence near Beijing's Western Hills
  let grid = 27.3
  let via_grid = @magnetic.magnetic_bearing(
    grid,
    @magnetic.grivation(decl, gamma),
  )
  let via_true = @magnetic.magnetic_bearing(grid + gamma, decl)
  assert_true((via_grid - via_true).abs() < 1.0e-12)
  assert_eq(
    @magnetic.true_bearing(@magnetic.magnetic_bearing(10.0, decl), decl),
    10.0,
  )
}
```

## Accuracy and tests

- **NOAA official test values, 100 points** (WMM2025COF.zip, double precision, 6 decimals, 2025.0–2029.5, 87°S–89°N, 0–98 km): Y and Z within 2.5e-6 nT (measured 2.2e-6), the file's rounding. X, H, and F within 1e-3 nT (measured 7.2e-4, 1.3e-8 of the field). The independent evaluation below reproduces our X rather than the file's, and the file agrees with both on the equator, so this small difference is in the file. Rates within 2e-6 nT/yr, D and I within their 2 printed decimals, their rates within 5.1e-7 deg/yr.
- **Report Table 6** (12 rounded values with grid variation): every element within its rounding.
- **Report Table 3b** (the worked example, 10 decimals): X, Y, Z, H, F, D, I and all rates within half a unit in the last place.
- **Report Table 3d** (grivation on UPS north, UPS south, UTM zones 10 and 11): matches to 1e-9 degrees, with the UTM convergence computed by this library's `utm` package.
- **Independent 40-digit evaluation, 119 points** (`tools/gen_wmm_independent.py`: exact Legendre polynomials, numerical gradient of the potential in ECEF, projection on the geodetic frame; includes both poles on three meridians, −1 km and 850 km, 2025.0 and 2030.0, and nine Chinese cities): X, Y, Z within 1e-10 nT (measured 5.1e-11) and rates within 1e-13 nT/yr.
- Property tests: finite and continuous at the poles, periodic in longitude, rates equal to finite differences in time, H, F, D, I consistent with X, Y, Z, and the domain errors.

These tests check that the code evaluates the model correctly. How well the model matches the real field is a separate matter: the report's error model (Section 3.4) gives a declination uncertainty of about 0.29° where the horizontal field is strong and 0.41° averaged over the globe, growing near the magnetic poles.
