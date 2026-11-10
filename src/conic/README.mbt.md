# conic

圆锥投影：兰勃特等角圆锥（LCC）和阿尔伯斯等积圆锥（Albers），单/双标准纬线，任意椭球，以及 EPSG 风格的带假原点格网。

Package `ChrisTu117/moongeodesy/conic`. Angles are degrees, lengths metres. The command-line front end is `src/examples/cli` (`lcc`, `albers`); `src/examples/s6_conic_area` compares areas on China grids.

The projections are ports of GeographicLib's `LambertConformalConic` and `AlbersEqualArea` (C. F. F. Karney, MIT/X11 License). They rewrite Snyder's formulas with divided differences, so they stay accurate for standard parallels close together, near the equator (the Mercator and cylindrical equal-area limits), and near a pole (the polar stereographic and azimuthal equal-area limits). They work on any ellipsoid with `f < 1`, including a sphere and prolate ellipsoids.

## Projections

投影本身。坐标相对于原点纬度和传入的中央经线。

```moonbit nocheck
pub fn LambertConformalConic::new(a : Double, f : Double, stdlat1 : Double, stdlat2 : Double, k1? : Double) -> LambertConformalConic raise ConicError
pub fn LambertConformalConic::new_1sp(a : Double, f : Double, stdlat : Double, k0? : Double) -> LambertConformalConic raise ConicError
pub fn LambertConformalConic::forward(Self, lon0 : Double, lat : Double, lon : Double) -> ConicProjected
pub fn LambertConformalConic::reverse(Self, lon0 : Double, x : Double, y : Double) -> ConicGeographic
pub fn LambertConformalConic::origin_latitude(Self) -> Double
pub fn LambertConformalConic::central_scale(Self) -> Double
pub fn LambertConformalConic::cone_constant(Self) -> Double

pub fn AlbersEqualArea::new(a : Double, f : Double, stdlat1 : Double, stdlat2 : Double, k1? : Double) -> AlbersEqualArea raise ConicError
pub fn AlbersEqualArea::new_1sp(a : Double, f : Double, stdlat : Double, k0? : Double) -> AlbersEqualArea raise ConicError
pub fn AlbersEqualArea::forward(Self, lon0 : Double, lat : Double, lon : Double) -> ConicProjected
pub fn AlbersEqualArea::reverse(Self, lon0 : Double, x : Double, y : Double) -> ConicGeographic

pub(all) struct ConicProjected { x : Double; y : Double; gamma : Double; k : Double }
pub(all) struct ConicGeographic { lat : Double; lon : Double; gamma : Double; k : Double }
```

`a` is the equatorial radius and `f` the flattening. With two standard parallels the scale on both is `k1` (default 1); with one, the scale on it is `k0`. `x` and `y` are relative to the point at `origin_latitude()` on the central meridian `lon0`: for the conformal projection that is where the scale is smallest, sin(lat0) = n. `gamma` is the meridian convergence, the bearing of grid north clockwise from true north, n (lon − lon0). For the conformal projection `k` is the scale in every direction; for the equal-area projection it is the scale along the parallel, and the scale along the meridian is 1 / k, so areas are unchanged.

Construction raises `ConicError`: `InvalidEllipsoid` unless `a` is finite and positive and `f < 1`; `InvalidScale` unless the scale is finite and positive; `InvalidStandardLatitude` outside [−90, 90]; `IncompatibleStandardLatitudes` for opposite poles, or (conformal) one parallel at a pole and the other not. Projecting a latitude outside [−90, 90] gives NaN. For the conformal projection the pole opposite the cone's apex is at infinity: `forward` returns very large coordinates for it (about 1e25 m on the 25°/47° cone), and `reverse` takes them back to that pole. Near a pole the equal-area projection squeezes the meridians so much that points within about 1e-6° (10 cm) of the pole fall on the pole's own map coordinates to rounding; `reverse` returns the pole for such inputs, so a pole round-trips exactly, a point 1e-6° away comes back within 4e-6°, and points 1e-4° or more away within 1e-7°. (PROJ also returns the pole there, and is 1e-4° off at 89.9999°.)

```mbt check
///|
test "scale is 1 on the standard parallels" {
  let lcc = @conic.LambertConformalConic::new(
    6378137.0,
    1.0 / 298.257222101,
    25.0,
    47.0,
  )
  assert_true((lcc.forward(105.0, 25.0, 120.0).k - 1.0).abs() < 1.0e-14)
  assert_true((lcc.forward(105.0, 47.0, 80.0).k - 1.0).abs() < 1.0e-14)
  assert_true(lcc.central_scale() < 1.0)
  let p = lcc.forward(105.0, 39.9087, 116.3975)
  let back = lcc.reverse(105.0, p.x, p.y)
  assert_true((back.lat - 39.9087).abs() < 1.0e-12)
  assert_true((back.lon - 116.3975).abs() < 1.0e-12)
}
```

## Grids with a false origin

EPSG 风格格网：假原点纬度、中央经线、假东/假北。

```moonbit nocheck
pub fn ConicGrid::lambert_2sp(a : Double, f : Double, lat_origin~ : Double, lon_origin~ : Double, stdlat1~ : Double, stdlat2~ : Double, false_easting? : Double, false_northing? : Double) -> ConicGrid raise ConicError
pub fn ConicGrid::lambert_1sp(a : Double, f : Double, lat_origin~ : Double, lon_origin~ : Double, k0? : Double, false_easting? : Double, false_northing? : Double) -> ConicGrid raise ConicError
pub fn ConicGrid::albers(a : Double, f : Double, lat_origin~ : Double, lon_origin~ : Double, stdlat1~ : Double, stdlat2~ : Double, false_easting? : Double, false_northing? : Double) -> ConicGrid raise ConicError
pub fn ConicGrid::forward(Self, lat : Double, lon : Double) -> GridPoint
pub fn ConicGrid::reverse(Self, easting : Double, northing : Double) -> ConicGeographic
pub(all) struct GridPoint { easting : Double; northing : Double; gamma : Double; k : Double }
```

These follow EPSG methods 9802 (Lambert Conic Conformal 2SP), 9801 (1SP, natural origin on the standard parallel with scale `k0`), and 9822 (Albers Equal Area): the false origin is (`lat_origin`, `lon_origin`), where the easting and northing are `false_easting` and `false_northing`. Lengths are metres; convert feet yourself (international foot 0.3048 m, US survey foot 1200/3937 m). `lat_origin` may be a pole (Belgian Lambert 72 uses 90°N).

China has no single official conic grid. Common choices are the Albers and Lambert projections with standard parallels 25°N and 47°N, central meridian 105°E (or 110°E), and the origin on the equator, on CGCS2000 or Krasovsky.

```mbt check
///|
test "EPSG Guidance Note 7-2: Jamaica National Grid" {
  // Clarke 1866, Lambert Conic Conformal (1SP), natural origin 18°N 77°W.
  let a = 6378206.4
  let jamaica = @conic.ConicGrid::lambert_1sp(
    a,
    1.0 - 6356583.8 / a,
    lat_origin=18.0,
    lon_origin=-77.0,
    false_easting=250000.0,
    false_northing=150000.0,
  )
  let lat = 17.0 + 55.0 / 60.0 + 55.80 / 3600.0
  let lon = -(76.0 + 56.0 / 60.0 + 37.26 / 3600.0)
  let p = jamaica.forward(lat, lon)
  assert_true((p.easting - 255966.58).abs() < 0.005)
  assert_true((p.northing - 142493.51).abs() < 0.005)
}
```

## Accuracy and tests

The references below do not use the formulation we port; `tools/gen_conic.py` regenerates them.

| Reference | What | Asserted bound (measured) |
| --- | --- | --- |
| 50-digit mpmath evaluation of the closed forms in IOGP Guidance Note 7-2 (= Snyder, USGS PP 1395, ch. 14–15) | 16 projections, 1592 points from 89°S to 89°N, up to 170° from the central meridian: China CGCS2000 and Krasovsky, US, Europe, Australia, southern hemisphere, parallels 1e-4° apart, polar, near the equator, one parallel | x, y < 2.5e-8 m (1.9e-8); convergence < 1e-12° (7e-14); k < 1e-13 relative (1.2e-15); reverse of the exact x, y < 2.5e-7 m (1.2e-7) |
| PROJ 9.5 (pyproj 3.7), Snyder formulas in C | the same points | x, y < 2.5e-7 m (1.9e-7); on the 1e-4° parallels < 1 cm (6 mm), where PROJ's own n = log(m1/m2)/log(t1/t2) loses digits; convergence < 5e-8° and scales < 1e-6, PROJ's numerical derivatives |
| IOGP GIGS test dataset 2.1, procedures 5102, 5103, 5109 | ED50 France EuroLambert, NTF Lambert zone II (gradians, Paris meridian), Belgian Lambert 72 (origin at the pole), Utah North in feet and US survey feet, Australian Albers | LCC < 1 mm and 1e-8°; Albers < 5 mm and 7e-8°, the files' rounding (PROJ shows the same 4.5 mm on GIGS-5109-09); GIGS's own tolerances are 3–5 cm |
| EPSG Guidance Note 7-2 worked examples | Texas South Central (2SP, US survey feet), Jamaica (1SP) | < 0.005 of the printed unit, reverse < 1e-7° |
| Snyder, PP 1395 Appendix A | Albers and Lambert, Clarke 1866 | < 0.1 m (Snyder carries about 8 digits) |

Property tests: both poles round-trip exactly for Albers (seven projections, including polar and cylindrical ones) and for Lambert; the Lambert scale measured with this library's geodesics matches `k` north, east, and north-east to 1e-9; the Albers area of a 30° × 50° box over China equals the ellipsoidal area to 1e-6 (measured 1.8e-9) and the parallel and meridian scales multiply to 1; the scale is 1 on the standard parallels and smallest at the origin latitude; standard parallel 0 gives Mercator and 90 gives this library's UPS to 1e-8 m; parallels 1e-9° apart reproduce 50-digit differences of a few micrometres to 2e-8 m.

GIGS data is © International Association of Oil & Gas Producers (IOGP), used under the [GIGS Terms of Use](https://gigs.iogp.org/TermsOfUse) v2.1, which allow redistribution with this acknowledgement and not for profit. It is used only in `src/conic/gigs_data_test.mbt`.
