# utm

横轴墨卡托、极球面投影、UTM/UPS 与 MGRS。

Package `ChrisTu117/moongeodesy/utm`. Geographic angles are degrees. Easting and northing are metres. This package is WGS-84. A bad input raises `UtmError`.

## UTM and UPS

标准分带，含挪威与斯瓦尔巴特例。极区走 UPS。

```moonbit nocheck
pub fn standard_zone(lat : Double, lon : Double) -> Int
pub fn central_meridian(zone : Int) -> Double
pub fn from_latlon(lat : Double, lon : Double, zone? : Int) -> UtmCoord raise UtmError
pub fn to_latlon(
  zone : Int, north : Bool, easting : Double, northing : Double,
) -> Geographic raise UtmError
pub fn UtmCoord::to_string(Self) -> String
```

`UtmCoord` is `zone` (1..=60, or 0 for UPS), `north`, `easting`, `northing`, meridian convergence `gamma` (degrees), and scale `k`.

`standard_zone` returns 0 (UPS) when `lat >= 84` or `lat < -80`. Elsewhere it returns the UTM zone, including the Norway exception (band V, zone 32) and the Svalbard exceptions (band X, zones 31, 33, 35, 37). Longitude ±180 is zone 1.

`from_latlon` uses that zone unless `zone` is passed. `zone` must be 0..=60. Any other value raises `InvalidZone`. A latitude outside [-90, 90], including NaN, raises `InvalidLatitude`.

False offsets are included: UTM easting 500 km, southern UTM northing 10_000 km, UPS easting and northing 2_000 km.

Latitude exactly ±90 is inside the latitude check, but the UPS easting and northing are NaN. Use a point off the pole (UPS starts at 84°N and -80°S) when you need a finite grid coordinate.

`to_latlon` reverses that. A non-finite easting or northing, or one with absolute value >= 1e8 m, raises `OutOfRange`.

`Geographic` is `lat`, `lon`, `gamma` (degrees), and `k`.

```mbt check
///|
test "utm and ups" {
  let u = @utm.from_latlon(39.9087, 116.3975)
  assert_eq(u.zone, 50)
  assert_true(u.north)
  let g = @utm.to_latlon(u.zone, u.north, u.easting, u.northing)
  assert_true((g.lat - 39.9087).abs() < 1.0e-8)
  assert_true((g.lon - 116.3975).abs() < 1.0e-8)
  assert_eq(@utm.standard_zone(85.0, 10.0), 0)
  assert_eq(@utm.standard_zone(-80.5, 0.0), 0)
  assert_eq(@utm.central_meridian(50), 117.0)
  let ups = @utm.from_latlon(85.0, 10.0)
  assert_eq(ups.zone, 0)
  let bad_lat = try {
    ignore(@utm.from_latlon(91.0, 0.0))
    false
  } catch {
    _ => true
  }
  assert_true(bad_lat)
}
```

## MGRS

军用格网。精度 0..=11。字母 I 与 O 无效。

```moonbit nocheck
pub fn to_mgrs(lat : Double, lon : Double, precision? : Int) -> String raise UtmError
pub fn UtmCoord::to_mgrs(Self, lat : Double, precision? : Int) -> String raise UtmError
pub fn from_mgrs(input : String, centered? : Bool) -> MgrsPoint raise UtmError
```

`precision` is the number of digits on each of easting and northing. 5 is 1 m. 0 is the 100 km square only. 11 is 1 µm. The default is 5. A value outside 0..=11 raises `OutOfRange`.

`lat` selects the UTM band letter (C–X). UPS uses A, B, Y, and Z. Letters I and O are never used. A UTM easting outside 100 km..900 km, or a UPS point outside the polar grid, raises `OutOfRange`.

`from_mgrs` ignores spaces and tabs and ignores case. The default result is the centre of the square. `centered=false` returns the south-west corner. A bad zone, an odd digit count, more than 22 digits, or a letter outside the band or 100 km square raises `UtmError`.

`MgrsPoint` reports `zone`, `north`, `easting`, `northing`, `precision`, `lat`, and `lon`.

```mbt check
///|
test "mgrs" {
  let m = @utm.to_mgrs(39.9087, 116.3975)
  assert_eq(m, "50SMK4850217797")
  let p = @utm.from_mgrs("50S MK 48502 17797")
  assert_eq(p.zone, 50)
  assert_true(p.north)
  assert_true((p.lat - 39.9087).abs() < 1.0e-4)
  let odd = try {
    ignore(@utm.from_mgrs("50SMK485021779"))
    false
  } catch {
    _ => true
  }
  assert_true(odd)
}
```

## Projections

投影本体，不含假东偏与假北偏。

```moonbit nocheck
pub fn TransverseMercator::utm() -> TransverseMercator
pub fn TransverseMercator::new(a : Double, f : Double, k0 : Double) -> TransverseMercator
pub fn TransverseMercator::forward(
  Self, lon0 : Double, lat : Double, lon : Double,
) -> Projected
pub fn TransverseMercator::reverse(
  Self, lon0 : Double, x : Double, y : Double,
) -> Geographic

pub fn PolarStereographic::ups() -> PolarStereographic
pub fn PolarStereographic::new(a : Double, f : Double, k0 : Double) -> PolarStereographic
pub fn PolarStereographic::forward(
  Self, northp : Bool, lat : Double, lon : Double,
) -> Projected
pub fn PolarStereographic::reverse(
  Self, northp : Bool, x : Double, y : Double,
) -> Geographic
```

`TransverseMercator::utm()` is WGS-84 with `k0 = 0.9996`. `forward` takes the central meridian and the point, all in degrees. `Projected.x` and `Projected.y` are metres from that meridian and the equator. There is no false easting or northing.

The 6th-order Krüger series stays within a few nanometres for points within about 3900 km of the central meridian. Latitude ±90 is valid.

`PolarStereographic::ups()` is WGS-84 with `k0 = 0.994`. `northp=true` projects about the north pole. `x` and `y` are metres from the pole. Latitude exactly ±90 yields NaN easting and northing, including through `from_latlon`. A point a fraction of a degree off the pole is finite.

```mbt check
///|
test "transverse mercator and polar stereographic" {
  let tm = @utm.TransverseMercator::utm()
  let p = tm.forward(117.0, 39.9087, 116.3975)
  assert_true(p.x < 0.0)
  let back = tm.reverse(117.0, p.x, p.y)
  assert_true((back.lat - 39.9087).abs() < 1.0e-9)
  let ps = @utm.PolarStereographic::ups()
  let near = ps.forward(true, 89.9, 0.0)
  let back_ps = ps.reverse(true, near.x, near.y)
  assert_true((back_ps.lat - 89.9).abs() < 1.0e-9)
  assert_true(@utm.from_latlon(90.0, 0.0).easting.is_nan())
}
```
