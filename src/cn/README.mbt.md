# cn

WGS-84、GCJ-02（火星坐标）与 BD-09（百度）互转。

Package `ChrisTu117/moongeodesy/cn`. Every angle is decimal degrees. These functions do not raise. The command-line front end is `src/examples/cli` (`coord`).

`out_of_china` is the usual bounding box, not the national border: longitude 72.004°..137.8347°, latitude 0.8293°..55.8271°. Outside that box, GCJ-02 equals WGS-84. BD-09 is still applied outside the box.

## Types and the China box

经纬度与中国范围盒。盒外 GCJ-02 等于 WGS-84。

```moonbit nocheck
pub struct LatLon { lat : Double; lon : Double }
pub fn LatLon::new(lat : Double, lon : Double) -> LatLon
pub fn out_of_china(lat : Double, lon : Double) -> Bool
```

`LatLon::new` does not range-check. `Show` prints `LatLon(lat, lon)`.

```mbt check
///|
test "out of china" {
  let paris = @cn.LatLon::new(48.8566, 2.3522)
  assert_true(@cn.out_of_china(paris.lat, paris.lon))
  assert_eq(@cn.wgs84_to_gcj02(paris), paris)
  assert_eq(@cn.gcj02_to_wgs84(paris), paris)
  assert_false(@cn.out_of_china(39.9087, 116.3975))
  let bd = @cn.wgs84_to_bd09(paris)
  assert_true((bd.lat - paris.lat).abs() > 1.0e-6)
}
```

## WGS-84, GCJ-02, and BD-09

正向偏移与迭代精确逆。往返误差小于 1e-9°。

```moonbit nocheck
pub fn wgs84_to_gcj02(LatLon) -> LatLon
pub fn gcj02_to_wgs84(LatLon) -> LatLon
pub fn gcj02_to_wgs84_approx(LatLon) -> LatLon
pub fn gcj02_to_bd09(LatLon) -> LatLon
pub fn bd09_to_gcj02(LatLon) -> LatLon
pub fn bd09_to_gcj02_approx(LatLon) -> LatLon
pub fn wgs84_to_bd09(LatLon) -> LatLon
pub fn bd09_to_wgs84(LatLon) -> LatLon
```

`wgs84_to_gcj02` is the published Krasovsky algorithm (`a = 6378245` m). Inside the box it matches the usual open-source ciphers. Agreement with the official cipher is about 1 m. Outside the box it returns the input.

`gcj02_to_wgs84` and `bd09_to_gcj02` are fixed-point inverses. A round trip stays under 1e-9 degrees. `bd09_to_wgs84` is those two inverses in sequence. `wgs84_to_bd09` is the two forward steps.

`gcj02_to_wgs84_approx` is one subtraction of the GCJ-02 offset. The error can reach about 1–2 m. `bd09_to_gcj02_approx` is Baidu's published one-step inverse. Its error can reach about 1e-6 degrees. Prefer the exact inverses when a later geodesic or grid step needs the GPS point.

```mbt check
///|
test "china coordinate round trip" {
  let w = @cn.LatLon::new(39.9087, 116.3975)
  let g = @cn.wgs84_to_gcj02(w)
  let b = @cn.gcj02_to_bd09(g)
  assert_true((g.lat - w.lat).abs() > 1.0e-4)
  let back = @cn.bd09_to_wgs84(b)
  assert_true((back.lat - w.lat).abs() < 1.0e-9)
  assert_true((back.lon - w.lon).abs() < 1.0e-9)
  let approx = @cn.gcj02_to_wgs84_approx(g)
  assert_true((approx.lat - w.lat).abs() < 5.0e-5)
  assert_true((approx.lon - w.lon).abs() < 5.0e-5)
}
```
