# cn

WGS-84、GCJ-02（火星坐标）与 BD-09（百度）互转。

Package `ChrisTu117/moongeodesy/cn`. Every angle is decimal degrees. These functions do not raise. The command-line front end is `src/examples/cli` (`coord`).

`out_of_china` is true outside the usual bounding box (longitude 72.004°..137.8347°, latitude 0.8293°..55.8271°) and on the land of other countries inside it, such as Seoul, Hanoi, and Ulaanbaatar. There GCJ-02 equals WGS-84. China, including Hong Kong, Macau, and Taiwan, and the sea inside the box get the GCJ-02 offset. The borders come from Natural Earth 1:10m admin-0 countries, China point of view (public domain), simplified to about 1 km near China (`tools/gen_china_boundary.py`); within about 2 km of a land border either answer is possible. `out_of_china_bbox` is the box-only test used by `coordtransform`, `gcoord`, and `coord-convert`. BD-09 is still applied everywhere. Tests cover the box edges, foreign capitals, Chinese cities including Hong Kong, Macau, and Taipei, and round trips next to every land border.

## Types and the China box

经纬度与中国范围判断。范围盒外，以及盒内其他国家的陆地上（首尔、河内、乌兰巴托等），GCJ-02 等于 WGS-84，不做偏移。

```moonbit nocheck
pub struct LatLon { lat : Double; lon : Double }
pub fn LatLon::new(lat : Double, lon : Double) -> LatLon
pub fn out_of_china(lat : Double, lon : Double) -> Bool
pub fn out_of_china_bbox(lat : Double, lon : Double) -> Bool
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
  let seoul = @cn.LatLon::new(37.5665, 126.978)
  assert_true(@cn.out_of_china(seoul.lat, seoul.lon))
  assert_false(@cn.out_of_china_bbox(seoul.lat, seoul.lon))
  assert_eq(@cn.wgs84_to_gcj02(seoul), seoul)
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

`wgs84_to_gcj02` is the published Krasovsky algorithm (`a = 6378245` m). Inside China it matches the usual open-source ciphers. Agreement with the official cipher is about 1 m. Where `out_of_china` is true it returns the input.

`gcj02_to_wgs84` and `bd09_to_gcj02` are fixed-point inverses. A round trip stays under 1e-9 degrees; for `gcj02_to_wgs84` this holds for every WGS-84 point inside China, including points next to the box edges and land borders. Next to a boundary the forward map is not one-to-one: when a GCJ-02 value has one source inside China and another just outside, the inside one is returned. `bd09_to_wgs84` is those two inverses in sequence. `wgs84_to_bd09` is the two forward steps.

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
