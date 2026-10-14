///|
/// Scenario 1: a delivery app gets a shop from Gaode/Tencent (GCJ-02) and a
/// customer from Baidu Maps (BD-09). It needs the true ground distance.
/// Run: `moon run src/examples/s1_cn_distance`
fn main {
  let g = @geodesic.Geodesic::wgs84()
  // Shop near Shenzhen Civic Center as reported by Gaode (GCJ-02).
  let shop_gcj = @cn.LatLon::new(22.543099, 114.057868)
  // Customer near Shenzhen Bay Park as reported by Baidu (BD-09).
  let cust_bd = @cn.LatLon::new(22.522806, 113.951238)
  let shop = @cn.gcj02_to_wgs84(shop_gcj)
  let cust = @cn.bd09_to_wgs84(cust_bd)
  println("shop     GCJ-02 \{shop_gcj} -> WGS-84 \{shop}")
  println("customer BD-09  \{cust_bd} -> WGS-84 \{cust}")
  let right = g.inverse(shop.lat, shop.lon, cust.lat, cust.lon)
  let wrong = g.inverse(shop_gcj.lat, shop_gcj.lon, cust_bd.lat, cust_bd.lon)
  println("geodesic distance after conversion : \{right.s12} m")
  println("distance mixing raw GCJ-02 + BD-09 : \{wrong.s12} m")
  println(
    "error avoided                      : \{(wrong.s12 - right.s12).abs()} m",
  )
  // How far each raw coordinate sits from the true WGS-84 position.
  let off_shop = g.inverse(shop.lat, shop.lon, shop_gcj.lat, shop_gcj.lon).s12
  let off_cust = g.inverse(cust.lat, cust.lon, cust_bd.lat, cust_bd.lon).s12
  println(
    "GCJ-02 offset of shop: \{off_shop} m; BD-09 offset of customer: \{off_cust} m",
  )
  // Round trip check: WGS-84 -> GCJ-02 -> WGS-84 (iterative inverse).
  let back = @cn.gcj02_to_wgs84(@cn.wgs84_to_gcj02(shop))
  let rt = g.inverse(shop.lat, shop.lon, back.lat, back.lon).s12
  println("round-trip residual: \{rt} m")
}
