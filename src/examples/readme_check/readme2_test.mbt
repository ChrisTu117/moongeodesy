// Compiles and runs the code snippets from README.md (keep in sync).

///|
test "readme zh" {
  let g = @geodesic.Geodesic::wgs84()
  let r = g.inverse(39.9087, 116.3975, 31.2304, 121.4737) // 北京 → 上海：r.s12 ≈ 1066.68 km
  let gcj = @cn.wgs84_to_gcj02(@cn.LatLon::new(39.9087, 116.3975)) // GPS → 高德
  let bd = @cn.gcj02_to_bd09(gcj) // 高德 → 百度
  let mgrs = @utm.to_mgrs(39.9087, 116.3975) // "50SMK4850217797"
  let rh = @geodesic.Rhumb::wgs84()
  let h = rh.inverse(51.4778, -0.0015, 35.6762, 139.6503) // 伦敦 → 东京，恒向线
  let ix = @geodesic.Intersect::wgs84()
  let p = ix.closest(0.0, 0.0, 45.0, 45.0, 10.0, 135.0) // 两条测地线的最近交点
  let lc = @geodesic.LocalCartesian::wgs84(39.9087, 116.3975) // 以天安门为原点的局部 ENU
  let enu = lc.forward(39.9097, 116.3985, 50.0) // enu.x 东, enu.y 北, enu.z 天
  ignore((r, bd, mgrs, h, p, enu))
}
