///|
test "transverse Mercator forward/reverse vs PROJ" {
  let tm = @utm.TransverseMercator::utm()
  let mut max_xy = 0.0
  let mut max_xy_near = 0.0
  let mut max_gamma = 0.0
  let mut max_k = 0.0
  let mut max_rt = 0.0
  for c in tm_cases {
    let (lat, lon, lon0, x, y, gamma, k) = c
    let p = tm.forward(lon0, lat, lon)
    let e = @cmp.maximum((p.x - x).abs(), (p.y - y).abs())
    if e > max_xy {
      max_xy = e
    }
    if (lon - lon0).abs() <= 3.0 && e > max_xy_near {
      max_xy_near = e
    }
    max_gamma = @cmp.maximum(max_gamma, (p.gamma - gamma).abs())
    max_k = @cmp.maximum(max_k, (p.k - k).abs())
    let g = tm.reverse(lon0, p.x, p.y)
    let mut dlon = (g.lon - lon).abs()
    if dlon > 180.0 {
      dlon = 360.0 - dlon
    }
    let rt = @cmp.maximum((g.lat - lat).abs(), dlon)
    max_rt = @cmp.maximum(max_rt, rt)
    // reverse gives the same convergence and scale
    assert_true((g.gamma - p.gamma).abs() < 1.0e-9)
    assert_true((g.k - p.k).abs() < 1.0e-12)
  }
  println(
    "TM vs PROJ: max |dxy| \{max_xy} m (|dlon|<=3: \{max_xy_near} m), max dgamma \{max_gamma} deg, max dk \{max_k}, round trip \{max_rt} deg",
  )
  assert_true(max_xy_near < 1.0e-8) // 10 nm
  assert_true(max_xy < 1.0e-8)
  assert_true(max_gamma < 1.0e-8)
  assert_true(max_k < 1.0e-9) // PROJ's factors are numerical derivatives
  assert_true(max_rt < 1.0e-12)
}
