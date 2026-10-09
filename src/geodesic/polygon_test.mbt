///|
fn near(a : Double, b : Double, rel_tol : Double) -> Bool {
  (a - b).abs() <= rel_tol * @cmp.maximum(1.0, b.abs())
}

///|
test "polygon area and perimeter match geographiclib" {
  let g = @geodesic.Geodesic::wgs84()
  for c in polygon_cases {
    let (name, pts, per, area, area_r) = c
    let p = g.polygon()
    for pt in pts {
      p.add_point(pt.0, pt.1)
    }
    let r = p.compute()
    assert_eq(r.num, pts.length())
    assert_true(near(r.perimeter, per, 1.0e-13), msg="perimeter \{name}")
    // area accuracy is limited by the ellipsoid area (5.1e14 m^2) for pole-enclosing cases
    assert_true(
      (r.area - area).abs() < 1.0e-14 * g.ellipsoid_area() + 1.0e-6,
      msg="area \{name}: \{r.area} vs \{area}",
    )
    let rr = p.compute(reverse=true, sign=false)
    assert_true(
      (rr.area - area_r).abs() < 1.0e-14 * g.ellipsoid_area() + 1.0e-6,
      msg="reverse area \{name}",
    )
  }
}

///|
test "add_edge, test_point, test_edge, polyline" {
  let g = @geodesic.Geodesic::wgs84()
  let p = g.polygon()
  p.add_point(30.0, 100.0)
  p.add_edge(90.0, 500000.0)
  p.add_edge(0.0, 300000.0)
  p.add_edge(270.0, 400000.0)
  let r = p.compute()
  assert_eq(r.num, edge_case.0)
  assert_true(near(r.perimeter, edge_case.1, 1.0e-13))
  assert_true(near(r.area, edge_case.2, 1.0e-12))
  let tp = p.test_point(31.0, 99.5)
  assert_eq(tp.num, test_point_case.0)
  assert_true(near(tp.perimeter, test_point_case.1, 1.0e-13))
  assert_true(near(tp.area, test_point_case.2, 1.0e-12))
  let te = p.test_edge(200.0, 100000.0)
  assert_eq(te.num, test_edge_case.0)
  assert_true(near(te.perimeter, test_edge_case.1, 1.0e-13))
  assert_true(near(te.area, test_edge_case.2, 1.0e-12))
  // test_* must not modify the polygon
  assert_eq(p.compute().num, 4)
  let l = g.polygon(polyline=true)
  for c in polygon_cases {
    if c.0 == "china_rough" {
      for pt in c.1 {
        l.add_point(pt.0, pt.1)
      }
    }
  }
  let lr = l.compute()
  assert_true(near(lr.perimeter, polyline_china_length, 1.0e-13))
  assert_true(lr.area.is_nan())
  l.clear()
  assert_eq(l.compute().num, 0)
}

///|
test "polygon_area convenience and degenerate polygons" {
  let g = @geodesic.Geodesic::wgs84()
  let r = g.polygon_area([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])
  assert_true(r.area > 0.0) // east, north, west: counter-clockwise → positive
  let cw = g.polygon_area([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
  assert_true((cw.area + r.area).abs() < 1.0e-3)
  assert_eq(g.polygon_area([]).area, 0.0)
  assert_eq(g.polygon_area([(1.0, 2.0)]).perimeter, 0.0)
}
