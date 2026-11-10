// Polygon area and perimeter against tools/gen_polygon_independent.py:
// closed forms (octant, lunes, hemispheres), L'Huilier's theorem on a
// sphere, and Green's theorem with mpmath quadrature, not geographiclib.

///|
fn polygon_indep_geod(code : Int32) -> @geodesic.Geodesic {
  if code == 0 {
    return @geodesic.Geodesic::wgs84()
  }
  let f = if code == 1 { -1.0 / 298.257223563 } else { 0.0 }
  let a = if code == 1 { 6378137.0 } else { 6371000.0 }
  @geodesic.Geodesic::new(a, f) catch {
    _ => abort("bad ellipsoid")
  }
}

///|
test "polygon area matches closed forms, L'Huilier and Green's theorem" {
  let mut d_closed = 0.0
  let mut d_sphere = 0.0
  let mut d_green = 0.0
  let mut d_per = 0.0
  let mut n = 0
  let mut n_closed = 0
  let mut n_sphere = 0
  let mut n_green = 0
  for c in polygon_indep_cases {
    let (name, code, pts, per, area, how) = c
    let g = polygon_indep_geod(code)
    let total = g.ellipsoid_area()
    let r = g.polygon_area(pts)
    n = n + 1
    assert_eq(r.num, pts.length())
    // Areas are defined modulo the ellipsoid's area; +A/2 and -A/2 are the
    // same hemisphere.
    let mut d = (r.area - area).abs() % total
    d = @cmp.minimum(d, total - d) / total
    if how == "closed form" {
      n_closed = n_closed + 1
      d_closed = indep_max(d_closed, d)
    } else if how == "L'Huilier" {
      n_sphere = n_sphere + 1
      d_sphere = indep_max(d_sphere, d)
    } else {
      n_green = n_green + 1
      d_green = indep_max(d_green, d)
    }
    d_per = indep_max(d_per, (r.perimeter - per).abs() / per)
    assert_true(d < 1.0e-14, msg="\{name}: \{r.area} vs \{area}")
  }
  println(
    "polygon vs independent: n=\{n} max |darea|/A closed forms=\{d_closed} sphere=\{d_sphere} Green=\{d_green}; max relative perimeter err=\{d_per}",
  )
  assert_eq(n, 41)
  assert_eq(n_closed, 26) // README: 13 closed forms on each of two ellipsoids
  assert_eq(n_sphere, 7) // README: 7 spherical triangles
  assert_eq(n_green, 8) // README: the 8 polygons of polygon_data_test.mbt
  assert_true(d_closed < 2.0e-16) // README: area < 2e-16 of the Earth's area
  assert_true(d_sphere < 2.0e-16)
  assert_true(d_green < 2.0e-16)
  assert_true(d_per < 5.0e-16) // README: perimeter < 5e-16 relative
}
