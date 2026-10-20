///|
let ix_sphere : @geodesic.Intersect = @geodesic.Intersect::new(
  @geodesic.Geodesic::new(6378137.0, 0.0) catch {
    _ => abort("unreachable")
  },
) catch {
  _ => abort("unreachable")
}

///|
let ix_prolate : @geodesic.Intersect = @geodesic.Intersect::new(
  @geodesic.Geodesic::new(6378137.0, -@geodesic.wgs84_f) catch {
    _ => abort("unreachable")
  },
) catch {
  _ => abort("unreachable")
}

///|
let ix_grs80 : @geodesic.Intersect = @geodesic.Intersect::new(
  @geodesic.Geodesic::grs80(),
) catch {
  _ => abort("unreachable")
}

///|
fn ix_max(errs : Array[Double]) -> Double {
  let mut m = 0.0
  for e in errs {
    if e > m || e.is_nan() {
      m = e
    }
  }
  m
}

///|
fn ix_abs(a : Double, b : Double) -> Double {
  if a == b {
    0.0
  } else {
    (a - b).abs()
  }
}

///|
fn ix_pick(code : Double) -> @geodesic.Intersect {
  if code == 0.0 {
    @geodesic.Intersect::wgs84()
  } else if code == 1.0 {
    ix_sphere
  } else if code == 2.0 {
    ix_prolate
  } else {
    ix_grs80
  }
}

///|
fn ix_geod(code : Double) -> @geodesic.Geodesic {
  if code == 0.0 {
    @geodesic.Geodesic::wgs84()
  } else if code == 1.0 {
    @geodesic.Geodesic::new(6378137.0, 0.0) catch {
      _ => abort("unreachable")
    }
  } else if code == 2.0 {
    @geodesic.Geodesic::new(6378137.0, -@geodesic.wgs84_f) catch {
      _ => abort("unreachable")
    }
  } else {
    @geodesic.Geodesic::grs80()
  }
}

///|
fn ix_meet(
  g : @geodesic.Geodesic,
  lat1 : Double,
  lon1 : Double,
  azi1 : Double,
  s1 : Double,
  lat2 : Double,
  lon2 : Double,
  azi2 : Double,
  s2 : Double,
) -> Double {
  let p1 = g.direct(lat1, lon1, azi1, s1)
  let p2 = g.direct(lat2, lon2, azi2, s2)
  if p1.lat2.is_nan() || p2.lat2.is_nan() {
    return 1.0e9
  }
  g.distance(p1.lat2, p1.lon2, p2.lat2, p2.lon2)
}

///|
test "closest matches IntersectTool" {
  let dx = []
  let dy = []
  let sep = []
  for row in ix_closest {
    let ix = ix_pick(row[0])
    let r = ix.closest(row[1], row[2], row[3], row[4], row[5], row[6])
    let ex = ix_abs(r.x, row[7])
    let ey = ix_abs(r.y, row[8])
    // Two intersections can share the same L1 distance. Either is a
    // closest point. Require the distances to match, and the signed
    // position when this run picks the same intersection.
    let l1 = @geodesic.Intersect::dist(r.x, r.y)
    let l1r = @geodesic.Intersect::dist(row[7], row[8])
    assert_true(l1 <= l1r + 1.0e-6)
    assert_true((l1 - l1r).abs() < 1.0e-4)
    if ex < 1.0e-3 && ey < 1.0e-3 {
      dx.push(ex)
      dy.push(ey)
      assert_eq(r.c, row[9].to_int())
    } else {
      println(
        "intersect closest tie code=\{row[0]} got=\{r.x} \{r.y} ref=\{row[7]} \{row[8]}",
      )
    }
    sep.push(
      ix_meet(
        ix_geod(row[0]),
        row[1],
        row[2],
        row[3],
        r.x,
        row[4],
        row[5],
        row[6],
        r.y,
      ),
    )
  }
  println(
    "intersect closest n=\{ix_closest.length()} max |dx| m = \{ix_max(dx)}, max |dy| m = \{ix_max(dy)}, max separation m = \{ix_max(sep)}",
  )
  assert_true(ix_max(dx) < 1.0e-6)
  assert_true(ix_max(dy) < 1.0e-6)
  assert_true(ix_max(sep) < 1.0e-6)
}

///|
test "offset closest matches IntersectTool" {
  let dx = []
  let dy = []
  for row in ix_offset {
    let ix = ix_pick(row[0])
    let r = ix.closest(
      row[1],
      row[2],
      row[3],
      row[4],
      row[5],
      row[6],
      x0=row[7],
      y0=row[8],
    )
    dx.push(ix_abs(r.x, row[9]))
    dy.push(ix_abs(r.y, row[10]))
    assert_eq(r.c, row[11].to_int())
    let d = @geodesic.Intersect::dist(r.x, r.y, x0=row[7], y0=row[8])
    let dref = @geodesic.Intersect::dist(row[9], row[10], x0=row[7], y0=row[8])
    assert_true((d - dref).abs() < 1.0e-6)
  }
  println(
    "intersect offset n=\{ix_offset.length()} max |dx| m = \{ix_max(dx)}, max |dy| m = \{ix_max(dy)}",
  )
  assert_true(ix_max(dx) < 1.0e-6)
  assert_true(ix_max(dy) < 1.0e-6)
}

///|
test "next matches IntersectTool" {
  let dx = []
  let dy = []
  let sep = []
  for row in ix_next {
    let ix = ix_pick(row[0])
    let r = ix.next(row[1], row[2], row[3], row[4])
    dx.push(ix_abs(r.x, row[5]))
    dy.push(ix_abs(r.y, row[6]))
    assert_eq(r.c, row[7].to_int())
    sep.push(
      ix_meet(
        ix_geod(row[0]),
        row[1],
        row[2],
        row[3],
        r.x,
        row[1],
        row[2],
        row[4],
        r.y,
      ),
    )
    assert_true(r.x.abs() + r.y.abs() > 1.0)
  }
  println(
    "intersect next n=\{ix_next.length()} max |dx| m = \{ix_max(dx)}, max |dy| m = \{ix_max(dy)}, max separation m = \{ix_max(sep)}",
  )
  assert_true(ix_max(dx) < 1.0e-6)
  assert_true(ix_max(dy) < 1.0e-6)
  assert_true(ix_max(sep) < 1.0e-6)
}

///|
test "segment matches IntersectTool" {
  let dx = []
  let dy = []
  let sep = []
  for row in ix_segment {
    let ix = ix_pick(row[0])
    let g = ix_geod(row[0])
    let r = ix.segment(
      row[1],
      row[2],
      row[3],
      row[4],
      row[5],
      row[6],
      row[7],
      row[8],
    )
    let ex = ix_abs(r.x, row[9])
    let ey = ix_abs(r.y, row[10])
    let sx = g.distance(row[1], row[2], row[3], row[4])
    let sy = g.distance(row[5], row[6], row[7], row[8])
    let l1 = @geodesic.Intersect::dist(r.x, r.y, x0=sx / 2.0, y0=sy / 2.0)
    let l1r = @geodesic.Intersect::dist(
      row[9],
      row[10],
      x0=sx / 2.0,
      y0=sy / 2.0,
    )
    assert_true(l1 <= l1r + 1.0e-4)
    assert_true((l1 - l1r).abs() < 1.0e-3)
    let near = row[9].abs() < 1.0e-6 ||
      (row[9] - sx).abs() < 1.0e-6 ||
      row[10].abs() < 1.0e-6 ||
      (row[10] - sy).abs() < 1.0e-6
    if ex < 1.0e-3 && ey < 1.0e-3 {
      dx.push(ex)
      dy.push(ey)
      assert_eq(r.c, row[11].to_int())
      if !near {
        assert_eq(r.segmode, row[12].to_int())
      }
    } else {
      println(
        "intersect segment tie code=\{row[0]} got=\{r.x} \{r.y} k=\{r.segmode} ref=\{row[9]} \{row[10]} k=\{row[12]}",
      )
    }
    let kx = if r.x < 0.0 { -1 } else if r.x <= sx { 0 } else { 1 }
    let ky = if r.y < 0.0 { -1 } else if r.y <= sy { 0 } else { 1 }
    assert_eq(r.segmode, kx * 3 + ky)
    let px = g.inverse_line(row[1], row[2], row[3], row[4]).position(r.x)
    let py = g.inverse_line(row[5], row[6], row[7], row[8]).position(r.y)
    sep.push(g.distance(px.lat2, px.lon2, py.lat2, py.lon2))
  }
  println(
    "intersect segment n=\{ix_segment.length()} max |dx| m = \{ix_max(dx)}, max |dy| m = \{ix_max(dy)}, max separation m = \{ix_max(sep)}",
  )
  assert_true(ix_max(dx) < 1.0e-6)
  assert_true(ix_max(dy) < 1.0e-6)
  assert_true(ix_max(sep) < 1.0e-6)
}

///|
test "equatorial coincident segments match GeographicLib intersecttest" {
  // a = 180/pi, so one degree along the equator is one metre.
  // Cases and expected (x, y) are from GeographicLib tests/intersecttest.cpp.
  let ix = @geodesic.Intersect::new(
    @geodesic.Geodesic::new(180.0 / @geomath.pi, 0.0) catch {
      _ => abort("unreachable")
    },
  )
  let rows : Array[(Double, Double, Double, Double, Double, Double)] = [
    (0.0, 40.0, -20.0, -10.0, -5.0, 15.0),
    (0.0, 40.0, -10.0, 0.0, 0.0, 10.0),
    (0.0, 40.0, -8.0, 2.0, 1.0, 9.0),
    (0.0, 40.0, -2.0, 8.0, 4.0, 6.0),
    (0.0, 40.0, 0.0, 10.0, 5.0, 5.0),
    (0.0, 40.0, 2.0, 12.0, 7.0, 5.0),
    (0.0, 40.0, 15.0, 25.0, 20.0, 5.0),
    (0.0, 40.0, 30.0, 40.0, 35.0, 5.0),
    (0.0, 40.0, 32.0, 42.0, 36.0, 4.0),
    (0.0, 40.0, 38.0, 48.0, 39.0, 1.0),
    (0.0, 40.0, 40.0, 50.0, 40.0, 0.0),
    (0.0, 40.0, 50.0, 60.0, 45.0, -5.0),
    (40.0, 0.0, -20.0, -10.0, 45.0, 15.0),
    (40.0, 0.0, -10.0, 0.0, 40.0, 10.0),
    (40.0, 0.0, -8.0, 2.0, 39.0, 9.0),
    (40.0, 0.0, -2.0, 8.0, 36.0, 6.0),
    (40.0, 0.0, 0.0, 10.0, 35.0, 5.0),
    (40.0, 0.0, 2.0, 12.0, 33.0, 5.0),
    (40.0, 0.0, 15.0, 25.0, 20.0, 5.0),
    (40.0, 0.0, 30.0, 40.0, 5.0, 5.0),
    (40.0, 0.0, 32.0, 42.0, 4.0, 4.0),
    (40.0, 0.0, 38.0, 48.0, 1.0, 1.0),
    (40.0, 0.0, 40.0, 50.0, 0.0, 0.0),
    (40.0, 0.0, 50.0, 60.0, -5.0, -5.0),
  ]
  let mut worst = 0.0
  for row in rows {
    let (lonx1, lonx2, lony1, lony2, px, py) = row
    let a = ix.segment(0.0, lonx1, 0.0, lonx2, 0.0, lony1, 0.0, lony2)
    let b = ix.segment(0.0, lony1, 0.0, lony2, 0.0, lonx1, 0.0, lonx2)
    let e = ix_max([
      ix_abs(a.x, px),
      ix_abs(a.y, py),
      ix_abs(b.x, py),
      ix_abs(b.y, px),
    ])
    if e > worst {
      worst = e
    }
  }
  println("intersect equatorial max |d| m = \{worst}")
  assert_true(worst < 1.0e-6)
}

///|
test "intersection basic properties" {
  let ix = @geodesic.Intersect::wgs84()
  let g = @geodesic.Geodesic::wgs84()
  let cross = ix.closest(12.0, 40.0, 25.0, 12.0, 40.0, 110.0)
  assert_true(cross.x.abs() < 1.0e-6)
  assert_true(cross.y.abs() < 1.0e-6)
  assert_eq(cross.c, 0)
  let nxt = ix.next(12.0, 40.0, 25.0, 110.0)
  assert_true(nxt.x.abs() + nxt.y.abs() > 1.0e6)
  let back_x = g.direct(12.0, 40.0, 90.0, -100000.0)
  let fore_x = g.direct(12.0, 40.0, 90.0, 150000.0)
  let back_y = g.direct(12.0, 40.0, 10.0, -80000.0)
  let fore_y = g.direct(12.0, 40.0, 10.0, 120000.0)
  let s = ix.segment(
    back_x.lat2,
    back_x.lon2,
    fore_x.lat2,
    fore_x.lon2,
    back_y.lat2,
    back_y.lon2,
    fore_y.lat2,
    fore_y.lon2,
  )
  assert_true((s.x - 100000.0).abs() < 1.0e-4)
  assert_true((s.y - 80000.0).abs() < 1.0e-4)
  assert_eq(s.segmode, 0)
  assert_eq(s.c, 0)
  assert_eq(@geodesic.Intersect::dist(3.0, -4.0), 7.0)
  assert_eq(@geodesic.Intersect::dist(3.0, -4.0, x0=1.0, y0=-1.0), 5.0)
}

///|
test "too eccentric ellipsoid is rejected" {
  let g = @geodesic.Geodesic::new(6378137.0, -1.0 / 3.0)
  let bad = try {
    ignore(@geodesic.Intersect::new(g))
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
