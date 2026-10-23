///|
fn lc_max(errs : Array[Double]) -> Double {
  let mut m = 0.0
  for e in errs {
    if e > m || e.is_nan() {
      m = e
    }
  }
  m
}

///|
fn lc_abs_err(a : Double, b : Double) -> Double {
  if (a.is_nan() && b.is_nan()) || a == b {
    0.0
  } else {
    (a - b).abs()
  }
}

///|
fn lc_ang_err(a : Double, b : Double) -> Double {
  if a.is_nan() && b.is_nan() {
    return 0.0
  }
  let mut d = (a - b).abs()
  if d > 180.0 {
    d = 360.0 - d
  }
  d
}

///|
fn earth_for_tag(tag : Int32) -> @geodesic.Geocentric {
  match tag {
    0 => @geodesic.Geocentric::wgs84()
    1 =>
      @geodesic.Geocentric::new(6371000.0, 0.0) catch {
        _ => abort("sphere")
      }
    2 =>
      @geodesic.Geocentric::new(6378137.0, -1.0 / 298.257223563) catch {
        _ => abort("prolate")
      }
    3 => @geodesic.Geocentric::grs80()
    _ => abort("bad tag")
  }
}

///|
test "geocentric forward matches GeographicLib reference" {
  let errs = []
  for row in geocentric_forward {
    let (tag, lat, lon, h, rx, ry, rz) = row
    let p = earth_for_tag(tag).forward(lat, lon, h)
    errs.push(lc_abs_err(p.x, rx))
    errs.push(lc_abs_err(p.y, ry))
    errs.push(lc_abs_err(p.z, rz))
  }
  println("geocentric forward max |ΔXYZ| m = \{lc_max(errs)}")
  assert_true(lc_max(errs) < 1.0e-7)
}

///|
test "geocentric reverse round-trips forward" {
  let dxyz = []
  let dpos = []
  for row in geocentric_forward {
    let (tag, lat, lon, h, _, _, _) = row
    let earth = earth_for_tag(tag)
    let p = earth.forward(lat, lon, h)
    let g = earth.reverse(p.x, p.y, p.z)
    dpos.push(lc_ang_err(g.lat, lat) * 111319.5)
    dpos.push(
      lc_ang_err(g.lon, lon) *
      111319.5 *
      @math.cos(lat * 3.141592653589793 / 180.0).abs(),
    )
    dpos.push(lc_abs_err(g.h, h))
    let p2 = earth.forward(g.lat, g.lon, g.h)
    dxyz.push(lc_abs_err(p2.x, p.x))
    dxyz.push(lc_abs_err(p2.y, p.y))
    dxyz.push(lc_abs_err(p2.z, p.z))
  }
  println(
    "geocentric reverse max pos err m = \{lc_max(dpos)}, max |ΔXYZ| m = \{lc_max(dxyz)}",
  )
  assert_true(lc_max(dpos) < 1.0e-6)
  assert_true(lc_max(dxyz) < 1.0e-7)
}

///|
test "local cartesian forward matches GeographicLib reference" {
  let errs = []
  for row in local_cartesian_forward {
    let (lat0, lon0, h0, lat, lon, h, rx, ry, rz) = row
    let lc = @geodesic.LocalCartesian::wgs84(lat0, lon0, h0~)
    let p = lc.forward(lat, lon, h)
    errs.push(lc_abs_err(p.x, rx))
    errs.push(lc_abs_err(p.y, ry))
    errs.push(lc_abs_err(p.z, rz))
  }
  println("local cartesian forward max |Δxyz| m = \{lc_max(errs)}")
  assert_true(lc_max(errs) < 1.0e-7)
}

///|
test "local cartesian reverse round-trips forward" {
  let dpos = []
  let dxyz = []
  for row in local_cartesian_forward {
    let (lat0, lon0, h0, lat, lon, h, _, _, _) = row
    let lc = @geodesic.LocalCartesian::wgs84(lat0, lon0, h0~)
    let p = lc.forward(lat, lon, h)
    let g = lc.reverse(p.x, p.y, p.z)
    dpos.push(lc_ang_err(g.lat, lat) * 111319.5)
    dpos.push(
      lc_ang_err(g.lon, lon) *
      111319.5 *
      @math.cos(lat * 3.141592653589793 / 180.0).abs(),
    )
    dpos.push(lc_abs_err(g.h, h))
    let p2 = lc.forward(g.lat, g.lon, g.h)
    dxyz.push(lc_abs_err(p2.x, p.x))
    dxyz.push(lc_abs_err(p2.y, p.y))
    dxyz.push(lc_abs_err(p2.z, p.z))
  }
  println(
    "local cartesian reverse max pos err m = \{lc_max(dpos)}, max |Δxyz| m = \{lc_max(dxyz)}",
  )
  assert_true(lc_max(dpos) < 1.0e-6)
  assert_true(lc_max(dxyz) < 1.0e-7)
}

///|
test "local cartesian origin maps to zero" {
  let lc = @geodesic.LocalCartesian::wgs84(39.9087, 116.3975, h0=100.0)
  let p = lc.forward(39.9087, 116.3975, 100.0)
  assert_true(p.x.abs() < 1.0e-9)
  assert_true(p.y.abs() < 1.0e-9)
  assert_true(p.z.abs() < 1.0e-9)
  assert_eq(lc.latitude_origin(), 39.9087)
  assert_eq(lc.height_origin(), 100.0)
}

///|
test "local cartesian from_geodesic matches wgs84" {
  let g = @geodesic.Geodesic::wgs84()
  let a = @geodesic.LocalCartesian::from_geodesic(g, 48.0, 2.0)
  let b = @geodesic.LocalCartesian::wgs84(48.0, 2.0)
  let pa = a.forward(48.01, 2.01, 50.0)
  let pb = b.forward(48.01, 2.01, 50.0)
  assert_true(lc_abs_err(pa.x, pb.x) < 1.0e-9)
  assert_true(lc_abs_err(pa.y, pb.y) < 1.0e-9)
  assert_true(lc_abs_err(pa.z, pb.z) < 1.0e-9)
}

///|
test "geocentric rejects a bad ellipsoid" {
  let bad = try {
    ignore(@geodesic.Geocentric::new(0.0, 0.0))
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
