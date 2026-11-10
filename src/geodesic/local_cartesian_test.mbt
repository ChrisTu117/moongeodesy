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
fn lc_of(
  tag : Int32,
  lat0 : Double,
  lon0 : Double,
  h0 : Double,
) -> @geodesic.LocalCartesian {
  let e = earth_for_tag(tag)
  @geodesic.LocalCartesian::new(
    e.equatorial_radius(),
    e.flattening(),
    lat0,
    lon0,
    h0~,
  ) catch {
    _ => abort("bad ellipsoid")
  }
}

///|
/// Horizontal position error (m) between two geographic points.
fn lc_pos_err(
  lat : Double,
  lon : Double,
  lat2 : Double,
  lon2 : Double,
) -> Double {
  let dn = lc_ang_err(lat, lat2) * 111319.5
  let de = if lat.abs() == 90.0 {
    0.0 // longitude is arbitrary at a pole
  } else {
    lc_ang_err(lon, lon2) *
    111319.5 *
    @math.cos(lat * 3.141592653589793 / 180.0).abs()
  }
  if dn > de {
    dn
  } else {
    de
  }
}

///|
test "geocentric forward matches 50-digit values and PROJ" {
  // 660 points on WGS-84, a sphere, a prolate ellipsoid, and GRS80: both
  // poles, the equator, the antimeridian, heights from -5 km to 35 786 km.
  let errs = []
  for row in geocentric_forward {
    let p = earth_for_tag(row[0].to_int()).forward(row[1], row[2], row[3])
    errs.push(lc_abs_err(p.x, row[4]))
    errs.push(lc_abs_err(p.y, row[5]))
    errs.push(lc_abs_err(p.z, row[6]))
  }
  let perrs = []
  for row in geocentric_proj {
    let p = earth_for_tag(row[0].to_int()).forward(row[1], row[2], row[3])
    perrs.push(lc_abs_err(p.x, row[4]))
    perrs.push(lc_abs_err(p.y, row[5]))
    perrs.push(lc_abs_err(p.z, row[6]))
  }
  println(
    "geocentric forward: max |dXYZ| vs 50-digit = \{lc_max(errs)} m, vs PROJ = \{lc_max(perrs)} m",
  )
  assert_true(lc_max(errs) < 2.0e-8)
  assert_true(lc_max(perrs) < 2.0e-8)
}

///|
test "geocentric reverse matches 50-digit values" {
  let dpos = []
  let dh = []
  for row in geocentric_reverse {
    let g = earth_for_tag(row[0].to_int()).reverse(row[1], row[2], row[3])
    dpos.push(lc_pos_err(row[4], row[5], g.lat, g.lon))
    dh.push(lc_abs_err(g.h, row[6]))
  }
  println(
    "geocentric reverse: max position error = \{lc_max(dpos)} m, max |dh| = \{lc_max(dh)} m",
  )
  assert_true(lc_max(dpos) < 2.0e-8)
  assert_true(lc_max(dh) < 2.0e-8)
}

///|
test "local cartesian forward matches 50-digit values and PROJ topocentric" {
  let errs = []
  for row in local_cartesian_forward {
    let lc = lc_of(row[0].to_int(), row[1], row[2], row[3])
    let p = lc.forward(row[4], row[5], row[6])
    errs.push(lc_abs_err(p.x, row[7]))
    errs.push(lc_abs_err(p.y, row[8]))
    errs.push(lc_abs_err(p.z, row[9]))
  }
  let perrs = []
  for row in local_cartesian_proj {
    let lc = lc_of(row[0].to_int(), row[1], row[2], row[3])
    let p = lc.forward(row[4], row[5], row[6])
    perrs.push(lc_abs_err(p.x, row[7]))
    perrs.push(lc_abs_err(p.y, row[8]))
    perrs.push(lc_abs_err(p.z, row[9]))
  }
  println(
    "local cartesian forward: max |dxyz| vs 50-digit = \{lc_max(errs)} m, vs PROJ = \{lc_max(perrs)} m",
  )
  assert_true(lc_max(errs) < 2.0e-8)
  assert_true(lc_max(perrs) < 2.0e-8)
}

///|
test "local cartesian reverse matches 50-digit values" {
  let dpos = []
  let dh = []
  for row in local_cartesian_forward {
    let lc = lc_of(row[0].to_int(), row[1], row[2], row[3])
    let g = lc.reverse(row[7], row[8], row[9])
    dpos.push(lc_pos_err(row[4], row[5], g.lat, g.lon))
    dh.push(lc_abs_err(g.h, row[6]))
  }
  println(
    "local cartesian reverse: max position error = \{lc_max(dpos)} m, max |dh| = \{lc_max(dh)} m",
  )
  assert_true(lc_max(dpos) < 2.0e-8)
  assert_true(lc_max(dh) < 2.0e-8)
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
