// Behaviours promised in README.mbt.md that had no test.

///|
fn utm_err(f : () -> Unit raise @utm.UtmError) -> @utm.UtmError? {
  try f() catch {
    e => Some(e)
  } noraise {
    _ => None
  }
}

///|
test "standard_zone boundaries" {
  assert_eq(@utm.standard_zone(84.0, 0.0), 0)
  assert_true(@utm.standard_zone(83.999, 0.0) != 0)
  assert_true(@utm.standard_zone(-80.0, 0.0) != 0)
  assert_eq(@utm.standard_zone(-80.000001, 0.0), 0)
  assert_eq(@utm.standard_zone(78.0, 25.0), 35)
  assert_eq(@utm.standard_zone(78.0, 32.9), 35)
  assert_eq(@utm.standard_zone(78.0, 33.0), 37)
  assert_eq(@utm.standard_zone(78.0, 8.9), 31)
  assert_eq(@utm.standard_zone(78.0, 9.0), 33)
  assert_eq(@utm.standard_zone(78.0, 20.9), 33)
  assert_eq(@utm.standard_zone(78.0, 21.0), 35)
  assert_eq(@utm.standard_zone(78.0, 40.0), 37)
}

///|
test "bad zones, latitudes, and grid coordinates raise" {
  assert_true(
    utm_err(() => ignore(@utm.from_latlon(40.0, 116.0, zone=61)))
    is Some(@utm.InvalidZone(_)),
  )
  assert_true(
    utm_err(() => ignore(@utm.from_latlon(@double.not_a_number, 116.0)))
    is Some(@utm.InvalidLatitude(_)),
  )
  assert_true(
    utm_err(() => ignore(@utm.to_latlon(50, true, 1.0e8, 4.0e6)))
    is Some(@utm.OutOfRange(_)),
  )
  assert_true(
    utm_err(() => ignore(@utm.to_latlon(50, true, @double.infinity, 4.0e6)))
    is Some(@utm.OutOfRange(_)),
  )
  // UPS: the pole is (2 000 000, 2 000 000) and maps back to ±90.
  let n = @utm.to_latlon(0, true, 2.0e6, 2.0e6)
  assert_eq(n.lat, 90.0)
  let s = @utm.to_latlon(0, false, 2.0e6, 2.0e6)
  assert_eq(s.lat, -90.0)
  // Transverse Mercator at the pole is finite.
  let p = @utm.TransverseMercator::utm().forward(117.0, 90.0, 0.0)
  assert_true(!p.x.is_nan() && !p.y.is_nan())
}

///|
test "MGRS precisions, centring, tabs, and limits" {
  // Precision 11 (1 micrometre) round trips.
  let m = @utm.to_mgrs(39.9087, 116.3975, precision=11)
  let back = @utm.from_mgrs(m, centered=false)
  assert_true((back.lat - 39.9087).abs() < 1.0e-10)
  assert_true(
    utm_err(() => ignore(@utm.to_mgrs(39.9087, 116.3975, precision=-1)))
    is Some(@utm.OutOfRange(_)),
  )
  // Default: centre of the square; centered=false: the south-west corner.
  let c = @utm.from_mgrs("50S MK 48502 17797")
  assert_eq(c.easting, 448502.5)
  assert_eq(c.northing, 4417797.5)
  let sw = @utm.from_mgrs("50S MK 48502 17797", centered=false)
  assert_eq(sw.easting, 448502.0)
  assert_eq(sw.northing, 4417797.0)
  // Tabs and lower case are ignored like spaces.
  let t = @utm.from_mgrs("50s\tmk\t48502\t17797")
  assert_eq(t.easting, c.easting)
  // More than 22 digits raises.
  assert_true(
    utm_err(() => {
      ignore(@utm.from_mgrs("50SMK" + "123456789012" + "123456789012"))
    })
    is Some(@utm.InvalidMgrs(_)),
  )
  // Easting outside 100-900 km raises.
  let u = @utm.from_latlon(40.0, 117.0, zone=50)
  let far = { ..u, easting: 950000.0, }
  assert_true(
    utm_err(() => ignore(far.to_mgrs(40.0, precision=5)))
    is Some(@utm.OutOfRange(_)),
  )
}
