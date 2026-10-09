///|
test "UTM matches PROJ (EPSG:326xx/327xx)" {
  let mut worst = 0.0
  for c in utm_cases {
    let (lat, lon, zone, north, e, n) = c
    let u = @utm.from_latlon(lat, lon)
    assert_eq(u.zone, zone)
    assert_eq(u.north, north == 1.0)
    worst = @cmp.maximum(
      worst,
      @cmp.maximum((u.easting - e).abs(), (u.northing - n).abs()),
    )
    let g = @utm.to_latlon(u.zone, u.north, u.easting, u.northing)
    assert_true((g.lat - lat).abs() < 1.0e-11 && (g.lon - lon).abs() < 1.0e-11)
  }
  println("UTM vs PROJ worst \{worst} m")
  assert_true(worst < 1.0e-8)
}

///|
test "UPS matches PROJ (EPSG:32661/32761)" {
  let mut worst = 0.0
  for c in ups_cases {
    let (lat, lon, north, e, n) = c
    let u = @utm.from_latlon(lat, lon)
    assert_eq(u.zone, 0)
    assert_eq(u.north, north == 1.0)
    worst = @cmp.maximum(
      worst,
      @cmp.maximum((u.easting - e).abs(), (u.northing - n).abs()),
    )
    let g = @utm.to_latlon(0, u.north, u.easting, u.northing)
    let mut dlon = (g.lon - lon).abs()
    if dlon > 180.0 {
      dlon = 360.0 - dlon
    }
    assert_true((g.lat - lat).abs() < 1.0e-11)
    assert_true(dlon * @math.cos(lat * 3.141592653589793 / 180.0) < 1.0e-11)
  }
  println("UPS vs PROJ worst \{worst} m")
  assert_true(worst < 1.0e-8)
}

///|
test "standard zones: regular, Norway, Svalbard, polar" {
  assert_eq(@utm.standard_zone(39.9, 116.4), 50) // Beijing
  assert_eq(@utm.standard_zone(22.54, 114.06), 50) // Shenzhen
  assert_eq(@utm.standard_zone(31.23, 121.47), 51) // Shanghai
  assert_eq(@utm.standard_zone(60.39, 5.32), 32) // Bergen (Norway exception)
  assert_eq(@utm.standard_zone(78.22, 15.65), 33) // Longyearbyen (Svalbard)
  assert_eq(@utm.standard_zone(78.0, 8.0), 31)
  assert_eq(@utm.standard_zone(78.0, 35.0), 37)
  assert_eq(@utm.standard_zone(0.0, 180.0), 1)
  assert_eq(@utm.standard_zone(0.0, -180.0), 1)
  assert_eq(@utm.standard_zone(85.0, 0.0), 0)
  assert_eq(@utm.standard_zone(-80.5, 0.0), 0)
  assert_eq(@utm.central_meridian(50), 117.0)
}

///|
test "errors and formatting" {
  assert_true(raises(() => ignore(@utm.from_latlon(91.0, 0.0))))
  assert_true(raises(() => ignore(@utm.to_latlon(61, true, 500000.0, 0.0))))
  let u = @utm.from_latlon(0.0, 3.0)
  assert_eq(u.to_string(), "31N 500000 0")
}
