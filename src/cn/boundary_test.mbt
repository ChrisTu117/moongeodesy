// Tests for the China / foreign-land boundary used by GCJ-02.
// City coordinates are city centres, each well away (>= 10 km) from any land
// border, so the ~1 km simplification of the borders cannot flip them.

///|
/// Wraps a round-trip check so failures name the point.
fn check_round_trip(lat : Double, lon : Double) -> Unit raise {
  let w = @cn.LatLon::new(lat, lon)
  let g = @cn.wgs84_to_gcj02(w)
  let back = @cn.gcj02_to_wgs84(g)
  // The inverse always returns a WGS-84 source of g.
  let again = @cn.wgs84_to_gcj02(back)
  assert_true(
    close(again.lat, g.lat, 1.0e-9) && close(again.lon, g.lon, 1.0e-9),
    msg="source of g at \{lat},\{lon}",
  )
  // Inside China the source is the original point.
  if !@cn.out_of_china(lat, lon) {
    assert_true(
      close(back.lat, lat, 1.0e-9) && close(back.lon, lon, 1.0e-9),
      msg="round trip at \{lat},\{lon}",
    )
  }
}

///|
/// Acceptance feedback: coordinates abroad (Seoul, Hanoi, ...) were offset
/// because only the bounding box was checked.
test "foreign cities are not offset" {
  let cities : Array[(String, Double, Double)] = [
    ("Seoul", 37.5665, 126.978),
    ("Busan", 35.1796, 129.0756),
    ("Pyongyang", 39.0392, 125.7625),
    ("Osaka", 34.6937, 135.5023),
    ("Hanoi", 21.0285, 105.8542),
    ("Ho Chi Minh City", 10.8231, 106.6297),
    ("Vientiane", 17.9757, 102.6331),
    ("Bangkok", 13.7563, 100.5018),
    ("Phnom Penh", 11.5564, 104.9282),
    ("Yangon", 16.8409, 96.1735),
    ("Kuala Lumpur", 3.139, 101.6869),
    ("Singapore", 1.3521, 103.8198),
    ("Manila", 14.5995, 120.9842),
    ("Ulaanbaatar", 47.8864, 106.9057),
    ("Vladivostok", 43.1198, 131.8869),
    ("Almaty", 43.222, 76.8512),
    ("Bishkek", 42.8746, 74.5698),
    ("Kathmandu", 27.7172, 85.324),
    ("Thimphu", 27.4728, 89.639),
    ("Dhaka", 23.8103, 90.4125),
    ("New Delhi", 28.6139, 77.209),
    ("Islamabad", 33.6844, 73.0479),
  ]
  for city in cities {
    let (name, lat, lon) = city
    let p = @cn.LatLon::new(lat, lon)
    assert_false(@cn.out_of_china_bbox(lat, lon), msg="\{name} is in the box")
    assert_true(@cn.out_of_china(lat, lon), msg=name)
    assert_true(@cn.wgs84_to_gcj02(p) == p, msg="\{name} forward")
    assert_true(@cn.gcj02_to_wgs84(p) == p, msg="\{name} inverse")
    assert_true(@cn.gcj02_to_wgs84_approx(p) == p, msg="\{name} approx")
    let back = @cn.bd09_to_wgs84(@cn.wgs84_to_bd09(p))
    assert_true(
      close(back.lat, lat, 1.0e-9) && close(back.lon, lon, 1.0e-9),
      msg="\{name} BD-09 round trip",
    )
  }
}

///|
test "Chinese cities are offset and round-trip" {
  let cities : Array[(String, Double, Double)] = [
    ("Beijing", 39.9087, 116.3975),
    ("Shanghai", 31.2304, 121.4737),
    ("Shenzhen", 22.5431, 114.0579),
    ("Kashgar", 39.4704, 75.9898),
    ("Urumqi", 43.8256, 87.6168),
    ("Lhasa", 29.65, 91.1),
    ("Kunming", 25.0389, 102.7183),
    ("Nanning", 22.817, 108.3665),
    ("Harbin", 45.8038, 126.5349),
    ("Hohhot", 40.8424, 111.749),
    ("Haikou", 20.044, 110.1983),
    ("Sanya", 18.2528, 109.5119),
    ("Hong Kong", 22.3193, 114.1694),
    ("Macau", 22.1987, 113.5439),
    ("Taipei", 25.033, 121.5654),
    ("Kaohsiung", 22.6273, 120.3014),
    ("Yongxing Island", 16.834, 112.338),
  ]
  for city in cities {
    let (name, lat, lon) = city
    let w = @cn.LatLon::new(lat, lon)
    assert_false(@cn.out_of_china(lat, lon), msg=name)
    let g = @cn.wgs84_to_gcj02(w)
    assert_true(
      (g.lat - lat).abs() > 1.0e-4 || (g.lon - lon).abs() > 1.0e-4,
      msg="\{name} offset",
    )
    let back = @cn.gcj02_to_wgs84(g)
    assert_true(
      close(back.lat, lat, 1.0e-9) && close(back.lon, lon, 1.0e-9),
      msg="\{name} exact inverse",
    )
    let approx = @cn.gcj02_to_wgs84_approx(g)
    assert_true(
      close(approx.lat, lat, 5.0e-5) && close(approx.lon, lon, 5.0e-5),
      msg="\{name} approximate inverse",
    )
  }
}

///|
/// Acceptance feedback: next to the box edges the exact inverse was off by
/// about 0.0055° when the GCJ-02 value crossed an edge.
test "round trip next to the box edges" {
  for i in 0..<400 {
    let t = i.to_double() / 400.0
    check_round_trip(0.8293 + 0.02 * t, 105.0 + 4.0 * t) // south
    check_round_trip(55.8271 - 0.02 * t, 120.0 + 17.0 * t) // north
    check_round_trip(15.0 + 10.0 * t, 72.004 + 0.02 * t) // west
    check_round_trip(25.0 + 10.0 * t, 137.8347 - 0.02 * t) // east
  }
}

///|
/// A coarse grid over the whole box, inside and outside China.
test "round trip over the whole box" {
  for i in 0..<=55 {
    for j in 0..<=66 {
      check_round_trip(0.83 + i.to_double(), 72.01 + j.to_double())
    }
  }
}
