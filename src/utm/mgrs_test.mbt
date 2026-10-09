///|
test "MGRS encoding matches GeoTrans" {
  let mut mismatches = 0
  for c in mgrs_cases {
    let (lat, lon, prec, expected) = c
    let s = @utm.to_mgrs(lat, lon, precision=prec)
    if s != expected {
      mismatches += 1
      println("MGRS mismatch: \{lat}, \{lon} -> \{s} (GeoTrans \{expected})")
    }
  }
  println(
    "MGRS: \{mgrs_cases.length() - mismatches}/\{mgrs_cases.length()} identical to GeoTrans",
  )
  assert_eq(mismatches, 0)
}

///|
test "MGRS decode round trip" {
  for c in mgrs_cases {
    let (lat, lon, prec, expected) = c
    let d = @utm.from_mgrs(expected, centered=false)
    assert_eq(d.precision, prec)
    // re-encoding the south-west corner gives the same string
    let again = @utm.from_latlon(d.lat, d.lon, zone=d.zone)
    let s = again.to_mgrs(d.lat, precision=prec) catch { _ => "?" }
    let _ = lat
    let _ = lon
    // (for 100 km / 10 km squares the centre may fall in a neighbouring band)
    if s != expected && prec >= 2 {
      // corner may round down across a metre boundary; the centre must not
      let m = @utm.from_mgrs(expected)
      let again2 = @utm.from_latlon(m.lat, m.lon, zone=m.zone)
      assert_eq(again2.to_mgrs(m.lat, precision=prec), expected)
    }
  }
}

///|
test "MGRS examples and errors" {
  let d = @utm.from_mgrs("50S MK 48502 17797")
  assert_eq(d.zone, 50)
  assert_true(d.north)
  assert_true(
    (d.lat - 39.9087).abs() < 1.0e-5 && (d.lon - 116.3975).abs() < 1.0e-5,
  )
  assert_eq(@utm.from_mgrs("50smk4850217797").easting, d.easting)
  assert_true((try? @utm.from_mgrs("50SMK485021779")) is Err(_)) // odd digits
  assert_true((try? @utm.from_mgrs("61SMK")) is Err(_)) // bad zone
  assert_true((try? @utm.from_mgrs("50IMK")) is Err(_)) // I is not a band
  assert_true((try? @utm.to_mgrs(0.0, 0.0, precision=12)) is Err(_))
  assert_eq(@utm.to_mgrs(0.0, 0.0, precision=0), "31NAA")
}
