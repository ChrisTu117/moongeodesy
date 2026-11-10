///|
fn expect_ok(args : Array[String]) -> String {
  match execute(args) {
    Ok(text) => text
    Err(text) => abort("expected ok, got \{text}")
  }
}

///|
fn expect_err(args : Array[String]) -> String {
  match execute(args) {
    Err(text) => text
    Ok(text) => abort("expected error, got \{text}")
  }
}

///|
test "geod inverse matches wgs84" {
  let text = expect_ok([
    "geod", "inverse", "39.9087", "116.3975", "31.2304", "121.4737",
  ])
  let r = @geodesic.Geodesic::wgs84().inverse(
    39.9087, 116.3975, 31.2304, 121.4737,
  )
  assert_eq(
    text,
    "distance_m: \{r.s12}\nazimuth1_deg: \{r.azi1}\nazimuth2_deg: \{r.azi2}",
  )
}

///|
test "geod direct matches wgs84" {
  let text = expect_ok([
    "geod", "direct", "22.5431", "114.0579", "90", "1000000",
  ])
  let d = @geodesic.Geodesic::wgs84().direct(22.5431, 114.0579, 90.0, 1000000.0)
  assert_eq(text, "lat2: \{d.lat2}\nlon2: \{d.lon2}\nazimuth2_deg: \{d.azi2}")
}

///|
test "geod inverse accepts a negative longitude" {
  let text = expect_ok(["geod", "inverse", "0", "0", "0", "-1"])
  let r = @geodesic.Geodesic::wgs84().inverse(0.0, 0.0, 0.0, -1.0)
  assert_eq(
    text,
    "distance_m: \{r.s12}\nazimuth1_deg: \{r.azi1}\nazimuth2_deg: \{r.azi2}",
  )
}

///|
test "rhumb inverse and direct match wgs84" {
  let inv = expect_ok([
    "rhumb", "inverse", "51.4778", "-0.0015", "35.6762", "139.6503",
  ])
  let h = @geodesic.Rhumb::wgs84().inverse(51.4778, -0.0015, 35.6762, 139.6503)
  assert_eq(
    inv,
    "distance_m: \{h.s12}\nazimuth_deg: \{h.azi12}\narea_m2: \{h.area12}",
  )
  let dir = expect_ok(["rhumb", "direct", "51.4778", "-0.0015", "60", "1000000"])
  let d = @geodesic.Rhumb::wgs84().direct(51.4778, -0.0015, 60.0, 1000000.0)
  assert_eq(dir, "lat2: \{d.lat2}\nlon2: \{d.lon2}\narea_m2: \{d.area12}")
}

///|
test "intersect closest matches wgs84" {
  let text = expect_ok([
    "intersect", "closest", "0", "0", "45", "45", "10", "135",
  ])
  let cross = @geodesic.Intersect::wgs84().closest(
    0.0, 0.0, 45.0, 45.0, 10.0, 135.0,
  ) catch {
    e => abort("\{e.to_repr()}")
  }
  assert_eq(text, "x_m: \{cross.x}\ny_m: \{cross.y}\nc: \{cross.c}")
}

///|
test "enu forward and reverse match LocalCartesian" {
  let fwd = expect_ok([
    "enu", "forward", "39.9087", "116.3975", "39.9097", "116.3985", "50",
  ])
  let lc = @geodesic.LocalCartesian::wgs84(39.9087, 116.3975)
  let enu = lc.forward(39.9097, 116.3985, 50.0)
  assert_eq(fwd, "east_m: \{enu.x}\nnorth_m: \{enu.y}\nup_m: \{enu.z}")
  let default_h = expect_ok([
    "enu", "forward", "39.9087", "116.3975", "39.9097", "116.3985",
  ])
  let enu0 = lc.forward(39.9097, 116.3985, 0.0)
  assert_eq(default_h, "east_m: \{enu0.x}\nnorth_m: \{enu0.y}\nup_m: \{enu0.z}")
  let rev = expect_ok([
    "enu", "reverse", "39.9087", "116.3975", "100", "200", "50",
  ])
  let back = lc.reverse(100.0, 200.0, 50.0)
  assert_eq(rev, "lat: \{back.lat}\nlon: \{back.lon}\nh_m: \{back.h}")
}

///|
test "coord conversions match cn" {
  let p = @cn.LatLon::new(39.9087, 116.3975)
  let cases : Array[(String, @cn.LatLon)] = [
    ("wgs2gcj", @cn.wgs84_to_gcj02(p)),
    ("gcj2wgs", @cn.gcj02_to_wgs84(p)),
    ("gcj2bd", @cn.gcj02_to_bd09(p)),
    ("bd2gcj", @cn.bd09_to_gcj02(p)),
    ("wgs2bd", @cn.wgs84_to_bd09(p)),
    ("bd2wgs", @cn.bd09_to_wgs84(p)),
  ]
  for pair in cases {
    let (op, q) = pair
    let text = expect_ok(["coord", op, "39.9087", "116.3975"])
    assert_eq(text, "lat: \{q.lat}\nlon: \{q.lon}")
  }
}

///|
test "runtime prefix is ignored" {
  let direct = expect_ok(["geod", "inverse", "0", "0", "0", "1"])
  let wasm = expect_ok(["/tmp/cli.wasm", "geod", "inverse", "0", "0", "0", "1"])
  assert_eq(wasm, direct)
  let js = expect_ok([
    "/usr/bin/node", "/tmp/cli.js", "geod", "inverse", "0", "0", "0", "1",
  ])
  assert_eq(js, direct)
}

///|
test "help lists every command" {
  let text = expect_ok([])
  let lines = [
    "geod inverse <lat1> <lon1> <lat2> <lon2>", "geod direct <lat1> <lon1> <azi1> <s12>",
    "rhumb inverse <lat1> <lon1> <lat2> <lon2>", "rhumb direct <lat1> <lon1> <azi12> <s12>",
    "intersect closest <latx> <lonx> <azix> <laty> <lony> <aziy>", "enu forward <lat0> <lon0> <lat> <lon> [h]",
    "enu reverse <lat0> <lon0> <x> <y> <z>", "coord wgs2gcj <lat> <lon>", "coord gcj2wgs <lat> <lon>",
    "coord gcj2bd <lat> <lon>", "coord bd2gcj <lat> <lon>", "coord wgs2bd <lat> <lon>",
    "coord bd2wgs <lat> <lon>", "mag <lat> <lon> <date> [h]", "lcc forward <stdlat1> <stdlat2> <lat0> <lon0> <lat> <lon>",
    "lcc reverse <stdlat1> <stdlat2> <lat0> <lon0> <x> <y>", "albers forward <stdlat1> <stdlat2> <lat0> <lon0> <lat> <lon>",
    "albers reverse <stdlat1> <stdlat2> <lat0> <lon0> <x> <y>",
  ]
  for line in lines {
    assert_true(text.contains(line), msg=line)
  }
  assert_eq(expect_ok(["help"]), text)
  assert_eq(expect_ok(["-h"]), text)
  assert_eq(expect_ok(["--help"]), text)
  let prefixed = expect_ok(["/tmp/cli.wasm", "--help"])
  assert_eq(prefixed, text)
}

///|
test "unknown command explains usage" {
  let text = expect_err(["nope"])
  assert_true(text.contains("unknown command"))
  assert_true(text.contains("geod inverse"))
  assert_true(text.contains("rhumb inverse"))
  assert_true(text.contains("intersect closest"))
  assert_true(text.contains("enu forward"))
  assert_true(text.contains("coord wgs2gcj"))
  assert_true(text.contains("coord bd2wgs"))
}

///|
test "missing arguments and bad numbers are errors" {
  let missing = expect_err(["geod", "inverse", "1", "2"])
  assert_true(missing.contains("geod inverse"))
  let bad = expect_err(["coord", "wgs2gcj", "north", "116"])
  assert_true(bad.contains("north"))
  let unknown = expect_err(["coord", "wgs2utm", "1", "2"])
  assert_true(unknown.contains("wgs2utm"))
  let extra = expect_err(["geod", "direct", "0", "0", "90", "1", "2"])
  assert_true(extra.contains("geod direct"))
  let rhumb_bad = expect_err(["rhumb", "inverse", "1", "2"])
  assert_true(rhumb_bad.contains("rhumb inverse"))
  let enu_bad = expect_err(["enu", "reverse", "0", "0", "1"])
  assert_true(enu_bad.contains("enu reverse"))
}

///|
test "mag matches the WMM2025 model for a date and a decimal year" {
  let wmm = @magnetic.MagneticModel::wmm2025()
  let year = @magnetic.decimal_year(2026, 11, 10)
  let e = wmm.field(39.9087, 116.3975, 0.0, year)
  let text = expect_ok(["mag", "39.9087", "116.3975", "2026-11-10"])
  assert_true(text.has_prefix("declination_deg: \{e.decl}\n"))
  assert_true(text.contains("h_nt: \{e.h}\n"))
  assert_true(text.has_suffix("compass: normal"))
  let e2 = wmm.field(-80.0, 240.0, 100000.0, 2027.5)
  let text2 = expect_ok(["mag", "-80", "240", "2027.5", "100000"])
  assert_true(text2.has_prefix("declination_deg: \{e2.decl}\n"))
  let pole = expect_ok(["mag", "89", "-121", "2025.0", "28000"])
  assert_true(pole.contains("blackout"))
}

///|
test "mag rejects bad input" {
  assert_true(expect_err(["mag", "39.9", "116.4"]).contains("mag needs"))
  assert_true(
    expect_err(["mag", "39.9", "116.4", "2026-02-30"]).contains("not a date"),
  )
  assert_true(
    expect_err(["mag", "39.9", "116.4", "2031-01-01"]).contains(
      "DateOutOfRange",
    ),
  )
  assert_true(expect_err(["mag", "x", "116.4", "2026.0"]).contains("`x`"))
}

///|
test "lcc and albers forward and reverse match the library" {
  let a = 6378137.0
  let f = 1.0 / 298.257223563
  let aea = @conic.ConicGrid::albers(
    a,
    f,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  )
  let p = aea.forward(39.9087, 116.3975)
  assert_eq(
    expect_ok([
      "albers", "forward", "25", "47", "0", "105", "39.9087", "116.3975",
    ]),
    "easting_m: \{p.easting}\nnorthing_m: \{p.northing}\nconvergence_deg: \{p.gamma}\nscale: \{p.k}",
  )
  let lcc = @conic.ConicGrid::lambert_2sp(
    a,
    f,
    lat_origin=0.0,
    lon_origin=105.0,
    stdlat1=25.0,
    stdlat2=47.0,
  )
  let q = lcc.reverse(1000000.0, 4000000.0)
  assert_eq(
    expect_ok(["lcc", "reverse", "25", "47", "0", "105", "1000000", "4000000"]),
    "lat: \{q.lat}\nlon: \{q.lon}\nconvergence_deg: \{q.gamma}\nscale: \{q.k}",
  )
  assert_true(expect_err(["lcc", "sideways"]).contains("forward or reverse"))
  assert_true(
    expect_err(["albers", "forward", "25", "47"]).contains("6 numbers"),
  )
  assert_true(
    expect_err(["lcc", "forward", "95", "47", "0", "105", "40", "116"]).contains(
      "InvalidStandardLatitude",
    ),
  )
}
