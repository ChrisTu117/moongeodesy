///|
/// Result of one CLI invocation. `Ok` is printed to stdout. `Err` is printed
/// to stderr and the process exits with status 1.
enum Run {
  Ok(String)
  Err(String)
}

///|
fn usage() -> String {
  let text =
    #|moongeodesy
    #|
    #|Usage:
    #|  geod inverse <lat1> <lon1> <lat2> <lon2>
    #|  geod direct <lat1> <lon1> <azi1> <s12>
    #|  rhumb inverse <lat1> <lon1> <lat2> <lon2>
    #|  rhumb direct <lat1> <lon1> <azi12> <s12>
    #|  intersect closest <latx> <lonx> <azix> <laty> <lony> <aziy>
    #|  enu forward <lat0> <lon0> <lat> <lon> [h]
    #|  enu reverse <lat0> <lon0> <x> <y> <z>
    #|  coord wgs2gcj <lat> <lon>
    #|  coord gcj2wgs <lat> <lon>
    #|  coord gcj2bd <lat> <lon>
    #|  coord bd2gcj <lat> <lon>
    #|  coord wgs2bd <lat> <lon>
    #|  coord bd2wgs <lat> <lon>
    #|
    #|geod / rhumb / intersect / enu use the WGS-84 ellipsoid.
    #|Distances and ENU offsets are metres. Angles are degrees. Heights are
    #|ellipsoidal metres (default 0 for enu forward).
    #|coord converts one point. gcj2wgs, bd2gcj, and bd2wgs are exact inverses.
    #|
    #|  moon run src/examples/cli -- geod inverse 39.9087 116.3975 31.2304 121.4737
    #|  moon run src/examples/cli -- rhumb inverse 51.4778 -0.0015 35.6762 139.6503
    #|  moon run src/examples/cli -- intersect closest 0 0 45 45 10 135
    #|  moon run src/examples/cli -- enu forward 39.9087 116.3975 39.9097 116.3985 50
    #|  moon run src/examples/cli -- coord wgs2gcj 39.9087 116.3975
  text
}

///|
fn die(message : String) -> Run {
  Err("error: \{message}\n\n\{usage()}")
}

///|
fn is_runtime_prefix(token : String) -> Bool {
  token == "node" || token.contains("/") || token.contains("\\")
}

///|
/// Drop the executable path (and `node` / the script path on JavaScript).
fn user_args(argv : Array[String]) -> Array[String] {
  let mut start = 0
  while start < argv.length() && is_runtime_prefix(argv[start]) {
    start = start + 1
  }
  let out : Array[String] = []
  for i in start..<argv.length() {
    out.push(argv[i])
  }
  out
}

///|
fn parse_num(token : String) -> Double? {
  Some(@string.parse_double(token)) catch {
    _ => None
  }
}

///|
fn bad_token(tokens : Array[String]) -> String {
  for t in tokens {
    if parse_num(t) is None {
      return t
    }
  }
  ""
}

///|
fn rest_from(args : Array[String], start : Int32) -> Array[String] {
  let out : Array[String] = []
  for i in start..<args.length() {
    out.push(args[i])
  }
  out
}

///|
fn execute(argv : Array[String]) -> Run {
  let args = user_args(argv)
  if args.length() == 0 ||
    args[0] == "help" ||
    args[0] == "-h" ||
    args[0] == "--help" {
    return Ok(usage())
  }
  if args[0] == "geod" {
    return execute_geod(args)
  }
  if args[0] == "rhumb" {
    return execute_rhumb(args)
  }
  if args[0] == "intersect" {
    return execute_intersect(args)
  }
  if args[0] == "enu" {
    return execute_enu(args)
  }
  if args[0] == "coord" {
    return execute_coord(args)
  }
  die("unknown command `\{args[0]}`")
}

///|
fn execute_geod(args : Array[String]) -> Run {
  if args.length() < 2 {
    return die("geod needs a subcommand: inverse or direct")
  }
  let rest = rest_from(args, 2)
  match args[1] {
    "inverse" => geod_inverse(rest)
    "direct" => geod_direct(rest)
    other => die("unknown geod command `\{other}`")
  }
}

///|
fn geod_inverse(tokens : Array[String]) -> Run {
  if tokens.length() != 4 {
    return die("geod inverse needs 4 numbers: <lat1> <lon1> <lat2> <lon2>")
  }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
    ) {
    (Some(lat1), Some(lon1), Some(lat2), Some(lon2)) => {
      let r = @geodesic.Geodesic::wgs84().inverse(lat1, lon1, lat2, lon2)
      Ok(
        "distance_m: \{r.s12}\nazimuth1_deg: \{r.azi1}\nazimuth2_deg: \{r.azi2}",
      )
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn geod_direct(tokens : Array[String]) -> Run {
  if tokens.length() != 4 {
    return die("geod direct needs 4 numbers: <lat1> <lon1> <azi1> <s12>")
  }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
    ) {
    (Some(lat1), Some(lon1), Some(azi1), Some(s12)) => {
      let d = @geodesic.Geodesic::wgs84().direct(lat1, lon1, azi1, s12)
      Ok("lat2: \{d.lat2}\nlon2: \{d.lon2}\nazimuth2_deg: \{d.azi2}")
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn execute_rhumb(args : Array[String]) -> Run {
  if args.length() < 2 {
    return die("rhumb needs a subcommand: inverse or direct")
  }
  let rest = rest_from(args, 2)
  match args[1] {
    "inverse" => rhumb_inverse(rest)
    "direct" => rhumb_direct(rest)
    other => die("unknown rhumb command `\{other}`")
  }
}

///|
fn rhumb_inverse(tokens : Array[String]) -> Run {
  if tokens.length() != 4 {
    return die("rhumb inverse needs 4 numbers: <lat1> <lon1> <lat2> <lon2>")
  }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
    ) {
    (Some(lat1), Some(lon1), Some(lat2), Some(lon2)) => {
      let r = @geodesic.Rhumb::wgs84().inverse(lat1, lon1, lat2, lon2)
      Ok("distance_m: \{r.s12}\nazimuth_deg: \{r.azi12}\narea_m2: \{r.area12}")
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn rhumb_direct(tokens : Array[String]) -> Run {
  if tokens.length() != 4 {
    return die("rhumb direct needs 4 numbers: <lat1> <lon1> <azi12> <s12>")
  }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
    ) {
    (Some(lat1), Some(lon1), Some(azi12), Some(s12)) => {
      let d = @geodesic.Rhumb::wgs84().direct(lat1, lon1, azi12, s12)
      Ok("lat2: \{d.lat2}\nlon2: \{d.lon2}\narea_m2: \{d.area12}")
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn execute_intersect(args : Array[String]) -> Run {
  if args.length() < 2 {
    return die("intersect needs a subcommand: closest")
  }
  let rest = rest_from(args, 2)
  match args[1] {
    "closest" => intersect_closest(rest)
    other => die("unknown intersect command `\{other}`")
  }
}

///|
fn intersect_closest(tokens : Array[String]) -> Run {
  if tokens.length() != 6 {
    return die(
      "intersect closest needs 6 numbers: <latx> <lonx> <azix> <laty> <lony> <aziy>",
    )
  }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
      parse_num(tokens[4]),
      parse_num(tokens[5]),
    ) {
    (
      Some(lat_x),
      Some(lon_x),
      Some(azi_x),
      Some(lat_y),
      Some(lon_y),
      Some(azi_y),
    ) => {
      let cross = @geodesic.Intersect::wgs84().closest(
        lat_x, lon_x, azi_x, lat_y, lon_y, azi_y,
      ) catch {
        e => return die("intersect closest failed: \{e.to_repr()}")
      }
      Ok("x_m: \{cross.x}\ny_m: \{cross.y}\nc: \{cross.c}")
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn execute_enu(args : Array[String]) -> Run {
  if args.length() < 2 {
    return die("enu needs a subcommand: forward or reverse")
  }
  let rest = rest_from(args, 2)
  match args[1] {
    "forward" => enu_forward(rest)
    "reverse" => enu_reverse(rest)
    other => die("unknown enu command `\{other}`")
  }
}

///|
fn enu_forward(tokens : Array[String]) -> Run {
  if tokens.length() != 4 && tokens.length() != 5 {
    return die(
      "enu forward needs 4 or 5 numbers: <lat0> <lon0> <lat> <lon> [h]",
    )
  }
  let h_tok = if tokens.length() == 5 { tokens[4] } else { "0" }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
      parse_num(h_tok),
    ) {
    (Some(lat0), Some(lon0), Some(lat), Some(lon), Some(h)) => {
      let lc = @geodesic.LocalCartesian::wgs84(lat0, lon0)
      let enu = lc.forward(lat, lon, h)
      Ok("east_m: \{enu.x}\nnorth_m: \{enu.y}\nup_m: \{enu.z}")
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn enu_reverse(tokens : Array[String]) -> Run {
  if tokens.length() != 5 {
    return die("enu reverse needs 5 numbers: <lat0> <lon0> <x> <y> <z>")
  }
  match
    (
      parse_num(tokens[0]),
      parse_num(tokens[1]),
      parse_num(tokens[2]),
      parse_num(tokens[3]),
      parse_num(tokens[4]),
    ) {
    (Some(lat0), Some(lon0), Some(x), Some(y), Some(z)) => {
      let lc = @geodesic.LocalCartesian::wgs84(lat0, lon0)
      let back = lc.reverse(x, y, z)
      Ok("lat: \{back.lat}\nlon: \{back.lon}\nh_m: \{back.h}")
    }
    _ => die("not a number: `\{bad_token(tokens)}`")
  }
}

///|
fn execute_coord(args : Array[String]) -> Run {
  if args.length() != 4 {
    return die("coord needs <op> <lat> <lon>")
  }
  let op = args[1]
  match (parse_num(args[2]), parse_num(args[3])) {
    (None, _) => die("not a number: `\{args[2]}`")
    (_, None) => die("not a number: `\{args[3]}`")
    (Some(lat), Some(lon)) =>
      match convert(op, lat, lon) {
        Some(text) => Ok(text)
        None => die("unknown coord command `\{op}`")
      }
  }
}

///|
fn coord_line(lat : Double, lon : Double) -> String {
  "lat: \{lat}\nlon: \{lon}"
}

///|
fn convert(op : String, lat : Double, lon : Double) -> String? {
  let p = @cn.LatLon::new(lat, lon)
  match op {
    "wgs2gcj" => {
      let q = @cn.wgs84_to_gcj02(p)
      Some(coord_line(q.lat, q.lon))
    }
    "gcj2wgs" => {
      let q = @cn.gcj02_to_wgs84(p)
      Some(coord_line(q.lat, q.lon))
    }
    "gcj2bd" => {
      let q = @cn.gcj02_to_bd09(p)
      Some(coord_line(q.lat, q.lon))
    }
    "bd2gcj" => {
      let q = @cn.bd09_to_gcj02(p)
      Some(coord_line(q.lat, q.lon))
    }
    "wgs2bd" => {
      let q = @cn.wgs84_to_bd09(p)
      Some(coord_line(q.lat, q.lon))
    }
    "bd2wgs" => {
      let q = @cn.bd09_to_wgs84(p)
      Some(coord_line(q.lat, q.lon))
    }
    _ => None
  }
}
