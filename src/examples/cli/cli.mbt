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
    #|  coord wgs2gcj <lat> <lon>
    #|  coord gcj2wgs <lat> <lon>
    #|  coord gcj2bd <lat> <lon>
    #|  coord bd2gcj <lat> <lon>
    #|  coord wgs2bd <lat> <lon>
    #|  coord bd2wgs <lat> <lon>
    #|
    #|geod uses the WGS-84 ellipsoid. Distances are metres. Angles are degrees.
    #|coord converts one point. gcj2wgs, bd2gcj, and bd2wgs are exact inverses.
    #|
    #|  moon run src/examples/cli -- geod inverse 39.9087 116.3975 31.2304 121.4737
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
  let rest : Array[String] = []
  for i in 2..<args.length() {
    rest.push(args[i])
  }
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
