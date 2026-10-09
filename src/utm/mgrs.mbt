// Military Grid Reference System (MGRS), following the conventions of
// GeographicLib::MGRS (Charles Karney, MIT/X11 License) / NGA.SIG.0012.

///|
let tile : Double = 100000.0

///|
let latband : String = "CDEFGHJKLMNPQRSTUVWX"

///|
let upsband : String = "ABYZ"

///|
let utmcols : Array[String] = ["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"]

///|
let utmrow : String = "ABCDEFGHJKLMNPQRSTUV"

///|
let upscols : Array[String] = [
  "JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ",
]

///|
let upsrows : Array[String] = ["ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP"]

///|
/// Minimum UPS easting (in 100 km units) for bands A, B, Y, Z.
let ups_mineasting : Array[Int] = [8, 20, 13, 20]

///|
/// Minimum UPS northing (in 100 km units), south and north.
let ups_minnorthing : Array[Int] = [8, 13]

///|
let max_precision : Int = 11

///|
fn char_at(s : String, i : Int) -> Char {
  s.to_array()[i]
}

///|
fn index_of(s : String, c : Char) -> Int {
  for i, ch in s.to_array() {
    if ch == c {
      return i
    }
  }
  -1
}

///|
fn pow10(n : Int) -> Int64 {
  let mut r = 1L
  for _ in 0..<n {
    r = r * 10L
  }
  r
}

///|
fn pad_digits(v : Int64, width : Int) -> String {
  let s = v.to_string()
  let buf = StringBuilder::new()
  for _ in 0..<(width - s.length()) {
    buf.write_char('0')
  }
  buf.write_string(s)
  buf.to_string()
}

///|
/// Encode a UTM/UPS coordinate as MGRS. `precision` is the number of digits
/// per coordinate: 5 → 1 m, 4 → 10 m, ..., 0 → 100 km square only, up to
/// 11 (1 µm). The latitude `lat` is used to pick the UTM latitude band.
pub fn UtmCoord::to_mgrs(
  self : UtmCoord,
  lat : Double,
  precision? : Int = 5,
) -> String raise UtmError {
  if precision < 0 || precision > max_precision {
    raise OutOfRange("MGRS precision must be in 0..=11")
  }
  let mult = 1000000.0 // micrometres
  let xi = (self.easting * mult).floor().to_int64()
  let yi = (self.northing * mult).floor().to_int64()
  let tile_u = 100000L * 1000000L
  let xh = (xi / tile_u).to_int()
  let yh = (yi / tile_u).to_int()
  let buf = StringBuilder::new()
  if self.zone == 0 {
    let z = (if self.north { 2 } else { 0 }) +
      (if self.easting >= 2000000.0 { 1 } else { 0 })
    let cols = upscols[z]
    let rows = upsrows[if self.north { 1 } else { 0 }]
    let ci = xh - ups_mineasting[z]
    let ri = yh - ups_minnorthing[if self.north { 1 } else { 0 }]
    if ci < 0 || ci >= cols.length() || ri < 0 || ri >= rows.length() {
      raise OutOfRange("UPS coordinate outside the MGRS polar grid")
    }
    buf.write_char(char_at(upsband, z))
    buf.write_char(char_at(cols, ci))
    buf.write_char(char_at(rows, ri))
  } else {
    if xh < 1 || xh > 8 {
      raise OutOfRange("UTM easting outside 100 km..900 km")
    }
    let band = latitude_band(lat)
    buf.write_string(pad_digits(self.zone.to_int64(), 2))
    buf.write_char(char_at(latband, band + 10))
    buf.write_char(char_at(utmcols[(self.zone - 1) % 3], xh - 1))
    let shift = if self.zone % 2 == 0 { 5 } else { 0 }
    buf.write_char(char_at(utmrow, (yh + shift) % 20))
  }
  if precision > 0 {
    let div = pow10(max_precision - precision)
    let ex = xi % tile_u / div
    let ny = yi % tile_u / div
    buf.write_string(pad_digits(ex, precision))
    buf.write_string(pad_digits(ny, precision))
  }
  buf.to_string()
}

///|
/// Encode geographic coordinates (degrees) as MGRS using the standard zone.
pub fn to_mgrs(
  lat : Double,
  lon : Double,
  precision? : Int = 5,
) -> String raise UtmError {
  from_latlon(lat, lon).to_mgrs(lat, precision~)
}

///|
/// A decoded MGRS reference.
pub(all) struct MgrsPoint {
  zone : Int
  north : Bool
  easting : Double
  northing : Double
  /// Digits per coordinate in the input (0..=11).
  precision : Int
  lat : Double
  lon : Double
} derive(Debug)

///|
pub extend MgrsPoint with @debug.Debug::{to_repr}

///|
/// Decode an MGRS string (spaces allowed, case-insensitive). By default the
/// centre of the referenced grid square is returned; with `centered=false`
/// its south-west corner.
pub fn from_mgrs(
  input : String,
  centered? : Bool = true,
) -> MgrsPoint raise UtmError {
  let chars : Array[Char] = []
  for c in input.to_array() {
    if c == ' ' || c == '\t' {
      continue
    }
    let u = if c >= 'a' && c <= 'z' {
      (c.to_int() - 32).unsafe_to_char()
    } else {
      c
    }
    chars.push(u)
  }
  let len = chars.length()
  let mut p = 0
  let mut zone = 0
  while p < len && p < 2 && chars[p] >= '0' && chars[p] <= '9' {
    zone = zone * 10 + (chars[p].to_int() - '0'.to_int())
    p += 1
  }
  if p == 0 && len > 0 && chars[0] >= '0' && chars[0] <= '9' {
    raise InvalidMgrs(input)
  }
  if p > 0 && (zone < 1 || zone > 60) {
    raise InvalidZone(zone)
  }
  if len - p < 3 {
    raise InvalidMgrs("\{input}: expected band and 100 km square letters")
  }
  let bandc = chars[p]
  let colc = chars[p + 1]
  let rowc = chars[p + 2]
  p += 3
  let ndig = len - p
  if ndig % 2 != 0 || ndig / 2 > max_precision {
    raise InvalidMgrs("\{input}: need an even number of digits (at most 22)")
  }
  let prec = ndig / 2
  let mut ed = 0L
  let mut nd = 0L
  for i in 0..<prec {
    let c1 = chars[p + i]
    let c2 = chars[p + prec + i]
    if !(c1 >= '0' && c1 <= '9' && c2 >= '0' && c2 <= '9') {
      raise InvalidMgrs("\{input}: non-digit in coordinates")
    }
    ed = ed * 10L + (c1.to_int() - '0'.to_int()).to_int64()
    nd = nd * 10L + (c2.to_int() - '0'.to_int()).to_int64()
  }
  let unit = pow10(5 - @cmp.minimum(prec, 5)).to_double() /
    (if prec > 5 { pow10(prec - 5).to_double() } else { 1.0 })
  let off = if centered { unit / 2.0 } else { 0.0 }
  let mut east = ed.to_double() * unit + off
  let mut north_v = nd.to_double() * unit + off
  if zone == 0 {
    let z = index_of(upsband, bandc)
    if z < 0 {
      raise InvalidMgrs("\{input}: UPS band must be A, B, Y or Z")
    }
    let northp = z >= 2
    let ci = index_of(upscols[z], colc)
    let ri = index_of(upsrows[if northp { 1 } else { 0 }], rowc)
    if ci < 0 || ri < 0 {
      raise InvalidMgrs("\{input}: invalid UPS 100 km square")
    }
    east = east + (ci + ups_mineasting[z]).to_double() * tile
    north_v = north_v +
      (ri + ups_minnorthing[if northp { 1 } else { 0 }]).to_double() * tile
    let g = to_latlon(0, northp, east, north_v)
    return {
      zone: 0,
      north: northp,
      easting: east,
      northing: north_v,
      precision: prec,
      lat: g.lat,
      lon: g.lon,
    }
  }
  let band = index_of(latband, bandc)
  if band < 0 {
    raise InvalidMgrs("\{input}: invalid latitude band")
  }
  let ci = index_of(utmcols[(zone - 1) % 3], colc)
  let ri0 = index_of(utmrow, rowc)
  if ci < 0 || ri0 < 0 {
    raise InvalidMgrs("\{input}: invalid UTM 100 km square")
  }
  let northp = band >= 10
  let shift = if zone % 2 == 0 { 5 } else { 0 }
  let ri = (ri0 - shift + 20) % 20
  east = east + (ci + 1).to_double() * tile
  // The row letter repeats every 2000 km; pick the cycle that puts the
  // point inside the stated latitude band.
  let lat_lo = (band - 10).to_double() * 8.0
  let lat_hi = if band == 19 { 84.0 } else { lat_lo + 8.0 }
  let tol = 0.5 // degrees; squares straddling a band edge are allowed
  for k in 0..<10 {
    let n = north_v + (ri + 20 * k).to_double() * tile
    if n > 10000000.0 {
      break
    }
    let g = to_latlon(zone, northp, east, n)
    if g.lat >= lat_lo - tol && g.lat <= lat_hi + tol {
      return {
        zone,
        north: northp,
        easting: east,
        northing: n,
        precision: prec,
        lat: g.lat,
        lon: g.lon,
      }
    }
  }
  raise InvalidMgrs("\{input}: grid square does not intersect latitude band")
}
