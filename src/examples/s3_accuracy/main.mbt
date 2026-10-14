///|
/// Scenario 3: accuracy against GeographicLib's official test set.
/// Each case is solved both ways (inverse: distance; direct: endpoint) and
/// the errors are grouped by distance. Run: `moon run src/examples/s3_accuracy`
fn main {
  let g = @geodesic.Geodesic::wgs84()
  let labels = ["< 1 000 km", "1 000-10 000 km", "> 10 000 km", "all"]
  let n = [0, 0, 0, 0]
  let inv = [0.0, 0.0, 0.0, 0.0]
  let dir = [0.0, 0.0, 0.0, 0.0]
  for c in cases {
    let s = c[6]
    let bin = if s < 1.0e6 { 0 } else if s < 1.0e7 { 1 } else { 2 }
    let ei = (g.inverse(c[0], c[1], c[3], c[4]).s12 - s).abs()
    let d = g.direct(c[0], c[1], c[2], s)
    let ed = g.inverse(c[3], c[4], d.lat2, d.lon2).s12
    for b in [bin, 3] {
      n[b] = n[b] + 1
      if ei > inv[b] {
        inv[b] = ei
      }
      if ed > dir[b] {
        dir[b] = ed
      }
    }
  }
  println(
    "| distance band | cases | max inverse |Δs12| (nm) | max direct position error (nm) |",
  )
  println("|---|---|---|---|")
  for i in 0..<4 {
    println(
      "| \{labels[i]} | \{n[i]} | \{fmt(inv[i] * 1.0e9)} | \{fmt(dir[i] * 1.0e9)} |",
    )
  }
}

///|
fn fmt(x : Double) -> String {
  ((x * 100.0).round() / 100.0).to_string()
}
