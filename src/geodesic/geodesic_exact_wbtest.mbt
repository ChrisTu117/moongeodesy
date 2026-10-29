///|
test "elliptic integrals match GeographicLib" {
  let e0 = Elliptic::from_k2(0.0)
  let half = @geomath.pi / 2.0
  assert_true((e0.ec - half).abs() < 1.0e-15)
  assert_true((e0.dc - half / 2.0).abs() < 1.0e-15)
  assert_true((e0.hc - half / 2.0).abs() < 1.0e-15)
  let e = Elliptic::from_k2(0.5)
  assert_true((e.ec - 1.3506438810476751).abs() < 1.0e-14)
  assert_true((e.dc - 1.0068615925073927).abs() < 1.0e-14)
  let en = Elliptic::from_parts(-0.3, -0.2, 1.3, 1.2)
  assert_true((en.ec - 1.6826930354888963).abs() < 1.0e-14)
  assert_true((en.hc - 0.72535779960468094).abs() < 1.0e-14)
  let sn = 0.6
  let cn = (1.0 - sn * sn).sqrt()
  let dn = en.delta_amp(sn, cn)
  assert_true((en.delta_e(sn, cn, dn) - -0.031527922086467042).abs() < 1.0e-14)
  assert_true((en.delta_h(sn, cn, dn) - 0.52757104079648987).abs() < 1.0e-14)
  assert_true((en.einv(0.4) - 0.39698911164997308).abs() < 1.0e-14)
  let ctau = (1.0 - 0.09).sqrt()
  assert_true((en.delta_einv(0.3, ctau) - 0.02003625722829605).abs() < 1.0e-14)
}
