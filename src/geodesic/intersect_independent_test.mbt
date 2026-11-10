// Geodesic intersections against tools/gen_intersect_independent.py: cases
// with known answers (equator and meridians, poles, nodes, mirror images),
// lengths by mpmath quadrature, not GeographicLib.

///|
test "intersections with known answers match mpmath" {
  let mut dx = 0.0
  let mut n = 0
  for row in ix_indep_closest {
    let ix = ix_pick(row[0])
    let r = ix.closest(
      row[1],
      row[2],
      row[3],
      row[4],
      row[5],
      row[6],
      x0=row[7],
      y0=row[8],
    ) catch {
      _ => abort("closest failed")
    }
    n = n + 1
    assert_eq(r.c, 0)
    dx = indep_max(dx, (r.x - row[9]).abs())
    dx = indep_max(dx, (r.y - row[10]).abs())
  }
  for row in ix_indep_next {
    let ix = ix_pick(row[0])
    let r = ix.next(row[1], row[2], row[3], row[4]) catch {
      _ => abort("next failed")
    }
    n = n + 1
    assert_eq(r.c, 0)
    // The crossing half a period ahead and the one half a period behind
    // are equally close; either may be returned.
    dx = indep_max(dx, (r.x.abs() - row[5]).abs())
    dx = indep_max(dx, (r.y.abs() - row[6]).abs())
    if row[1] != 90.0 {
      assert_true(r.x * r.y > 0.0)
    }
  }
  for row in ix_indep_segment {
    let ix = ix_pick(row[0])
    let r = ix.segment(
      row[1],
      row[2],
      row[3],
      row[4],
      row[5],
      row[6],
      row[7],
      row[8],
    ) catch {
      _ => abort("segment failed")
    }
    n = n + 1
    assert_eq(r.c, 0)
    assert_eq(r.segmode.to_double(), row[11])
    dx = indep_max(dx, (r.x - row[9]).abs())
    dx = indep_max(dx, (r.y - row[10]).abs())
  }
  println("intersect vs known answers: n=\{n} max |dx|, |dy| m=\{dx}")
  assert_eq(n, 56)
  assert_true(dx < 1.0e-8) // README: < 10 nm
}
