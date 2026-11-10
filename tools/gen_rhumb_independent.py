"""Rhumb-line references from textbook formulas in mpmath (45 digits).

Independent of GeographicLib: tools/gen_rhumb.py uses RhumbSolve, which is
the code src/geodesic/rhumb.mbt is ported from. Here (see tools/mpgeod.py):

* inverse: azimuth = atan2(lambda12, psi2 - psi1) with the isometric
  latitude psi = asinh(tan phi) - e atanh(e sin phi); distance = meridian
  arc difference / cos(azimuth), the arc being the quadrature of the
  meridian radius; on an east-west line, the parallel length
  N(phi) cos(phi) |lambda12|;
* direct: the meridian arc advances by s12 cos(azimuth) (Newton on the
  arc integral gives lat2), and lambda12 = tan(azimuth) (psi2 - psi1);
  past a pole the arc is reflected and lon2 is undefined;
* area12, the area between the line and the equator, by Green's theorem:
  the integral of A(phi) dlambda along the line, where A(phi) is the area
  from the equator to the parallel per radian of longitude.

Ellipsoids: WGS-84 (code 0) and its prolate twin, f = -1/298.257223563
(code 1). Cases cover the poles, latitudes within 1e-4 deg of a pole,
the antimeridian, opposite meridians, east-west lines, and lines that spiral
into a pole (long_unroll).

Usage (from the repository root):
  pip install mpmath
  python3 tools/gen_rhumb_independent.py && moon fmt
"""
import random

import mpgeod as G

OUT = "src/geodesic/rhumb_independent_data_test.mbt"
ELL = {0: (6378137.0, 1 / 298.257223563), 1: (6378137.0, -1 / 298.257223563)}


def inverse_cases():
    lats = [-90.0, -89.9999, -89.0, -60.0, -45.5, -10.0, 0.0, 0.5, 30.0, 45.0,
            75.0, 89.5, 89.9999, 90.0]
    lons = [(0.0, 0.0), (10.0, 20.0), (179.0, -179.0), (-179.5, 179.5),
            (0.0, 180.0), (-170.0, 10.0), (45.0, -135.5), (100.0, 100.000001),
            (-30.0, 150.0), (170.0, -100.0), (5.0, -5.0)]
    cases = []
    k = 0
    for la1 in lats:
        for la2 in lats:
            if abs(la1) == 90 and abs(la2) == 90 and la1 != la2:
                continue
            lo1, lo2 = lons[k % len(lons)]
            k += 1
            cases.append((0, la1, lo1, la2, lo2))
    # East-west and almost east-west lines.
    for la in (0.0, 45.0, -60.0, 89.9, -89.99):
        cases.append((0, la, -10.0, la, 50.0))
        cases.append((0, la, 179.9, la, -179.9))
        cases.append((0, la, 0.0, la + 1e-9, 120.0))
    rnd = random.Random(20261111)
    for _ in range(60):
        cases.append((0, rnd.uniform(-90, 90), rnd.uniform(-180, 180),
                      rnd.uniform(-90, 90), rnd.uniform(-180, 180)))
    for _ in range(30):
        la1 = rnd.uniform(-89, 89)
        cases.append((1, la1, rnd.uniform(-180, 180), rnd.uniform(-89, 89),
                      rnd.uniform(-180, 180)))
    for la1, lo1, la2, lo2 in ((90.0, 0.0, 10.0, 25.0), (-45.0, 179.0, 45.0, -179.0),
                               (30.0, 0.0, 30.0, 90.0), (0.0, 0.0, 89.9999, 0.0)):
        cases.append((1, la1, lo1, la2, lo2))
    return [tuple(float(x) for x in c) for c in cases]


def direct_cases():
    cases = []
    lats = [-89.99, -60.0, -1.0, 0.0, 33.3, 60.0, 89.9]
    azis = [0.0, 30.0, 45.5, 89.999, 90.0, 90.001, 135.0, 180.0, -45.0, -90.0,
            -179.99]
    dists = [1.0, 1000.0, 100000.0, 5000000.0, 20000000.0, -3000000.0]
    k = 0
    for la in lats:
        for az in azis:
            s = dists[k % len(dists)]
            lo = [0.0, 179.5, -179.5, 100.0][k % 4]
            k += 1
            cases.append((0, 0, la, lo, az, s))
    # Antimeridian crossings.
    cases += [(0, 0, 10.0, 179.5, 90.0, 200000.0), (0, 0, -40.0, -179.9, -60.0, 300000.0),
              (0, 0, 0.0, 179.0, 90.0, 500000.0)]
    # Lines that spiral into a pole, with longitude unrolled.
    for la, az, s in ((80.0, 89.0, 1000000.0), (80.0, 89.0, 1100000.0),
                      (-85.0, 135.0, 500000.0), (60.0, 80.0, 1700000.0),
                      (0.0, 179.0, 40000.0), (89.0, 45.0, 157000.0)):
        cases.append((0, 1, la, 179.0, az, s))
    # Past a pole: lon2 is NaN, lat2 is reflected.
    cases += [(0, 0, 80.0, 10.0, 0.0, 3000000.0), (0, 0, -70.0, 0.0, 170.0, 5000000.0),
              (0, 0, 10.0, 0.0, 10.0, 30000000.0)]
    rnd = random.Random(20261112)
    for _ in range(60):
        cases.append((0, 0, rnd.uniform(-89, 89), rnd.uniform(-180, 180),
                      rnd.uniform(-180, 180), rnd.uniform(0, 1.5e7)))
    for _ in range(30):
        cases.append((1, 0, rnd.uniform(-89, 89), rnd.uniform(-180, 180),
                      rnd.uniform(-180, 180), rnd.uniform(0, 1.0e7)))
    return [tuple(float(x) for x in c) for c in cases]


def main():
    ells = {c: G.Ellipsoid(*ae) for c, ae in ELL.items()}
    for e in ells.values():
        e.self_check()
    inv_rows = []
    for code, la1, lo1, la2, lo2 in inverse_cases():
        azi, s12, area = ells[code].rhumb_inverse(la1, lo1, la2, lo2)
        inv_rows.append([code, la1, lo1, la2, lo2, azi, s12, area])
    dir_rows = []
    n_pole = 0
    for code, unroll, la1, lo1, az, s in direct_cases():
        lat2, lon2, area, lon12 = ells[code].rhumb_direct(la1, lo1, az, s, unroll=bool(unroll))
        n_pole += lon2 is None
        dir_rows.append([code, unroll, la1, lo1, az, s, lat2, lon2, area, lon12])
    out = [
        "// Generated by tools/gen_rhumb_independent.py. Do not edit.",
        "// Reference: textbook rhumb-line formulas evaluated with mpmath at 45 digits",
        "// (isometric latitude, meridian-arc quadrature, Green's theorem for the area).",
        "// Not from GeographicLib. Ellipsoid code: 0 WGS-84, 1 prolate f = -1/298.257223563.",
        "// Inverse columns: code lat1 lon1 lat2 lon2 azi12 s12 area12 (azi12 NaN: undefined).",
        "// Direct columns: code unroll lat1 lon1 azi12 s12 lat2 lon2 area12 lon12",
        "// (lon12 is the unrolled longitude change; lon2, area12, lon12 NaN: the line",
        "// passes a pole).",
        "",
        "///|",
        "let rhumb_indep_inverse : Array[FixedArray[Double]] = [",
    ]
    for r in inv_rows:
        out.append("  [" + ", ".join(G.lit(v) for v in r) + "],")
    out += ["]", "", "///|", "let rhumb_indep_direct : Array[FixedArray[Double]] = ["]
    for r in dir_rows:
        out.append("  [" + ", ".join(G.lit(v) for v in r) + "],")
    out += ["]", ""]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print("wrote %d inverse and %d direct rows (%d past a pole) to %s"
          % (len(inv_rows), len(dir_rows), n_pole, OUT))


if __name__ == "__main__":
    main()
