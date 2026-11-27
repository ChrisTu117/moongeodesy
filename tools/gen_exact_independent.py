"""GeodesicExact references by numerical quadrature in mpmath (45 digits).

Independent of GeographicLib: tools/gen_geodesic_exact.py uses GeodSolve -E,
the code src/geodesic/geodesic_exact.mbt is ported from. Here the geodesic
is reduced to Bessel's auxiliary sphere (Helmert 1880) and the two classical
integrals are evaluated by mpmath.quad (see tools/mpgeod.py):

  s12 = b * integral sqrt(1 + k^2 sin^2 sigma) dsigma,
  lambda12 = omega12 - f (2 - f) sin(alpha0)
             * integral dsigma / (1 + (1 - f) sqrt(1 + k^2 sin^2 sigma)),

with k^2 = e'^2 cos^2(alpha0). The direct problem solves the first integral
for sigma2 by Newton's method. The inverse problem solves the second for the
auxiliary longitude omega12 (Vincenty's 1975 iteration, with exact integrals
instead of series); every inverse solution is checked by running the direct
problem back to the end point. S12, the area between the geodesic and the
equator, is Green's theorem, the integral of A(phi) dlambda along the line.
GeodesicExact uses elliptic integrals in Carlson's form; nothing here does.

Ellipsoids: a = 6378137 m, oblate f = 1/298.257223563, 0.05, 0.1, 0.2, 0.3,
0.5 and prolate f = -0.05, -0.1, -0.2, -0.3, -0.5. Distances are fractions of
each ellipsoid's quarter meridian Q, so the cases scale with the shape.

--check-ode also integrates the geodesic equations directly in latitude,
longitude and azimuth (classical RK4 in double precision, 200 000 steps, no
auxiliary sphere) for 15 cases from f = -0.5 to 0.5 and requires agreement
with the quadrature to 1e-10 deg. It checks the formulas, not the digits,
and takes a few minutes.

Usage (from the repository root):
  pip install mpmath
  python3 tools/gen_exact_independent.py && moon fmt
  python3 tools/gen_exact_independent.py --check-ode
"""
import math
import random
import sys

import mpmath as mp

import mpgeod as G

OUT = "src/geodesic/geodesic_exact_independent_data_test.mbt"
A = 6378137.0
FLATTENINGS = [1 / 298.257223563, 0.05, 0.1, 0.2, 0.3, 0.5,
               -0.05, -0.1, -0.2, -0.3, -0.5]

# Direct: lat1, lon1, azi1, s12 as a multiple of Q (or metres when absolute).
DIRECT = [
    (40.6, -73.8, 45.0, 1.0, False),
    (-33.0, 151.0, 270.0, 0.5, False),
    (10.0, 20.0, 123.0, 0.8, False),
    (0.0, 0.0, 90.0, 0.3, False),       # along the equator
    (0.0, 0.0, 89.999, 1.9, False),     # just off the equator, long
    (80.0, 0.0, 0.0001, 0.3, False),    # passes 1 m-ish from the pole
    (80.0, 0.0, 0.0, 0.3, False),       # a meridian over the pole
    (-60.0, 179.5, 60.0, 0.2, False),   # across the antimeridian
    (20.0, 0.0, 30.0, 3.5, False),      # most of the way round twice
    (45.0, 10.0, 135.0, -0.6, False),   # negative distance
    (-0.5, 179.0, 179.5, 1.2, False),
    (1.0, 2.0, 45.0, 0.15, True),       # 15 cm
    (89.9, 0.0, 89.95, 20.0, True),     # 20 m next to the pole
]

# Inverse: |lon12| <= 60 deg and away from conjugate points, so the solution
# found from omega12 = lambda12 is the shortest geodesic.
INVERSE = [
    (40.6, -73.8, 1.5, -20.0),
    (10.0, 20.0, -30.0, 70.0),
    (60.0, 30.0, -10.0, 80.0),
    (0.0, 0.0, 0.5, 50.0),
    (-45.0, 170.0, 30.0, -160.0),       # across the antimeridian
    (5.0, 10.0, 5.0, 60.0),             # same latitude
    (-20.0, 0.0, 20.0, 0.0),            # a meridian
    (1.0, 2.0, 1.000001, 2.000001),     # 15 cm
    (89.9, 0.0, 89.9, 0.1),
    (-70.0, -100.0, -75.0, -55.0),
]


def rk4_direct(a, f, lat1, lon1, azi1, s12, n=200000):
    """Integrate dphi/ds = cos(alpha)/M, dlambda/ds = sin(alpha)/(N cos phi),
    dalpha/ds = sin(alpha) tan(phi)/N with classical RK4 (doubles)."""
    e2 = f * (2 - f)

    def F(y):
        sp, cp = math.sin(y[0]), math.cos(y[0])
        w = 1 - e2 * sp * sp
        M = a * (1 - e2) / w ** 1.5
        N = a / math.sqrt(w)
        return (math.cos(y[2]) / M, math.sin(y[2]) / (N * cp), math.sin(y[2]) * sp / (N * cp))

    y = (math.radians(lat1), math.radians(lon1), math.radians(azi1))
    h = s12 / n
    for _ in range(n):
        k1 = F(y)
        k2 = F(tuple(y[i] + h / 2 * k1[i] for i in range(3)))
        k3 = F(tuple(y[i] + h / 2 * k2[i] for i in range(3)))
        k4 = F(tuple(y[i] + h * k3[i] for i in range(3)))
        y = tuple(y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(3))
    return tuple(math.degrees(v) for v in y)


def check_ode():
    worst = 0.0
    for f in (0.5, 0.2, 1 / 298.257223563, -0.2, -0.5):
        E = G.Ellipsoid(A, f)
        Q = float(E.quarter_meridian())
        for lat1, lon1, azi1, k in ((10.0, 20.0, 123.0, 0.8), (-33.0, 151.0, 270.0, 0.5),
                                    (40.6, -73.8, 45.0, 1.0)):
            d = E.geod_direct(lat1, lon1, azi1, k * Q, area=False)
            o = rk4_direct(A, f, lat1, lon1, azi1, k * Q)
            dz = abs(float(d["azi2"]) - o[2]) % 360
            errs = (abs(float(d["lat2"]) - o[0]), abs(float(d["lon2_unrolled"]) - o[1]),
                    min(dz, 360 - dz))
            worst = max(worst, *errs)
            print("f=%.4f %s: dlat %.1e, dlon %.1e, dazi %.1e deg" % ((f, (lat1, lon1, azi1, k)) + errs))
    print("largest difference %.1e deg" % worst)
    assert worst < 1e-10


def main():
    if "--check-ode" in sys.argv:
        check_ode()
        return
    rnd = random.Random(20261112)
    rows = []
    n_dir = n_inv = 0
    for code, f in enumerate(FLATTENINGS):
        E = G.Ellipsoid(A, f)
        E.self_check()
        Q = E.quarter_meridian()
        cases = list(DIRECT)
        for _ in range(6):
            cases.append((rnd.uniform(-85, 85), rnd.uniform(-180, 180),
                          rnd.uniform(-180, 180), rnd.uniform(-2.0, 2.0), False))
        for lat1, lon1, azi1, s, absolute in cases:
            s12 = float(s if absolute else s * Q)
            d = E.geod_direct(lat1, lon1, azi1, s12)
            rows.append([code, 1, lat1, lon1, azi1, s12, d["lat2"], d["lon2"], d["azi2"],
                         d["S12"]])
            n_dir += 1
        icases = list(INVERSE)
        for _ in range(5):
            la1 = rnd.uniform(-80, 80)
            lo1 = rnd.uniform(-180, 180)
            icases.append((la1, lo1, rnd.uniform(-80, 80), lo1 + rnd.uniform(-60, 60)))
        for lat1, lon1, lat2, lon2 in icases:
            v = E.geod_inverse(lat1, lon1, lat2, lon2)
            # Check: the direct problem from the solution returns to point 2.
            back = E.geod_direct(lat1, lon1, v["azi1"], v["s12"], area=False)
            assert abs(back["lat2"] - G.D(lat2)) < 1e-20, (f, lat1, lon1, lat2, lon2)
            assert abs(G.wrap180(back["lon2"] - G.D(lon2))) < 1e-20, (f, lat1, lon1, lat2, lon2)
            rows.append([code, 0, lat1, lon1, lat2, lon2, v["s12"], v["azi1"], v["azi2"],
                         v["S12"]])
            n_inv += 1
    out = [
        "// Generated by tools/gen_exact_independent.py. Do not edit.",
        "// Reference: numerical quadrature (mpmath, 45 digits) of the classical",
        "// distance and longitude integrals on the auxiliary sphere; not GeographicLib.",
        "// a = 6378137 m; the flattening of code k is exact_indep_f[k].",
        "// Direct rows: code 1 lat1 lon1 azi1 s12 lat2 lon2 azi2 S12",
        "// Inverse rows: code 0 lat1 lon1 lat2 lon2 s12 azi1 azi2 S12",
        "// S12 NaN: undefined (a meridian over a pole).",
        "",
        "///|",
        "let exact_indep_f : Array[Double] = [" + ", ".join(G.lit(f) for f in FLATTENINGS) + "]",
        "",
        "///|",
        "let exact_indep : Array[FixedArray[Double]] = [",
    ]
    for r in rows:
        out.append("  [" + ", ".join(G.lit(v) for v in r) + "],")
    out += ["]", ""]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print("wrote %d direct and %d inverse rows to %s" % (n_dir, n_inv, OUT))


if __name__ == "__main__":
    main()
