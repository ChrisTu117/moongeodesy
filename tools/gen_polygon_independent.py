"""Polygon area references that do not come from GeographicLib.

tools/gen_polygon.py uses Python geographiclib, which src/geodesic/polygon.mbt
is ported from. This script uses three independent sources:

1. Closed forms on the ellipsoid (WGS-84 and its prolate twin), where every
   edge is a meridian or an arc of the equator (both are geodesics; an
   equatorial arc is the shortest path when it spans less than (1 - f) 180
   deg): the octant is 1/8 of the ellipsoid's area, a lune from the equator
   to a pole spanning dlon is dlon / 360 of a hemisphere, a lune from pole to
   pole is dlon / 360 of the ellipsoid, and the equator bounds a hemisphere.
   The ellipsoid's area is 2 pi a^2 (1 + (1 - e^2) atanh(e) / e), checked in
   mpgeod.Ellipsoid.self_check against the surface-of-revolution integral.
   Perimeters are equator arcs a * dlon plus meridian quadrants by quadrature.
2. Spherical triangles (f = 0): area R^2 E with the spherical excess E from
   L'Huilier's theorem, sides from the haversine formula.
3. The eight polygons of tools/gen_polygon.py on WGS-84, by Green's theorem:
   area = -sum over edges of the integral of A(phi) dlambda, plus half the
   ellipsoid when the boundary winds round the pole an odd number of times,
   reduced to (-area/2, area/2]. Each edge is solved by quadrature on the
   auxiliary sphere (tools/mpgeod.py, as in gen_exact_independent.py); a
   vertex at a pole adds A(pole) times the turn in longitude there.

All values at 45 digits with mpmath; counter-clockwise is positive, as in
Geodesic::polygon.

Usage (from the repository root):
  pip install mpmath
  python3 tools/gen_polygon_independent.py && moon fmt
"""
import mpmath as mp

import mpgeod as G

OUT = "src/geodesic/polygon_independent_data_test.mbt"
WGS84 = (6378137.0, 1 / 298.257223563)
PROLATE = (6378137.0, -1 / 298.257223563)
SPHERE = (6371000.0, 0.0)

# The eight polygons of tools/gen_polygon.py.
POLYS = {
    "square_1deg": [(0, 0), (0, 1), (1, 1), (1, 0)],
    "clockwise_square": [(0, 0), (1, 0), (1, 1), (0, 1)],
    "antimeridian_triangle": [(10, 179), (10, -179), (12, 180)],
    "north_pole_cap": [(89, 0), (89, 90), (89, 180), (89, 270)],
    "antarctica_like": [(-72.9, -74), (-71.9, -102), (-74.9, -102), (-74.3, -131),
                        (-77.5, -163), (-77.4, 163), (-71.7, 172), (-65.9, 140),
                        (-65.7, 113), (-66.6, 88), (-66.9, 59), (-69.8, 25), (-70.0, -4),
                        (-71.0, -14), (-77.3, -33), (-77.9, -46), (-74.7, -61)],
    "china_rough": [(53.5, 123.3), (48.4, 135.0), (39.8, 124.3), (21.5, 108.0),
                    (21.4, 101.8), (28.0, 97.4), (27.9, 88.8), (35.6, 78.0), (40.3, 73.6),
                    (45.2, 82.5), (49.1, 87.3), (42.8, 96.4), (42.6, 111.0), (49.9, 116.7)],
    "planimeter_reference": [(0, 0), (0, 90), (90, 0)],
    "big_quadrant": [(0, 0), (0, 90), (0, 180), (0, -90)],
}


def reduce_signed(area, total):
    """Reduce to (-total/2, total/2], as Geodesic::polygon does with sign."""
    a = mp.fmod(area, total)
    if a > total / 2:
        a -= total
    elif a <= -total / 2:
        a += total
    return a


def closed_forms(E, tag):
    """Polygons of meridians and equator arcs with exact areas."""
    A = E.total_area()
    Q = E.quarter_meridian()
    a = E.a
    out = []
    # Octant: equator 0..90, meridian 90 to the pole, meridian 0 back.
    out.append(("octant" + tag, [(0, 0), (0, 90), (90, 0)], a * mp.pi / 2 + 2 * Q, A / 8))
    # Lunes from the equator to the north pole (counter-clockwise).
    for lon1, lon2, plon in ((0, 1, 0), (20, 50, 123), (-60, 60, 0), (10, 189, -45),
                             (170, 190, 180)):
        dl = mp.mpf(lon2 - lon1)
        pts = [(0, G.wrap180(mp.mpf(lon1))), (0, G.wrap180(mp.mpf(lon2))), (90, plon)]
        out.append(("north_lune_%g%s" % (lon2 - lon1, tag), pts, a * G.rad(dl) + 2 * Q,
                    dl / 360 * A / 2))
    # Lunes from the equator to the south pole (counter-clockwise: west along
    # the equator, then down to the pole and back up).
    for lon1, lon2 in ((30, 75), (-179, -1)):
        dl = mp.mpf(lon2 - lon1)
        pts = [(0, lon2), (0, lon1), (-90, 0)]
        out.append(("south_lune_%g%s" % (lon2 - lon1, tag), pts, a * G.rad(dl) + 2 * Q,
                    dl / 360 * A / 2))
    # Pole-to-pole lunes: down meridian lon1, up meridian lon2 (east of it).
    for lon1, lon2 in ((0, 90), (-100, -40), (150, 210)):
        dl = mp.mpf(lon2 - lon1)
        pts = [(90, 0), (0, G.wrap180(mp.mpf(lon1))), (-90, 0), (0, G.wrap180(mp.mpf(lon2)))]
        out.append(("pole_lune_%g%s" % (lon2 - lon1, tag), pts, 4 * Q, dl / 360 * A))
    # The equator, eastward: the northern hemisphere.
    out.append(("equator" + tag, [(0, 0), (0, 90), (0, 180), (0, -90)], 2 * mp.pi * a, A / 2))
    # Westward round the equator from 45: the southern hemisphere (same area).
    out.append(("equator_west" + tag, [(0, 45), (0, -45), (0, -135), (0, 135)], 2 * mp.pi * a,
                reduce_signed(A / 2, A)))
    return out


def unit(lat, lon):
    p, l = G.rad(G.D(lat)), G.rad(G.D(lon))
    return mp.matrix([mp.cos(p) * mp.cos(l), mp.cos(p) * mp.sin(l), mp.sin(p)])


def sphere_triangle(R, pts):
    """Perimeter and signed area of a spherical triangle (minor-arc sides)."""
    v = [unit(*p) for p in pts]

    def side(u, w):
        return 2 * mp.asin(mp.norm(u - w) / 2)

    a, b, c = side(v[1], v[2]), side(v[2], v[0]), side(v[0], v[1])
    s = (a + b + c) / 2
    t = mp.tan(s / 2) * mp.tan((s - a) / 2) * mp.tan((s - b) / 2) * mp.tan((s - c) / 2)
    E = 4 * mp.atan(mp.sqrt(t))
    triple = (v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1])
              - v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0])
              + v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0]))
    sign = 1 if triple > 0 else -1
    return R * (a + b + c), sign * R * R * E


SPHERE_TRIANGLES = {
    "sphere_small": [(0, 0), (0, 1), (1, 0)],
    "sphere_clockwise": [(10, 10), (20, 10), (15, 20)],
    "sphere_big": [(-30, -40), (10, 60), (50, -20)],
    "sphere_antimeridian": [(-5, 170), (8, -165), (25, 175)],
    "sphere_pole_cap": [(10, 0), (10, 120), (10, -120)],
    "sphere_south_pole": [(-60, 0), (-70, -100), (-65, 130)],
    "sphere_thin": [(0, 0), (0.001, 60), (0, 120)],
}


def green_polygon(E, pts):
    """Perimeter and signed area of a geodesic polygon by Green's theorem."""
    A = E.total_area()
    n = len(pts)
    per = mp.mpf(0)
    integral = mp.mpf(0)   # sum of the integral of A(phi) dlambda
    turn = mp.mpf(0)       # net change of longitude (radians)
    for i in range(n):
        la1, lo1 = pts[i]
        la2, lo2 = pts[(i + 1) % n]
        p1, p2 = abs(la1) == 90, abs(la2) == 90
        if p1 or p2:
            assert not (p1 and p2)
            # A meridian from or to the pole: no longitude change along it.
            other = la2 if p1 else la1
            per += abs(E.quarter_meridian() - E.meridian_arc(G.rad(G.D(abs(other)))))
            continue
        if la1 == 0 and la2 == 0:
            # Along the equator, the shortest path while |dlon| < (1 - f) 180.
            d = G.rad(G.ang_diff(G.D(lo1), G.D(lo2)))
            assert abs(d) < (1 - E.f) * mp.pi
            per += E.a * abs(d)
            turn += d
            continue
        r = E.geod_inverse(la1, lo1, la2, lo2)
        per += r["s12"]
        integral += r["S12"]
        turn += r["lam12"]
    # Turns at polar vertices: from the meridian of the previous vertex to the
    # meridian of the next one, the short way.
    for i in range(n):
        la, _ = pts[i]
        if abs(la) == 90:
            lin = pts[i - 1][1]
            lout = pts[(i + 1) % n][1]
            d = G.rad(G.ang_diff(G.D(lin), G.D(lout)))
            integral += E.area_fn(G.rad(G.D(la))) * d
            turn += d
    k = mp.nint(turn / (2 * mp.pi))
    assert abs(turn - 2 * mp.pi * k) < mp.mpf(10) ** -30
    area = -integral + (A / 2 if int(k) % 2 else 0)
    return per, reduce_signed(area, A)


def main():
    rows = []
    for code, (ae, tag) in enumerate(((WGS84, ""), (PROLATE, "_prolate"))):
        E = G.Ellipsoid(*ae)
        E.self_check()
        for name, pts, per, area in closed_forms(E, tag):
            rows.append((name, code, pts, per, area, "closed form"))
    S = G.Ellipsoid(*SPHERE)
    for name, pts in SPHERE_TRIANGLES.items():
        per, area = sphere_triangle(S.a, pts)
        rows.append((name, 2, pts, per, area, "L'Huilier"))
    W = G.Ellipsoid(*WGS84)
    for name, pts in POLYS.items():
        per, area = green_polygon(W, pts)
        rows.append((name, 0, pts, per, area, "Green"))
    # Check the Green's-theorem code on the closed forms and the triangles too.
    for name, code, pts, per, area, how in rows:
        if how == "Green":
            continue
        E = G.Ellipsoid(*(WGS84, PROLATE, SPHERE)[code])
        gp, ga = green_polygon(E, pts)
        assert abs(gp - per) < mp.mpf(10) ** -25 * per, (name, gp, per)
        assert abs(ga - area) < mp.mpf(10) ** -25 * E.total_area(), (name, ga, area)

    def pt(p):
        return "(%s, %s)" % (G.lit(p[0]), G.lit(p[1]))

    out = [
        "// Generated by tools/gen_polygon_independent.py. Do not edit.",
        "// Reference: closed forms (octant, lunes, hemispheres), L'Huilier's theorem on",
        "// a sphere, and Green's theorem with edges by quadrature (mpmath, 45 digits).",
        "// Not from geographiclib. Ellipsoid code: 0 WGS-84, 1 prolate",
        "// f = -1/298.257223563 (a = 6378137 m), 2 sphere R = 6371000 m.",
        "// Columns: name, code, vertices (lat, lon), perimeter (m), area (m^2,",
        "// counter-clockwise positive, in (-A/2, A/2]), source.",
        "",
        "///|",
        "let polygon_indep_cases : Array[(String, Int, Array[(Double, Double)], Double, Double, String)] = [",
    ]
    for name, code, pts, per, area, how in rows:
        out.append('  ("%s", %d, [%s], %s, %s, "%s"),'
                   % (name, code, ", ".join(pt(p) for p in pts), G.lit(per), G.lit(area), how))
    out += ["]", ""]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print("wrote %d polygons to %s" % (len(rows), OUT))


if __name__ == "__main__":
    main()
