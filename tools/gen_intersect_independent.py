"""Geodesic-intersection references with known answers (mpmath, 45 digits).

tools/gen_intersect.py uses IntersectTool, the code src/geodesic/intersect.mbt
is ported from. These cases have answers that follow from geometry, with the
lengths from tools/mpgeod.py (meridian arcs and the auxiliary-sphere integrals
by quadrature), not from GeographicLib:

* the equator meets a meridian at the meridian's longitude: x is the
  equator arc a * dlon, y the meridian arc from the start latitude
  (also with an offset that picks the crossing on the far side of the earth);
* two meridians meet at a pole: x and y are the meridian arcs to it;
* a geodesic meets the equator at its nodes, the points at auxiliary arc
  sigma = k pi; distance and longitude come from the two classical integrals
  (Intersect::closest with the equator, and Intersect::next from a node:
  the next node half a period away, on WGS-84 only; on the prolate
  ellipsoid nodes are more than 180 deg apart and the candidates come
  within a few km of each other in |x| + |y|);
* two geodesics that are mirror images in a meridian, or in the equator,
  meet on the mirror at equal distances x = y; the distance is found by
  solving the direct problem for the longitude of the mirror (or the node);

Ellipsoid code: 0 WGS-84, 2 prolate f = -1/298.257223563 (a = 6378137 m), as
in intersect_data_test.mbt. For each closest case the generator checks that
every other crossing it can see is at least 10 km farther in |x| + |y|.

Usage (from the repository root):
  pip install mpmath
  python3 tools/gen_intersect_independent.py && moon fmt
"""
import mpmath as mp

import mpgeod as G

OUT = "src/geodesic/intersect_independent_data_test.mbt"
ELL = {0: (6378137.0, 1 / 298.257223563), 2: (6378137.0, -1 / 298.257223563)}
MARGIN = 10000  # metres between the answer and the next-closest crossing


def arc(E, lat):
    """Signed meridian arc from the equator to lat (degrees)."""
    return E.meridian_arc(G.rad(G.D(lat)))


def eq_arc(E, dlon):
    return E.a * G.rad(G.D(dlon))


def closest_check(cands, x0=0, y0=0):
    """The first candidate must be the closest by a margin."""
    l1 = [abs(x - x0) + abs(y - y0) for x, y in cands]
    best = l1[0]
    assert all(v > best + MARGIN for v in l1[1:]), (cands, l1)


def equator_meridian(E, lonX, aziX, latY, lonY, aziY):
    """X along the equator from lonX (aziX = +-90), Y along meridian lonY from
    latY (aziY = 0 or 180). Returns x, y of the crossing at lonY."""
    sx = 1 if aziX == 90 else -1
    sy = 1 if aziY == 180 else -1     # +1: heading south
    Q = E.quarter_meridian()
    x = sx * eq_arc(E, G.ang_diff(G.D(lonX), G.D(lonY)))
    y = sy * arc(E, latY)
    # Other crossing: lonY + 180, past a pole, either way along X.
    far = G.ang_diff(G.D(lonX), G.D(lonY) + 180)
    yfar = [y + 2 * Q * sgn for sgn in (1, -1)]
    xfar = [sx * eq_arc(E, far), sx * eq_arc(E, far - 360 if far > 0 else far + 360)]
    closest_check([(x, y)] + [(a, b) for a in xfar for b in yfar]
                  + [(x + 2 * mp.pi * E.a * k, y) for k in (-1, 1)])
    return x, y


def meridians_at_pole(E, latX, lonX, aziX, latY, lonY, aziY):
    """Two meridians (azimuth 0 or 180) meet at the poles."""
    Q = E.quarter_meridian()

    def to_pole(lat, azi, north):
        m = arc(E, lat)
        d = (Q - m) if north else (Q + m)       # distance to that pole
        heading_north = azi == 0
        return d if heading_north == north else -d

    cands = []
    for north in (True, False):
        cands.append((to_pole(latX, aziX, north), to_pole(latY, aziY, north)))
    cands.sort(key=lambda c: abs(c[0]) + abs(c[1]))
    closest_check(cands)
    return cands[0]


def equator_nodes(E, lat1, lon1, azi1, lonE):
    """X general, Y the equator eastward from lonE: the closest node of X."""
    nodes = E.geod_nodes(lat1, lon1, azi1, range(-3, 4))
    cands = []
    for _, s, lon in nodes:
        dl = G.ang_diff(G.D(lonE), lon)
        for k in (-1, 0, 1):
            cands.append((s, eq_arc(E, dl + 360 * k)))
    cands.sort(key=lambda c: abs(c[0]) + abs(c[1]))
    closest_check(cands)
    return cands[0]


def mirror_meridian(E, lat, lon0, d, azi):
    """X from (lat, lon0 - d, azi), Y from (lat, lon0 + d, -azi): they meet on
    meridian lon0 at x = y."""
    lonX = lon0 - d
    s_guess = E.a * G.rad(G.D(d)) / mp.sin(G.rad(G.D(azi))) / mp.cos(G.rad(G.D(lat)))
    s = E.geod_s_at_lon(lat, lonX, azi, lon0, s_guess)
    return s, s


def mirror_equator(E, lat, lon, azi):
    """X from (lat, lon, azi), Y from (-lat, lon, 180 - azi): they meet where X
    crosses the equator, at x = y."""
    nodes = E.geod_nodes(lat, lon, azi, range(-2, 3))
    ahead = sorted((s for _, s, _ in nodes if s > 0))
    behind = sorted((-s for _, s, _ in nodes if s < 0))
    s = ahead[0]
    assert s + MARGIN < min(ahead[1], behind[0]), (ahead, behind)
    return s, s


def node_next(E, lon0, aziY):
    """X the equator eastward from (0, lon0), Y from the same node with aziY
    (0 < aziY < 180). Candidates: the other nodes of Y, and the start point
    reached again after going round the equator (x = +-2 pi a, y = 0), each
    with x wrapped by +-2 pi a. A candidate and its mirror (-x, -y) are
    equally close; returns |x|, |y| of the closest pair."""
    nodes = E.geod_nodes(0, lon0, aziY, range(-3, 5))
    circ = 2 * mp.pi * E.a
    cands = []
    for _, s, lon in nodes:
        dx = E.a * G.rad(lon - G.D(lon0))
        for j in (-1, 0, 1):
            if abs(s) < 1 and j == 0:
                continue                      # the starting point itself
            cands.append((dx + circ * j, s))
    cands.sort(key=lambda c: abs(c[0]) + abs(c[1]))
    best = abs(cands[0][0]) + abs(cands[0][1])
    tied = [c for c in cands if abs(abs(c[0]) + abs(c[1]) - best) < 1e-20]
    assert all(abs(abs(c[0]) - abs(tied[0][0])) < 1e-20 for c in tied), tied
    rest = [abs(c[0]) + abs(c[1]) for c in cands if c not in tied]
    assert min(rest) > best + MARGIN, (cands[:4], best, min(rest))
    return abs(tied[0][0]), abs(tied[0][1])


def main():
    ells = {c: G.Ellipsoid(*ae) for c, ae in ELL.items()}
    closest = []
    nxt = []
    seg = []
    for code, E in ells.items():
        Q = E.quarter_meridian()
        for lonX, aziX, latY, lonY, aziY in ((10.0, 90.0, 40.0, 30.0, 180.0),
                                             (-170.0, 90.0, -35.0, 170.0, 0.0),
                                             (0.0, -90.0, 60.0, 5.0, 180.0),
                                             (100.0, 90.0, 89.5, 101.0, 180.0),
                                             (0.0, 90.0, -10.0, -30.0, 180.0)):
            x, y = equator_meridian(E, lonX, aziX, latY, lonY, aziY)
            closest.append([code, 0.0, lonX, aziX, latY, lonY, aziY, 0, 0, x, y])
        # The meridian as X and the equator as Y.
        x, y = equator_meridian(E, 20.0, 90.0, 45.0, 50.0, 180.0)
        closest.append([code, 45.0, 50.0, 180.0, 0.0, 20.0, 90.0, 0, 0, y, x])
        # Offset (x0, y0) at the far crossing: (0, -150), past the south pole.
        x = eq_arc(E, 210)
        y = arc(E, 40) + 2 * Q
        closest.append([code, 0.0, 0.0, 90.0, 40.0, 30.0, 180.0, x, y, x, y])
        closest.append([code, 0.0, 0.0, 90.0, 40.0, 30.0, 180.0, -eq_arc(E, 150), y,
                        -eq_arc(E, 150), y])
        for c in ((30.0, 10.0, 0.0, 50.0, 70.0, 0.0), (30.0, 10.0, 0.0, 50.0, 70.0, 180.0),
                  (-20.0, -100.0, 0.0, 75.0, 30.0, 0.0), (-30.0, 10.0, 180.0, -50.0, 70.0, 180.0)):
            x, y = meridians_at_pole(E, *c)
            closest.append([code, *c, 0, 0, x, y])
        for lat1, lon1, azi1 in ((30.0, 20.0, 120.0), (-50.0, -170.0, 60.0), (10.0, 179.0, 45.0),
                                 (70.0, -60.0, 100.0)):
            (_, s, lon), = E.geod_nodes(lat1, lon1, azi1, [1])
            lonE = float(G.wrap180(lon - 3))
            x, y = equator_nodes(E, lat1, lon1, azi1, lonE)
            closest.append([code, lat1, lon1, azi1, 0.0, lonE, 90.0, 0, 0, x, y])
        for lat, lon0, d, azi in ((20.0, 100.0, 10.0, 45.0), (-40.0, 180.0, 5.0, 120.0),
                                  (60.0, 0.0, 5.0, 80.0)):
            x, y = mirror_meridian(E, lat, lon0, d, azi)
            closest.append([code, lat, lon0 - d, azi, lat, float(G.wrap180(G.D(lon0 + d))), -azi,
                            0, 0, x, y])
        for lat, lon, azi in ((25.0, 40.0, 130.0), (-10.0, -120.0, 20.0)):
            x, y = mirror_equator(E, lat, lon, azi)
            closest.append([code, lat, lon, azi, -lat, lon, 180.0 - azi, 0, 0, x, y])
        for lon0, aziY in ((0.0, 30.0), (-20.0, 60.0), (170.0, 89.0), (45.0, 135.0)):
            if code == 2:
                # Prolate: nodes are more than 180 deg apart, and the next node
                # by the other way round the equator is within 1-8 km in
                # |x| + |y| of the start point after a full turn: no clear answer.
                continue
            x, y = node_next(E, lon0, aziY)
            nxt.append([code, 0.0, lon0, 90.0, aziY, x, y])
        # Segments: equator against meridian.
        seg.append([code, 0.0, 0.0, 0.0, 40.0, 20.0, 25.0, -20.0, 25.0,
                    eq_arc(E, 25), arc(E, 20), 0])
        seg.append([code, 0.0, -10.0, 0.0, 10.0, -30.0, 5.0, 10.0, 5.0,
                    eq_arc(E, 15), arc(E, 30), 0])
        # Past both ends: segmode 3 * 1 + 1.
        seg.append([code, 0.0, 0.0, 0.0, 40.0, 20.0, 50.0, 10.0, 50.0,
                    eq_arc(E, 50), arc(E, 20), 4])
        # Mirror images in a meridian.
        for la1, lo1, la2, lo2, lon0 in ((10.0, 90.0, 30.0, 110.0, 100.0),
                                         (-40.0, 170.0, -20.0, -170.0, 180.0)):
            v = E.geod_inverse(la1, lo1, la2, lo2)
            s = E.geod_s_at_lon(la1, lo1, v["azi1"], lon0, v["s12"] / 2)
            mlo1 = float(G.wrap180(2 * G.D(lon0) - G.D(lo1)))
            mlo2 = float(G.wrap180(2 * G.D(lon0) - G.D(lo2)))
            seg.append([code, la1, lo1, la2, lo2, la1, mlo1, la2, mlo2, s, s, 0])
    out = [
        "// Generated by tools/gen_intersect_independent.py. Do not edit.",
        "// Reference: cases with known answers (equator and meridians, poles, nodes,",
        "// mirror images), lengths by mpmath quadrature at 45 digits; not GeographicLib.",
        "// Ellipsoid code: 0 WGS-84, 2 prolate f = -1/298.257223563. Lengths in metres.",
        "// ix_indep_closest: code latX lonX aziX latY lonY aziY x0 y0 x y",
        "// ix_indep_next: code lat lon aziX aziY |x| |y| (x, y and -x, -y are equally close)",
        "// ix_indep_segment: code latX1 lonX1 latX2 lonX2 latY1 lonY1 latY2 lonY2 x y segmode",
        "",
    ]
    for name, rows in (("ix_indep_closest", closest), ("ix_indep_next", nxt),
                       ("ix_indep_segment", seg)):
        out += ["///|", "let %s : Array[FixedArray[Double]] = [" % name]
        for r in rows:
            out.append("  [" + ", ".join(G.lit(v) for v in r) + "],")
        out += ["]", ""]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print("wrote %d closest, %d next, %d segment rows to %s"
          % (len(closest), len(nxt), len(seg), OUT))


if __name__ == "__main__":
    main()
