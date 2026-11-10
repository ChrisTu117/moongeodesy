"""Geocentric (ECEF) and LocalCartesian (ENU) reference vectors.

Writes src/geodesic/local_cartesian_data_test.mbt.

The MoonBit code ports GeographicLib's Geocentric and LocalCartesian
(Vermeille's closed-form reverse, rotation matrices). The references here do
not use that code:

1. 50-digit mpmath values from the textbook definitions:
   forward  X = (N + h) cos(lat) cos(lon), Y = (N + h) cos(lat) sin(lon),
            Z = (N (1 - e2) + h) sin(lat), N = a / sqrt(1 - e2 sin^2 lat);
   reverse  by fixed-point iteration of tan(lat) = (Z + e2 N sin(lat)) / p
            until it no longer changes at 50 digits;
   ENU      the ECEF difference from the origin rotated by the standard
            east/north/up unit vectors.
2. PROJ (through pyproj): +proj=cart for ECEF and +proj=topocentric for ENU,
   on WGS-84, GRS80, and a sphere.

Needs: pip install mpmath pyproj
Usage: python3 tools/gen_local_cartesian.py
"""
import mpmath as mp
import pyproj

mp.mp.dps = 50
OUT = "src/geodesic/local_cartesian_data_test.mbt"

# tag: (a, f, PROJ ellipsoid args or None)
ELLIPSOIDS = {
    0: ("6378137", "1/298.257223563", "+ellps=WGS84"),
    1: ("6371000", "0", "+a=6371000 +b=6371000"),
    2: ("6378137", "-1/298.257223563", None),  # prolate: not in PROJ
    3: ("6378137", "1/298.257222101", "+ellps=GRS80"),
}


def ell(tag):
    a, f, _ = ELLIPSOIDS[tag]
    a = mp.mpf(a)
    num, _, den = f.partition("/")
    f = mp.mpf(num) / mp.mpf(den) if den else mp.mpf(num)
    return a, f * (2 - f)


def rad(d):
    return mp.mpf(d) * mp.pi / 180


def sincos_deg(d):
    # exact zeros at multiples of 90 degrees
    d = mp.mpf(d)
    if d % 180 == 0:
        return mp.mpf(0), mp.mpf(-1) if d % 360 else mp.mpf(1)
    if d % 180 == 90:
        return (mp.mpf(1) if d % 360 == 90 else mp.mpf(-1)), mp.mpf(0)
    r = rad(d)
    return mp.sin(r), mp.cos(r)


def mp_forward(tag, lat, lon, h):
    a, e2 = ell(tag)
    sp, cp = sincos_deg(lat)
    sl, cl = sincos_deg(lon)
    n = a / mp.sqrt(1 - e2 * sp * sp)
    h = mp.mpf(h)
    return (n + h) * cp * cl, (n + h) * cp * sl, (n * (1 - e2) + h) * sp


def mp_reverse(tag, x, y, z):
    a, e2 = ell(tag)
    x, y, z = mp.mpf(x), mp.mpf(y), mp.mpf(z)
    p = mp.sqrt(x * x + y * y)
    lon = mp.atan2(y, x) * 180 / mp.pi if p != 0 else mp.mpf(0)
    if p == 0:
        b = a * mp.sqrt(1 - e2)
        return (90 if z >= 0 else -90), lon, abs(z) - b
    phi = mp.atan2(z, p * (1 - e2))
    for _ in range(400):
        s = mp.sin(phi)
        n = a / mp.sqrt(1 - e2 * s * s)
        new = mp.atan2(z + e2 * n * s, p)
        if abs(new - phi) < mp.mpf(10) ** -45:
            phi = new
            break
        phi = new
    s, c = mp.sin(phi), mp.cos(phi)
    n = a / mp.sqrt(1 - e2 * s * s)
    h = p / c - n if abs(c) > abs(s) else z / s - n * (1 - e2)
    return phi * 180 / mp.pi, lon, h


def mp_enu(tag, lat0, lon0, h0, lat, lon, h):
    x0, y0, z0 = mp_forward(tag, lat0, lon0, h0)
    x, y, z = mp_forward(tag, lat, lon, h)
    dx, dy, dz = x - x0, y - y0, z - z0
    sp, cp = sincos_deg(lat0)
    sl, cl = sincos_deg(lon0)
    e = -sl * dx + cl * dy
    n = -sp * cl * dx - sp * sl * dy + cp * dz
    u = cp * cl * dx + cp * sl * dy + sp * dz
    return e, n, u


def lit(v):
    v = float(v)
    r = repr(v)
    if "e" in r:
        m, e = r.split("e")
        if "." not in m:
            m += ".0"
        r = m + "e" + e
    elif "." not in r:
        r += ".0"
    return r


def row(values):
    return "  [" + ", ".join(lit(v) for v in values) + "],"


LATS = [-90.0, -89.9999, -60.0, -45.0, -1.0e-9, 0.0, 30.0, 45.0, 60.0, 89.99, 90.0]
LONS = [-180.0, -120.0, 0.0, 45.5, 179.5]
HS = [-5000.0, 0.0, 8848.86, 400000.0, 35786000.0]

ORIGINS = [
    (39.9087, 116.3975, 0.0),
    (0.0, 0.0, 0.0),
    (-33.8688, 151.2093, 58.0),
    (89.9, -30.0, 2000.0),
    (-90.0, 45.0, 100.0),
    (22.5431, 114.0579, -10.0),
]
TARGETS = [
    (1.0e-4, 1.0e-4, 50.0),
    (-0.01, 0.02, -30.0),
    (1.0, -2.0, 1000.0),
    (10.0, 30.0, 0.0),
    (-45.0, 120.0, 20000.0),
]


def clamp_lat(x):
    return max(-90.0, min(90.0, x))


def main():
    gfwd, grev, gproj, lfwd, lproj = [], [], [], [], []
    for tag in sorted(ELLIPSOIDS):
        prj = ELLIPSOIDS[tag][2]
        cart = None
        if prj:
            cart = pyproj.Transformer.from_pipeline(f"+proj=cart {prj}")
        for i, lat in enumerate(LATS):
            for j, lon in enumerate(LONS):
                # three heights per point, rotating through HS
                for k in range(3):
                    h = HS[(i + j + k) % len(HS)]
                    x, y, z = mp_forward(tag, lat, lon, h)
                    gfwd.append([tag, lat, lon, h, x, y, z])
                    xr, yr, zr = float(x), float(y), float(z)
                    rl, rlon, rh = mp_reverse(tag, xr, yr, zr)
                    grev.append([tag, xr, yr, zr, rl, rlon, rh])
                    if cart is not None:
                        px, py, pz = cart.transform(lon, lat, h)
                        gproj.append([tag, lat, lon, h, px, py, pz])
        for lat0, lon0, h0 in ORIGINS:
            topo = None
            if prj:
                topo = pyproj.Transformer.from_pipeline(
                    f"+proj=pipeline +step +proj=cart {prj} "
                    f"+step +proj=topocentric {prj} +lat_0={lat0!r} +lon_0={lon0!r} +h_0={h0!r}"
                )
            for dlat, dlon, h in TARGETS:
                lat = clamp_lat(lat0 + dlat)
                lon = lon0 + dlon
                e, n, u = mp_enu(tag, lat0, lon0, h0, lat, lon, h)
                lfwd.append([tag, lat0, lon0, h0, lat, lon, h, e, n, u])
                if topo is not None:
                    pe, pn, pu = topo.transform(lon, lat, h)
                    lproj.append([tag, lat0, lon0, h0, lat, lon, h, pe, pn, pu])
    out = [
        "// Generated by tools/gen_local_cartesian.py. Do not edit.",
        f"// 50-digit mpmath values from the textbook definitions, and PROJ {pyproj.proj_version_str}",
        "// (+proj=cart, +proj=topocentric) via pyproj. Neither uses GeographicLib.",
        "// Ellipsoid tag: 0 WGS-84, 1 sphere a = 6371000, 2 prolate f = -1/298.257223563,",
        "// 3 GRS80 (PROJ has no prolate ellipsoid).",
        "",
        "///|",
        "/// Geocentric forward, 50 digits: tag, lat, lon, h, X, Y, Z.",
        "let geocentric_forward : Array[FixedArray[Double]] = [",
    ]
    out += [row(r) for r in gfwd]
    out += ["]", "", "///|",
            "/// Geocentric reverse, 50 digits, of the double X, Y, Z: tag, X, Y, Z, lat, lon, h.",
            "let geocentric_reverse : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in grev]
    out += ["]", "", "///|", "/// Geocentric forward by PROJ +proj=cart: tag, lat, lon, h, X, Y, Z.",
            "let geocentric_proj : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in gproj]
    out += ["]", "", "///|",
            "/// LocalCartesian forward, 50 digits: tag, lat0, lon0, h0, lat, lon, h, x, y, z.",
            "let local_cartesian_forward : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in lfwd]
    out += ["]", "", "///|",
            "/// LocalCartesian forward by PROJ +proj=topocentric: tag, lat0, lon0, h0, lat, lon, h, x, y, z.",
            "let local_cartesian_proj : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in lproj]
    out += ["]", ""]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print(f"wrote {OUT}: {len(gfwd)} forward, {len(grev)} reverse, {len(gproj)} PROJ ECEF, "
          f"{len(lfwd)} ENU, {len(lproj)} PROJ ENU")


if __name__ == "__main__":
    main()
