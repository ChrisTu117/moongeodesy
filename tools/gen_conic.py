"""Reference data for the Lambert conformal conic and Albers equal-area
projections (src/conic).

Writes
  src/conic/conic_data_test.mbt   high-precision and PROJ vectors
  src/conic/gigs_data_test.mbt    IOGP GIGS test procedures 5102, 5103, 5109

The MoonBit code ports GeographicLib's formulation (Karney's divided
differences). None of the references below use it:

1. A 50-digit mpmath evaluation of the closed-form formulas in IOGP
   Guidance Note 7-2 (EPSG methods 9801, 9802, 9822; the same as Snyder,
   USGS Professional Paper 1395, chapters 14 and 15): easting, northing,
   convergence n * (lon - lon0), and scale k = rho * n / (a * m).
2. PROJ (through pyproj), an independent C implementation of the same
   Snyder formulas in double precision: easting, northing, convergence,
   both scale factors, and the inverse.
3. IOGP GIGS conversion tests 5102 (LCC 1SP), 5103 (LCC 2SP) and 5109
   (Albers): official points with the EPSG parameters of each CRS.
   GIGS data (c) International Association of Oil & Gas Producers, used
   under the GIGS Terms of Use v2.1 (https://gigs.iogp.org/TermsOfUse):
   free redistribution with acknowledgement of IOGP's ownership; not to be
   distributed for profit. Download the ASCII files from
   https://github.com/IOGP-GIGS/GIGSTestDataset (folder
   "GIGS 5100 Conversion test data/ASCII").

Needs: pip install mpmath pyproj
Usage: python3 tools/gen_conic.py <dir with the GIGS 5102/5103/5109 files>
"""
import os
import sys

import mpmath as mp
import pyproj

mp.mp.dps = 50
OUT = "src/conic/conic_data_test.mbt"
GIGS_OUT = "src/conic/gigs_data_test.mbt"

GRS80 = (6378137.0, 1 / 298.257222101)
WGS84 = (6378137.0, 1 / 298.257223563)
KRASOVSKY = (6378245.0, 1 / 298.3)
CLARKE1866 = (6378206.4, 1 / 294.978698213898)
INTL1924 = (6378388.0, 1 / 297.0)
CLARKE1880IGN = (6378249.2, 1 / 293.4660212936269)

# kind: 0 = LCC 2SP, 1 = LCC 1SP, 2 = Albers.
# (name, kind, (a, f), lat0, lon0, stdlat1, stdlat2, k0)
SETUPS = [
    ("China CGCS2000 LCC 25/47", 0, GRS80, 0.0, 105.0, 25.0, 47.0, 1.0),
    ("China CGCS2000 Albers 25/47", 2, GRS80, 0.0, 105.0, 25.0, 47.0, 1.0),
    ("China Krasovsky Albers 25/47 110E", 2, KRASOVSKY, 0.0, 110.0, 25.0, 47.0, 1.0),
    ("US CONUS Albers", 2, GRS80, 23.0, -96.0, 29.5, 45.5, 1.0),
    ("US LCC 33/45", 0, GRS80, 39.0, -96.0, 33.0, 45.0, 1.0),
    ("Europe LCC 35/65", 0, GRS80, 52.0, 10.0, 35.0, 65.0, 1.0),
    ("Australian Albers -18/-36", 2, GRS80, 0.0, 132.0, -18.0, -36.0, 1.0),
    ("South LCC -30/-60", 0, WGS84, -45.0, -60.0, -30.0, -60.0, 1.0),
    ("LCC 1SP 18N k0 1", 1, CLARKE1866, 18.0, -77.0, 18.0, 18.0, 1.0),
    ("LCC 1SP 35S k0 0.9996", 1, WGS84, -35.0, 150.0, -35.0, -35.0, 0.9996),
    ("LCC close parallels 40/40.0001", 0, WGS84, 40.0, 0.0, 40.0, 40.0001, 1.0),
    ("LCC polar 70/85", 0, WGS84, 75.0, 0.0, 70.0, 85.0, 1.0),
    ("Albers polar 60/85", 2, WGS84, 70.0, 0.0, 60.0, 85.0, 1.0),
    ("LCC near equator 0.5/3", 0, WGS84, 0.0, 0.0, 0.5, 3.0, 1.0),
    ("Albers -5/10", 2, WGS84, 0.0, 30.0, -5.0, 10.0, 1.0),
    ("Albers one parallel 40", 2, WGS84, 40.0, -100.0, 40.0, 40.0, 1.0),
]

LATS = [-89.0, -75.0, -60.0, -45.0, -30.0, -15.0, 0.0, 12.5, 25.0, 36.0, 47.0, 60.0, 75.0, 89.0]
DLONS = [-170.0, -90.0, -30.0, -2.5, 0.0, 7.0, 45.0, 120.0]


def lat_range(kind, s1, s2):
    """Latitudes kept away from the pole opposite the cone's apex."""
    north = s1 + s2 >= 0
    if kind == 2:
        return [x for x in LATS if (x > -89.0 if north else x < 89.0)]
    return [x for x in LATS if (x >= -60.0 if north else x <= 60.0)]


# ---------------------------------------------------------------- mpmath --

def ell(a, f):
    a = mp.mpf(a)
    f = mp.mpf(f)
    e2 = f * (2 - f)
    return a, e2, mp.sqrt(e2)


def m_of(phi, e2):
    return mp.cos(phi) / mp.sqrt(1 - e2 * mp.sin(phi) ** 2)


def t_of(phi, e):
    if phi == mp.pi / 2:
        return mp.mpf(0)
    s = mp.sin(phi)
    return mp.tan(mp.pi / 4 - phi / 2) / ((1 - e * s) / (1 + e * s)) ** (e / 2)


def alpha_of(phi, e2, e):
    s = mp.sin(phi)
    return (1 - e2) * (s / (1 - e2 * s * s) - 1 / (2 * e) * mp.log((1 - e * s) / (1 + e * s)))


def rad(d):
    return mp.mpf(d) * mp.pi / 180


def mp_forward(setup, lat, lon):
    _, kind, (a0, f0), lat0, lon0, s1, s2, k0 = setup
    a, e2, e = ell(a0, f0)
    phi = rad(lat)
    phi0 = rad(lat0)
    p1, p2 = rad(s1), rad(s2)
    dlon = (mp.mpf(lon) - mp.mpf(lon0) + 180) % 360 - 180
    dlam = rad(dlon)
    if kind in (0, 1):
        if kind == 1 or s1 == s2:
            n = mp.sin(p1)
            F = m_of(p1, e2) / (n * t_of(p1, e) ** n)
            scale = mp.mpf(k0)
        else:
            m1, m2 = m_of(p1, e2), m_of(p2, e2)
            t1, t2 = t_of(p1, e), t_of(p2, e)
            n = (mp.log(m1) - mp.log(m2)) / (mp.log(t1) - mp.log(t2))
            F = m1 / (n * t1**n)
            scale = mp.mpf(1)
        # t^n for the southern hemisphere: t > 1 and n < 0, still fine.
        r = a * scale * F * t_of(phi, e) ** n
        r0 = a * scale * F * t_of(phi0, e) ** n
        theta = n * dlam
        x = r * mp.sin(theta)
        y = r0 - r * mp.cos(theta)
        k = r * n / (a * m_of(phi, e2))
    else:
        m1, m2 = m_of(p1, e2), m_of(p2, e2)
        a1, a2 = alpha_of(p1, e2, e), alpha_of(p2, e2, e)
        if s1 == s2:
            n = mp.sin(p1)
        else:
            n = (m1**2 - m2**2) / (a2 - a1)
        C = m1**2 + n * a1
        rho = a * mp.sqrt(C - n * alpha_of(phi, e2, e)) / n
        rho0 = a * mp.sqrt(C - n * alpha_of(phi0, e2, e)) / n
        theta = n * dlam
        x = rho * mp.sin(theta)
        y = rho0 - rho * mp.cos(theta)
        k = rho * n / (a * m_of(phi, e2))
    gamma = theta * 180 / mp.pi
    return x, y, gamma, k


# ------------------------------------------------------------------ PROJ --

def proj_string(setup):
    _, kind, (a, f), lat0, lon0, s1, s2, k0 = setup
    base = f"+a={a!r} +rf={1 / f!r} +x_0=0 +y_0=0 +units=m +no_defs"
    if kind == 0:
        return f"+proj=lcc +lat_1={s1!r} +lat_2={s2!r} +lat_0={lat0!r} +lon_0={lon0!r} {base}"
    if kind == 1:
        return f"+proj=lcc +lat_1={lat0!r} +lat_0={lat0!r} +lon_0={lon0!r} +k_0={k0!r} {base}"
    return f"+proj=aea +lat_1={s1!r} +lat_2={s2!r} +lat_0={lat0!r} +lon_0={lon0!r} {base}"


def lit(v):
    v = float(v)
    if v != v:
        return "@double.not_a_number"
    r = repr(v)
    if "e" in r:
        m, e = r.split("e")
        if "." not in m:
            m += ".0"
        r = m + "e" + e
    return r


def row(values):
    return "  [" + ", ".join(lit(v) for v in values) + "],"


def gen_main():
    setups, ref, prj = [], [], []
    for i, s in enumerate(SETUPS):
        _, kind, (a, f), lat0, lon0, s1, s2, k0 = s
        setups.append([kind, a, f, lat0, lon0, s1, s2, k0])
        p = pyproj.Proj(proj_string(s))
        for lat in lat_range(kind, s1, s2):
            for dl in DLONS:
                lon = lon0 + dl
                lon = (lon + 180.0) % 360.0 - 180.0
                x, y, g, k = mp_forward(s, lat, lon)
                ref.append([i, lat, lon, x, y, g, k])
                px, py = p(lon, lat)
                fac = p.get_factors(lon, lat)
                ilon, ilat = p(px, py, inverse=True)
                prj.append([i, lat, lon, px, py, fac.meridian_convergence,
                            fac.parallel_scale, fac.meridional_scale, ilat, ilon])
    out = [
        "// Generated by tools/gen_conic.py. Do not edit.",
        f"// PROJ {pyproj.proj_version_str} via pyproj {pyproj.__version__}; mpmath {mp.__version__} at 50 digits.",
        "",
        "///|",
        "/// Projections: kind (0 LCC 2SP, 1 LCC 1SP, 2 Albers), a, f, lat0, lon0,",
        "/// stdlat1, stdlat2, k0. Names (same order):",
    ]
    for i, s in enumerate(SETUPS):
        out.append(f"///   {i}: {s[0]}")
    out.append("let conic_setups : Array[FixedArray[Double]] = [")
    out += [row(r) for r in setups]
    out += ["]", "", "///|",
            "/// 50-digit closed-form (IOGP GN 7-2 / Snyder) values: setup, lat, lon,",
            "/// x, y, convergence (deg), k (scale along the parallel).",
            "let conic_mp : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in ref]
    out += ["]", "", "///|",
            "/// PROJ values: setup, lat, lon, x, y, meridian convergence (deg),",
            "/// parallel scale, meridional scale, inverse lat, inverse lon.",
            "let conic_proj : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in prj]
    out += ["]", ""]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print(f"wrote {OUT}: {len(setups)} projections, {len(ref)} points")


# ------------------------------------------------------------------ GIGS --

FT = 0.3048
FTUS = 1200.0 / 3937.0
GRAD = 0.9  # degrees per gradian
# file, kind, (a, f), lat0, lon0, s1, s2, k0, FE (m), FN (m), length unit (m),
# angle unit (deg); EPSG parameters as printed by pyproj for each CRS.
GIGS = [
    ("GIGS_conv_5102_LCC1_output_part1.txt", 1, INTL1924, 46.8, 2.337229166666667,
     46.8, 46.8, 0.99987742, 600000.0, 2200000.0, 1.0, 1.0),
    # NTF (Paris) / Lambert zone II: angles in gradians from the Paris meridian.
    ("GIGS_conv_5102_LCC1_output_part2.txt", 1, CLARKE1880IGN, 52.0 * GRAD, 0.0,
     52.0 * GRAD, 52.0 * GRAD, 0.99987742, 600000.0, 2200000.0, 1.0, GRAD),
    ("GIGS_conv_5103_LCC2_output_part1.txt", 0, INTL1924, 90.0, 4.367486666666666,
     51.16666723333333, 49.8333339, 1.0, 150000.013, 5400088.438, 1.0, 1.0),
    ("GIGS_conv_5103_LCC2_output_part2.txt", 0, GRS80, 40.333333333333336, -111.5,
     41.78333333333333, 40.71666666666667, 1.0, 1640419.948 * FT, 3280839.895 * FT, FT, 1.0),
    ("GIGS_conv_5103_LCC2_output_part3.txt", 0, GRS80, 40.333333333333336, -111.5,
     41.78333333333333, 40.71666666666667, 1.0, 1640416.6667 * FTUS, 3280833.3333 * FTUS,
     FTUS, 1.0),
    ("GIGS_conv_5109_Albers_output.txt", 2, GRS80, 0.0, 132.0, -18.0, -36.0, 1.0,
     0.0, 0.0, 1.0, 1.0),
]


def gen_gigs(folder):
    sets, pts = [], []
    for i, g in enumerate(GIGS):
        fname, kind, (a, f), lat0, lon0, s1, s2, k0, fe, fn, lu, au = g
        sets.append([kind, a, f, lat0, lon0, s1, s2, k0, fe, fn, lu, au])
        tol = None
        for line in open(os.path.join(folder, fname), encoding="utf-8"):
            if line.startswith("# Cartesian Tolerance:"):
                tol = line
            if line.startswith("#") or not line.strip():
                continue
            fld = line.rstrip("\n").split("\t")
            pid = int(fld[0].split("-")[-1])
            fwd = 1.0 if fld[6].strip() == "FORWARD" else 0.0
            pts.append([i, pid, float(fld[1]), float(fld[2]), float(fld[3]), float(fld[4]), fwd])
        if tol is None:
            sys.exit(f"{fname}: no tolerance header")
    out = [
        "// Generated by tools/gen_conic.py. Do not edit.",
        "// IOGP GIGS Test Dataset v2.1.0 (EPSG Dataset 10.019), conversion tests",
        "// 5102 (Lambert Conic Conformal 1SP), 5103 (2SP) and 5109 (Albers).",
        "// GIGS data (c) International Association of Oil & Gas Producers (IOGP),",
        "// used under the GIGS Terms of Use v2.1, https://gigs.iogp.org/TermsOfUse",
        "// (free redistribution with this acknowledgement; not for profit).",
        "// Source: https://github.com/IOGP-GIGS/GIGSTestDataset",
        "",
        "///|",
        "/// CRS: kind (0 LCC 2SP, 1 LCC 1SP, 2 Albers), a, f, lat0, lon0, stdlat1,",
        "/// stdlat2, k0, FE (m), FN (m), length unit (m), angle unit (deg). Files:",
    ]
    for i, g in enumerate(GIGS):
        out.append(f"///   {i}: {g[0]}")
    out.append("let gigs_sets : Array[FixedArray[Double]] = [")
    out += [row(r) for r in sets]
    out += ["]", "", "///|",
            "/// Points: set, point number, lat, lon (angle unit), easting, northing",
            "/// (length unit), 1 = FORWARD / 0 = REVERSE in the GIGS file.",
            "let gigs_points : Array[FixedArray[Double]] = ["]
    out += [row(r) for r in pts]
    out += ["]", ""]
    with open(GIGS_OUT, "w") as fh:
        fh.write("\n".join(out))
    print(f"wrote {GIGS_OUT}: {len(sets)} CRSs, {len(pts)} points")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    gen_main()
    gen_gigs(sys.argv[1])
