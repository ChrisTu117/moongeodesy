"""Independent high-precision reference for the WMM2025 field.

This is a second, deliberately different implementation of the model, used
only to generate test vectors (src/magnetic/wmm2025_independent_data_test.mbt):

  * Legendre functions from the exact rational coefficients of
    d^m/dx^m P_n(x) (Rodrigues' formula, Python Fraction), not a recurrence;
  * the field as minus the gradient of the potential in Earth-centred
    Cartesian coordinates, by mpmath numerical differentiation at 40 digits,
    not analytic derivative sums;
  * the north/east/down components by projecting that vector on the
    geodetic unit vectors, not by rotating through the angle phi' - phi.

The MoonBit code (src/magnetic/magnetic.mbt) shares none of these steps. Both
take only the coefficients and WGS-84 constants from the official model.

The points are the 100 official test inputs plus poles, the ends of the
height and date ranges, and Chinese cities.

Needs: pip install mpmath
Usage: python3 tools/gen_wmm_independent.py WMM.COF WMM2025_TestValues.txt
(the two files are described in tools/gen_wmm.py; the SHA-256 is checked)
"""
import hashlib
import sys
from fractions import Fraction
from math import comb, factorial

import mpmath as mp

mp.mp.dps = 40
NMAX = 12
OUT = "src/magnetic/wmm2025_independent_data_test.mbt"
SHA256 = {
    "cof": "dfa8597825af4e0b87ff4198a5b4fb661b3c49f4cd090cd0164e0259b075582f",
    "precise": "e6975b093dddeb6153e0b23cc418425c438167e7c5b1dd795da379cb654f5819",
}


def read_checked(path, key):
    raw = open(path, "rb").read()
    if hashlib.sha256(raw).hexdigest() != SHA256[key]:
        sys.exit(f"{path}: unexpected sha256")
    return raw.decode("ascii")


def read_cof(text):
    g, h, gd, hd = {}, {}, {}, {}
    lines = text.splitlines()
    epoch = mp.mpf(lines[0].split()[0])
    for line in lines[1:]:
        f = line.split()
        if not f or f[0].startswith("9999"):
            break
        n, m = int(f[0]), int(f[1])
        g[n, m], h[n, m], gd[n, m], hd[n, m] = (mp.mpf(x) for x in f[2:6])
    return epoch, g, h, gd, hd


def dpoly(n, m):
    """Exact coefficients {power: Fraction} of d^m/dx^m P_n(x)."""
    c = {}
    for k in range(n // 2 + 1):
        c[n - 2 * k] = Fraction((-1) ** k * comb(n, k) * comb(2 * n - 2 * k, n), 2**n)
    for _ in range(m):
        c = {p - 1: v * p for p, v in c.items() if p > 0}
    return c


POLY = {(n, m): dpoly(n, m) for n in range(NMAX + 1) for m in range(n + 1)}


def norm(n, m):
    if m == 0:
        return mp.mpf(1)
    return mp.sqrt(mp.mpf(2) * factorial(n - m) / factorial(n + m))


NORM = {(n, m): norm(n, m) for n in range(NMAX + 1) for m in range(n + 1)}


def schmidt(n, m, x):
    """Schmidt semi-normalized P(n, m)(x), no Condon-Shortley phase."""
    s = mp.fsum(mp.mpf(v.numerator) / v.denominator * x**p for p, v in POLY[n, m].items())
    return NORM[n, m] * (1 - x * x) ** (mp.mpf(m) / 2) * s


class Model:
    def __init__(self, text):
        self.epoch, self.g, self.h, self.gd, self.hd = read_cof(text)
        self.a = mp.mpf(6371200)

    def potential(self, gc, hc, x, y, z):
        r = mp.sqrt(x * x + y * y + z * z)
        sinp = z / r
        lam = mp.atan2(y, x)
        v = mp.mpf(0)
        for n in range(1, NMAX + 1):
            s = mp.mpf(0)
            for m in range(n + 1):
                s += (gc[n, m] * mp.cos(m * lam) + hc[n, m] * mp.sin(m * lam)) * schmidt(n, m, sinp)
            v += (self.a / r) ** (n + 1) * s
        return self.a * v

    def ned(self, gc, hc, lat, lon, hm):
        A = mp.mpf(6378137)
        f = 1 / mp.mpf("298.257223563")
        e2 = f * (2 - f)
        phi, lam = mp.radians(lat), mp.radians(lon)
        sp, cp = mp.sin(phi), mp.cos(phi)
        if abs(lat) == 90:
            cp = mp.mpf(0)
        sl, cl = mp.sin(lam), mp.cos(lam)
        N = A / mp.sqrt(1 - e2 * sp**2)
        X, Y, Z = (N + hm) * cp * cl, (N + hm) * cp * sl, (N * (1 - e2) + hm) * sp
        # B = -grad V by central differences with a fixed 1 mm step. (The
        # default step is so small at 40 digits that, on the polar axis, the
        # displaced point is still on the axis and the horizontal field is
        # lost.) Truncation error is about 1e-16 nT.
        step = mp.mpf("1e-3")
        bx = -mp.diff(lambda q: self.potential(gc, hc, q, Y, Z), X, h=step)
        by = -mp.diff(lambda q: self.potential(gc, hc, X, q, Z), Y, h=step)
        bz = -mp.diff(lambda q: self.potential(gc, hc, X, Y, q), Z, h=step)
        north = (-sp * cl, -sp * sl, cp)
        east = (-sl, cl, mp.mpf(0))
        down = (-cp * cl, -cp * sl, -sp)

        def dot(u):
            return u[0] * bx + u[1] * by + u[2] * bz

        return dot(north), dot(east), dot(down)

    def field(self, year, hkm, lat, lon):
        dt = mp.mpf(year) - self.epoch
        gc = {k: self.g[k] + dt * self.gd[k] for k in self.g}
        hc = {k: self.h[k] + dt * self.hd[k] for k in self.h}
        hm = mp.mpf(hkm) * 1000
        return self.ned(gc, hc, lat, lon, hm) + self.ned(self.gd, self.hd, lat, lon, hm)


def lit(v):
    r = repr(float(v))
    if "e" in r:
        m, e = r.split("e")
        if "." not in m:
            m += ".0"
        r = m + "e" + e
    return r


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    model = Model(read_checked(sys.argv[1], "cof"))
    points = []
    for line in read_checked(sys.argv[2], "precise").splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            f = [float(x) for x in line.split()]
            points.append((f[0], f[1], f[2], f[3]))
    # Poles (several meridians), ends of the height and date ranges, and
    # Chinese cities (Mohe, Harbin, Beijing, Shanghai, Shenzhen, Sanya,
    # Lhasa, Urumqi, Kashgar).
    for lat in (90.0, -90.0):
        for lon in (0.0, 90.0, -135.0):
            points.append((2026.0, 0.0, lat, lon))
    points += [
        (2025.0, -1.0, 45.0, 45.0),
        (2030.0, 850.0, -45.0, -45.0),
        (2030.0, 0.0, 0.0, 180.0),
        (2029.75, 400.0, 60.0, -150.0),
    ]
    for lat, lon in (
        (52.97, 122.54),
        (45.80, 126.53),
        (39.9087, 116.3975),
        (31.2304, 121.4737),
        (22.5431, 114.0579),
        (18.25, 109.51),
        (29.65, 91.13),
        (43.83, 87.62),
        (39.47, 75.99),
    ):
        points.append((2026.5, 0.0, lat, lon))
    rows = []
    for p in points:
        x, y, z, xd, yd, zd = model.field(*p)
        rows.append(list(p) + [x, y, z, xd, yd, zd])
    out = [
        "// Generated by tools/gen_wmm_independent.py. Do not edit.",
        "// Independent 40-digit evaluation of WMM2025 (exact Legendre polynomials,",
        "// numerical gradient of the potential in ECEF, projection on the geodetic",
        "// frame). Columns: year, height (km), lat, lon, X, Y, Z (nT),",
        "// Xdot, Ydot, Zdot (nT/yr).",
        "",
        "///|",
        "let wmm2025_independent : Array[FixedArray[Double]] = [",
    ]
    for r in rows:
        out.append("  [" + ", ".join(lit(v) for v in r) + "],")
    out.append("]")
    out.append("")
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print(f"wrote {OUT}: {len(rows)} points")


if __name__ == "__main__":
    main()
