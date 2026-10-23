"""LocalCartesian / Geocentric reference vectors.

Faithful transcription of GeographicLib Geocentric.cpp and LocalCartesian.cpp
(Charles Karney, MIT/X11). Cross-checked against PROJ topocentric on WGS-84.
The MoonBit port is tested against these vectors.

Usage:
  python3 tools/gen_local_cartesian.py
"""
from __future__ import annotations

import math
import sys

OUT = "src/geodesic/local_cartesian_data_test.mbt"

# Degree helpers matching GeographicLib::Math (reduce angle, then sin/cos).
def ang_normalize(x: float) -> float:
    y = math.fmod(x, 360.0)
    if y <= -180.0:
        y += 360.0
    elif y > 180.0:
        y -= 360.0
    # GeographicLib AngNormalize uses remainder; match Python geomath style
    # Use Math::AngNormalize: remainder(x, 360) then map -180 to +180
    y = math.remainder(x, 360.0)
    if y == -180.0:
        return 180.0
    return y


def lat_fix(x: float) -> float:
    if math.isnan(x):
        return x
    return max(-90.0, min(90.0, x))


def sincosd(x: float) -> tuple[float, float]:
    # Reduce to [-45, 45] by quarters like GeographicLib for accuracy near 90.
    r = math.remainder(x, 360.0)
    q = int(math.floor(r / 90.0 + 0.5))
    r -= 90.0 * q
    s = math.sin(math.radians(r))
    c = math.cos(math.radians(r))
    q &= 3
    if q == 0:
        return s, c
    if q == 1:
        return c, -s
    if q == 2:
        return -s, -c
    return -c, s


def atan2d(y: float, x: float) -> float:
    return math.degrees(math.atan2(y, x))


class Geocentric:
    def __init__(self, a: float, f: float):
        self.a = a
        self.f = f
        self.e2 = f * (2.0 - f)
        self.e2m = (1.0 - f) ** 2
        self.e2a = abs(self.e2)
        self.e4a = self.e2 * self.e2
        self.maxrad = 2.0 * a / sys.float_info.epsilon

    @staticmethod
    def WGS84() -> "Geocentric":
        return Geocentric(6378137.0, 1.0 / 298.257223563)

    @staticmethod
    def rotation(sphi: float, cphi: float, slam: float, clam: float) -> list[float]:
        M = [0.0] * 9
        M[0] = -slam
        M[3] = clam
        M[6] = 0.0
        M[1] = -clam * sphi
        M[4] = -slam * sphi
        M[7] = cphi
        M[2] = clam * cphi
        M[5] = slam * cphi
        M[8] = sphi
        return M

    def forward(self, lat: float, lon: float, h: float) -> tuple[float, float, float]:
        sphi, cphi = sincosd(lat_fix(lat))
        slam, clam = sincosd(lon)
        n = self.a / math.sqrt(1.0 - self.e2 * sphi * sphi)
        Z = (self.e2m * n + h) * sphi
        X = (n + h) * cphi
        Y = X * slam
        X = X * clam
        return X, Y, Z

    def reverse(self, X: float, Y: float, Z: float) -> tuple[float, float, float]:
        R = math.hypot(X, Y)
        slam = Y / R if R != 0.0 else 0.0
        clam = X / R if R != 0.0 else 1.0
        h = math.hypot(R, Z)
        if h > self.maxrad:
            R = math.hypot(X / 2.0, Y / 2.0)
            slam = (Y / 2.0) / R if R != 0.0 else 0.0
            clam = (X / 2.0) / R if R != 0.0 else 1.0
            H = math.hypot(Z / 2.0, R)
            sphi = (Z / 2.0) / H
            cphi = R / H
        elif self.e4a == 0.0:
            H = math.hypot(1.0 if h == 0.0 else Z, R)
            sphi = (1.0 if h == 0.0 else Z) / H
            cphi = R / H
            h = h - self.a
        else:
            p = (R / self.a) ** 2
            q = self.e2m * (Z / self.a) ** 2
            r = (p + q - self.e4a) / 6.0
            f = self.f
            if f < 0.0:
                p, q = q, p
            if not (self.e4a * q == 0.0 and r <= 0.0):
                S = self.e4a * p * q / 4.0
                r2 = r * r
                r3 = r * r2
                disc = S * (2.0 * r3 + S)
                u = r
                if disc >= 0.0:
                    T3 = S + r3
                    T3 += -math.sqrt(disc) if T3 < 0.0 else math.sqrt(disc)
                    T = math.copysign(abs(T3) ** (1.0 / 3.0), T3)
                    u += T + (r2 / T if T != 0.0 else 0.0)
                else:
                    ang = math.atan2(math.sqrt(-disc), -(S + r3))
                    u += 2.0 * r * math.cos(ang / 3.0)
                v = math.sqrt(u * u + self.e4a * q)
                uv = self.e4a * q / (v - u) if u < 0.0 else u + v
                w = max(0.0, self.e2a * (uv - q) / (2.0 * v))
                k = uv / (math.sqrt(uv + w * w) + w)
                k1 = k if f >= 0.0 else k - self.e2
                k2 = k + self.e2 if f >= 0.0 else k
                d = k1 * R / k2
                H = math.hypot(Z / k1, R / k2)
                sphi = (Z / k1) / H
                cphi = (R / k2) / H
                h = (1.0 - self.e2m / k1) * math.hypot(d, Z)
            else:
                zz = math.sqrt((self.e4a - p if f >= 0.0 else p) / self.e2m)
                xx = math.sqrt(self.e4a - p if f < 0.0 else p)
                H = math.hypot(zz, xx)
                sphi = zz / H
                cphi = xx / H
                if Z < 0.0:
                    sphi = -sphi
                h = -self.a * (self.e2m if f >= 0.0 else 1.0) * H / self.e2a
        lat = atan2d(sphi, cphi)
        lon = atan2d(slam, clam)
        return lat, lon, h


class LocalCartesian:
    def __init__(self, lat0: float, lon0: float, h0: float = 0.0, earth: Geocentric | None = None):
        self.earth = earth or Geocentric.WGS84()
        self.reset(lat0, lon0, h0)

    def reset(self, lat0: float, lon0: float, h0: float = 0.0) -> None:
        self.lat0 = lat_fix(lat0)
        self.lon0 = ang_normalize(lon0)
        self.h0 = h0
        self.x0, self.y0, self.z0 = self.earth.forward(self.lat0, self.lon0, self.h0)
        sphi, cphi = sincosd(self.lat0)
        slam, clam = sincosd(self.lon0)
        self.r = Geocentric.rotation(sphi, cphi, slam, clam)

    def forward(self, lat: float, lon: float, h: float) -> tuple[float, float, float]:
        xc, yc, zc = self.earth.forward(lat, lon, h)
        xc -= self.x0
        yc -= self.y0
        zc -= self.z0
        x = self.r[0] * xc + self.r[3] * yc + self.r[6] * zc
        y = self.r[1] * xc + self.r[4] * yc + self.r[7] * zc
        z = self.r[2] * xc + self.r[5] * yc + self.r[8] * zc
        return x, y, z

    def reverse(self, x: float, y: float, z: float) -> tuple[float, float, float]:
        xc = self.x0 + self.r[0] * x + self.r[1] * y + self.r[2] * z
        yc = self.y0 + self.r[3] * x + self.r[4] * y + self.r[5] * z
        zc = self.z0 + self.r[6] * x + self.r[7] * y + self.r[8] * z
        return self.earth.reverse(xc, yc, zc)


def lit(v: float) -> str:
    if v != v:
        return "@double.not_a_number"
    r = repr(float(v))
    if r == "inf":
        return "@double.infinity"
    if r == "-inf":
        return "@double.neg_infinity"
    if "e" in r or "E" in r:
        m, e = r.lower().split("e")
        if "." not in m:
            m += ".0"
        return m + "e" + e
    if "." not in r:
        return r + ".0"
    return r


def main() -> None:
    # Cross-check vs PROJ when available
    try:
        from pyproj import Transformer

        pipe = (
            "+proj=pipeline "
            "+step +proj=latlong +ellps=WGS84 "
            "+step +proj=cart +ellps=WGS84 "
            "+step +proj=topocentric +ellps=WGS84 +lat_0=39.9087 +lon_0=116.3975 +h_0=0"
        )
        t = Transformer.from_pipeline(pipe)
        lc = LocalCartesian(39.9087, 116.3975, 0.0)
        ref = lc.forward(39.9097, 116.3985, 50.0)
        proj = t.transform(116.3985, 39.9097, 50.0)
        err = max(abs(ref[i] - proj[i]) for i in range(3))
        print(f"PROJ cross-check max |Δ| = {err:.3e} m", file=sys.stderr)
        if err > 1e-7:
            raise SystemExit(f"PROJ disagreement too large: {err}")
    except ImportError:
        print("pyproj not installed; skipping PROJ cross-check", file=sys.stderr)

    earth = Geocentric.WGS84()
    sphere = Geocentric(6371000.0, 0.0)
    prolate = Geocentric(6378137.0, -1.0 / 298.257223563)
    grs80 = Geocentric(6378137.0, 1.0 / 298.257222101)

    # Geocentric forward cases: (a_tag, lat, lon, h, X, Y, Z)
    # a_tag: 0=WGS84, 1=sphere, 2=prolate, 3=GRS80
    geoc_fwd = []
    pts = [
        (0.0, 0.0, 0.0),
        (0.0, 90.0, 0.0),
        (90.0, 0.0, 0.0),
        (-90.0, 45.0, 0.0),
        (45.0, -45.0, 1000.0),
        (39.9087, 116.3975, 0.0),
        (39.9087, 116.3975, 100.0),
        (-33.8688, 151.2093, 50.0),
        (89.0, 0.0, -1000.0),
        (-89.0, 180.0, 10.0),
        (1e-8, 1e-8, 0.0),
        (48.0, 2.0, 100.0),
    ]
    for tag, earth_obj in [(0, earth), (1, sphere), (2, prolate), (3, grs80)]:
        for lat, lon, h in pts:
            X, Y, Z = earth_obj.forward(lat, lon, h)
            geoc_fwd.append((tag, lat, lon, h, X, Y, Z))

    # LocalCartesian: origins and nearby points
    origins = [
        (39.9087, 116.3975, 0.0),
        (0.0, 0.0, 0.0),
        (45.0, -120.0, 100.0),
        (-45.0, 170.0, 0.0),
        (89.5, 10.0, 0.0),
        (-89.5, -10.0, 50.0),
    ]
    deltas = [
        (0.0, 0.0, 0.0),
        (0.001, 0.001, 10.0),
        (-0.01, 0.02, -5.0),
        (0.1, -0.05, 100.0),
        (0.0, 0.5, 0.0),
        (0.5, 0.0, 1000.0),
    ]
    local_fwd = []
    for lat0, lon0, h0 in origins:
        lc = LocalCartesian(lat0, lon0, h0)
        for dlat, dlon, dh in deltas:
            lat = lat0 + dlat
            if abs(lat) > 90.0:
                continue
            lon = lon0 + dlon
            h = h0 + dh
            x, y, z = lc.forward(lat, lon, h)
            local_fwd.append((lat0, lon0, h0, lat, lon, h, x, y, z))

    # Also a few absolute targets far from origin (still useful for reverse)
    far = [
        (39.9087, 116.3975, 0.0, 31.2304, 121.4737, 0.0),
        (0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        (48.0, 2.0, 100.0, 51.5, -0.12, 50.0),
    ]
    for lat0, lon0, h0, lat, lon, h in far:
        lc = LocalCartesian(lat0, lon0, h0)
        x, y, z = lc.forward(lat, lon, h)
        local_fwd.append((lat0, lon0, h0, lat, lon, h, x, y, z))

    lines = []
    lines.append("// Generated by tools/gen_local_cartesian.py from GeographicLib Geocentric/LocalCartesian algorithm.")
    lines.append("// Cross-checked against PROJ topocentric on WGS-84. Do not edit.")
    lines.append("")
    lines.append("///|")
    lines.append("/// Geocentric Forward cases: (ellipsoid_tag, lat, lon, h, X, Y, Z).")
    lines.append("/// tag 0=WGS-84, 1=sphere a=6371000 f=0, 2=prolate f=-1/298.257223563, 3=GRS80.")
    lines.append("let geocentric_forward : Array[(Int32, Double, Double, Double, Double, Double, Double)] = [")
    for row in geoc_fwd:
        tag, lat, lon, h, X, Y, Z = row
        lines.append(
            f"  ({tag}, {lit(lat)}, {lit(lon)}, {lit(h)}, {lit(X)}, {lit(Y)}, {lit(Z)}),"
        )
    lines.append("]")
    lines.append("")
    lines.append("///|")
    lines.append("/// LocalCartesian Forward: (lat0, lon0, h0, lat, lon, h, x, y, z) on WGS-84.")
    lines.append("let local_cartesian_forward : Array[(Double, Double, Double, Double, Double, Double, Double, Double, Double)] = [")
    for row in local_fwd:
        lat0, lon0, h0, lat, lon, h, x, y, z = row
        lines.append(
            f"  ({lit(lat0)}, {lit(lon0)}, {lit(h0)}, {lit(lat)}, {lit(lon)}, {lit(h)}, {lit(x)}, {lit(y)}, {lit(z)}),"
        )
    lines.append("]")
    lines.append("")

    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT}: {len(geoc_fwd)} geocentric + {len(local_fwd)} local cases", file=sys.stderr)


if __name__ == "__main__":
    main()
