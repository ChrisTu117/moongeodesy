"""Generate src/cn/boundary_data.mbt: foreign land inside the GCJ-02 box.

GCJ-02 is applied inside the usual box (lon 72.004..137.8347, lat 0.8293..55.8271)
except on the land of other countries. This script takes Natural Earth's
1:10m admin-0 countries, China point of view (public domain), keeps every
country except China, Hong Kong, Macau and Scarborough Shoal, clips each outer
ring to the box, simplifies it, and writes the rings as MoonBit arrays.

Usage (stdlib only):
    python3 tools/gen_china_boundary.py [path/to/ne_10m_admin_0_countries_chn.geojson]
Without a path the file is downloaded from SRC_URL. Its SHA-256 is checked.
"""

import hashlib
import json
import os
import sys
import urllib.request

SRC_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "v5.1.2/geojson/ne_10m_admin_0_countries_chn.geojson"
)
SRC_SHA256 = "a13bf5f310fde87bc0a5f994f8ce9bd706cc198d8ee37d221e61c2546b945372"
INSIDE = {"CHN", "HKG", "MAC", "SCR"}  # treated as inside China (no exclusion)
W, S, E, N = 72.004, 0.8293, 137.8347, 55.8271
TOL = 0.01  # Douglas-Peucker tolerance near China, degrees (about 1.1 km)
TOL_FAR = 0.05  # tolerance more than about NEAR_CELLS cells from Chinese land
CELL = 0.25  # grid cell (degrees) used to decide "near China"
NEAR_CELLS = 2
MIN_AREA = 0.002  # drop clipped rings smaller than this, square degrees
OUT = os.path.join(os.path.dirname(__file__), "..", "src", "cn", "boundary_data.mbt")

# Points well away (>= 10 km) from any land border, used as a self-check.
FOREIGN = {
    "Seoul": (37.5665, 126.9780), "Busan": (35.1796, 129.0756),
    "Pyongyang": (39.0392, 125.7625), "Osaka": (34.6937, 135.5023),
    "Hanoi": (21.0285, 105.8542), "Ho Chi Minh City": (10.8231, 106.6297),
    "Vientiane": (17.9757, 102.6331), "Bangkok": (13.7563, 100.5018),
    "Phnom Penh": (11.5564, 104.9282), "Yangon": (16.8409, 96.1735),
    "Kuala Lumpur": (3.1390, 101.6869), "Singapore": (1.3521, 103.8198),
    "Manila": (14.5995, 120.9842), "Ulaanbaatar": (47.8864, 106.9057),
    "Vladivostok": (43.1198, 131.8869), "Almaty": (43.2220, 76.8512),
    "Bishkek": (42.8746, 74.5698), "Kathmandu": (27.7172, 85.3240),
    "Thimphu": (27.4728, 89.6390), "Dhaka": (23.8103, 90.4125),
    "New Delhi": (28.6139, 77.2090), "Islamabad": (33.6844, 73.0479),
}
CHINA = {
    "Beijing": (39.9087, 116.3975), "Shanghai": (31.2304, 121.4737),
    "Shenzhen": (22.5431, 114.0579), "Kashgar": (39.4704, 75.9898),
    "Urumqi": (43.8256, 87.6168), "Lhasa": (29.6500, 91.1000),
    "Kunming": (25.0389, 102.7183), "Nanning": (22.8170, 108.3665),
    "Harbin": (45.8038, 126.5349), "Hohhot": (40.8424, 111.7490),
    "Haikou": (20.0440, 110.1983), "Sanya": (18.2528, 109.5119),
    "Hong Kong": (22.3193, 114.1694), "Macau": (22.1987, 113.5439),
    "Taipei": (25.0330, 121.5654), "Kaohsiung": (22.6273, 120.3014),
    "Yongxing Island": (16.8340, 112.3380),
}


def load(path):
    if path is None:
        with urllib.request.urlopen(SRC_URL) as r:
            data = r.read()
    else:
        with open(path, "rb") as fh:
            data = fh.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != SRC_SHA256:
        sys.exit(f"SHA-256 mismatch: {digest}")
    return json.loads(data)


def clip(ring):
    """Sutherland-Hodgman clip of a closed ring against the box."""
    def edge(pts, inside, cross):
        out = []
        for i in range(len(pts)):
            cur, prev = pts[i], pts[i - 1]
            if inside(cur):
                if not inside(prev):
                    out.append(cross(prev, cur))
                out.append(cur)
            elif inside(prev):
                out.append(cross(prev, cur))
        return out

    def at_x(x):
        return lambda a, b: (x, a[1] + (b[1] - a[1]) * (x - a[0]) / (b[0] - a[0]))

    def at_y(y):
        return lambda a, b: (a[0] + (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]), y)

    pts = [tuple(p[:2]) for p in ring[:-1]]
    for inside, cross in (
        (lambda p: p[0] >= W, at_x(W)),
        (lambda p: p[0] <= E, at_x(E)),
        (lambda p: p[1] >= S, at_y(S)),
        (lambda p: p[1] <= N, at_y(N)),
    ):
        if not pts:
            break
        pts = edge(pts, inside, cross)
    return pts


def seg_dist(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    if dx == 0 and dy == 0:
        return ((p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2) ** 0.5
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)))
    return ((p[0] - a[0] - t * dx) ** 2 + (p[1] - a[1] - t * dy) ** 2) ** 0.5


def simplify(pts, tol):
    """Iterative Douglas-Peucker on an open polyline (keeps both ends).

    `tol(a, b, c)` returns the tolerance for a segment from a to b whose
    farthest point is c, so detail is kept only near Chinese land."""
    keep = [False] * len(pts)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        best, k = 0.0, -1
        for m in range(i + 1, j):
            d = seg_dist(pts[m], pts[i], pts[j])
            if d > best:
                best, k = d, m
        if k >= 0 and best > tol(pts[i], pts[j], pts[k]):
            keep[k] = True
            stack += [(i, k), (k, j)]
    return [p for p, kp in zip(pts, keep) if kp]


def simplify_ring(pts, tol):
    # Split the ring at its two mutually far points so both halves are open.
    far = max(range(len(pts)), key=lambda m: (pts[m][0] - pts[0][0]) ** 2 + (pts[m][1] - pts[0][1]) ** 2)
    a = simplify(pts[: far + 1], tol)
    b = simplify(pts[far:] + [pts[0]], tol)
    return a[:-1] + b[:-1]


def area(pts):
    return abs(sum(pts[i - 1][0] * p[1] - p[0] * pts[i - 1][1] for i, p in enumerate(pts))) / 2.0


def contains(pts, lon, lat):
    inside = False
    for i in range(len(pts)):
        (x1, y1), (x2, y2) = pts[i - 1], pts[i]
        if (y1 > lat) != (y2 > lat) and lon < x1 + (lat - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def near_china_cells(features):
    cells = set()
    for f in features:
        if f["properties"]["ADM0_A3"] not in INSIDE:
            continue
        geom = f["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            for lon, lat in (q[:2] for q in poly[0]):
                cx, cy = int(lon // CELL), int(lat // CELL)
                for dx in range(-NEAR_CELLS, NEAR_CELLS + 1):
                    for dy in range(-NEAR_CELLS, NEAR_CELLS + 1):
                        cells.add((cx + dx, cy + dy))
    return cells


def main():
    g = load(sys.argv[1] if len(sys.argv) > 1 else None)
    cells = near_china_cells(g["features"])

    def near(p):
        return (int(p[0] // CELL), int(p[1] // CELL)) in cells

    def tol(a, b, c):
        return TOL if near(a) or near(b) or near(c) else TOL_FAR

    rings, holes, raw = [], 0, 0
    for f in g["features"]:
        code = f["properties"]["ADM0_A3"]
        if code in INSIDE:
            continue
        geom = f["geometry"]
        polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polys:
            holes += len(poly) - 1
            pts = clip(poly[0])
            if len(pts) < 3:
                continue
            raw += len(pts)
            pts = simplify_ring(pts, tol)
            if len(pts) >= 3 and area(pts) >= MIN_AREA:
                rings.append((code, pts))

    def foreign(lat, lon):
        return any(contains(p, lon, lat) for _, p in rings)

    bad = [n for n, (la, lo) in FOREIGN.items() if not foreign(la, lo)]
    bad += [n for n, (la, lo) in CHINA.items() if foreign(la, lo)]
    if bad:
        sys.exit(f"self-check failed for: {bad}")

    starts, xy, boxes = [0], [], []
    for _, pts in rings:
        for lon, lat in pts:
            xy += [lon, lat]
        starts.append(len(xy) // 2)
        lons, lats = [p[0] for p in pts], [p[1] for p in pts]
        boxes += [min(lons), min(lats), max(lons), max(lats)]

    def arr(vals, fmt):
        lines, line = [], "  "
        for v in vals:
            item = fmt(v) + ", "
            if len(line) + len(item) > 80:
                lines.append(line.rstrip())
                line = "  "
            line += item
        lines.append(line.rstrip())
        return "\n".join(lines)

    def num(v):
        return f"{v:.3f}"

    countries = sorted({c for c, _ in rings})
    out = [
        "// Generated by tools/gen_china_boundary.py. Do not edit.",
        "// Source: Natural Earth 1:10m admin-0 countries, China point of view",
        "// (v5.1.2, public domain), https://www.naturalearthdata.com/",
        f"// Outer rings of {len(countries)} countries clipped to the GCJ-02 box,",
        f"// Douglas-Peucker {TOL} deg within ~{NEAR_CELLS * CELL} deg of Chinese land, {TOL_FAR} deg",
        f"// elsewhere; rings under {MIN_AREA} sq deg dropped.",
        f"// {len(rings)} rings, {len(xy) // 2} vertices; holes ignored: {holes}.",
        f"// Countries: {' '.join(countries)}",
        "",
        "///|",
        "/// Start index (in vertices) of each ring in `foreign_xy`; one extra entry at",
        "/// the end marks the total vertex count.",
        "let foreign_ring_start : FixedArray[Int] = [",
        arr(starts, str),
        "]",
        "",
        "///|",
        "/// Per ring: min lon, min lat, max lon, max lat (degrees).",
        "let foreign_ring_bbox : FixedArray[Double] = [",
        arr(boxes, num),
        "]",
        "",
        "///|",
        "/// Ring vertices as lon, lat pairs (degrees).",
        "let foreign_xy : FixedArray[Double] = [",
        arr(xy, num),
        "]",
        "",
    ]
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    print(f"rings={len(rings)} vertices={len(xy) // 2} (clipped {raw}) holes_ignored={holes}")


if __name__ == "__main__":
    main()
