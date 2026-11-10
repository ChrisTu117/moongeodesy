"""Generate src/cn/boundary_data.mbt: where GCJ-02 is NOT applied inside the box.

GCJ-02 is applied inside the usual box (lon 72.004..137.8347, lat
0.8293..55.8271) except on foreign land and in the sea next to foreign
coasts. This script builds that "foreign" region from Natural Earth (public
domain, China point of view):

1. Natural Earth 1:10m admin-0 countries, China point of view, v5.1.2: every
   country except China, Hong Kong, Macau and Scarborough Shoal.
2. Natural Earth 1:10m minor islands, v5.1.2: an islet is added when it is
   closer to a foreign country than to China and within 0.5 degree of that
   country (this adds islands such as Ly Son, Vietnam, that are too small for
   the countries layer; disputed features of the China point of view are
   already inside China's polygons and are never added).
3. The union is grown by BUFFER degrees (about 9 km), so coastal towns,
   harbours, and small islands are inside it even after simplification, and
   then Chinese land is cut out exactly. Near China's land borders the
   boundary is therefore China's border.
4. The region is simplified (TOL degrees, about 1 km), clipped to 1-degree
   cells, and written with a per-cell index: each cell is empty, full, or a
   list of rings tested with the even-odd rule.

A self-check runs before writing: coastal and island places abroad must be
foreign; Chinese cities, including coastal and island ones, must not be.

Needs: pip install shapely
Usage:
    python3 tools/gen_china_boundary.py [countries.geojson [minor_islands.geojson]]
Without paths the files are downloaded from Natural Earth's GitHub mirror.
Their SHA-256 is checked.
"""

import hashlib
import json
import math
import os
import sys
import urllib.request

from shapely.geometry import Point, box, shape
from shapely.ops import unary_union
from shapely.strtree import STRtree

BASE = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/v5.1.2/geojson/"
SOURCES = {
    "countries": (
        BASE + "ne_10m_admin_0_countries_chn.geojson",
        "a13bf5f310fde87bc0a5f994f8ce9bd706cc198d8ee37d221e61c2546b945372",
    ),
    "minor": (
        BASE + "ne_10m_minor_islands.geojson",
        "8c933ca7a4760256bdc46408355706e39764b0fa01c160c28888b90b4faec29f",
    ),
}
INSIDE = {"CHN", "HKG", "MAC", "SCR"}  # China, including Hong Kong and Macau
W, S, E, N = 72.004, 0.8293, 137.8347, 55.8271
BUFFER = 0.08  # degrees of sea added around foreign land (about 9 km)
TOL = 0.01  # simplification tolerance, degrees (about 1.1 km)
MINOR_REACH = 0.5  # an islet must be this close (degrees) to a foreign country
CELL_LON0, CELL_LAT0, NCOLS, NROWS = 72, 0, 66, 56  # 1-degree cells
OUT = os.path.join(os.path.dirname(__file__), "..", "src", "cn", "boundary_data.mbt")

# Self-check. Foreign: capitals, and coastal and island places that the
# previous (simplified-coast) data got wrong. China: inland, coastal and
# island cities, Hong Kong, Macau, Taiwan, and the South China Sea.
FOREIGN = {
    "Seoul": (37.5665, 126.9780), "Busan": (35.1796, 129.0756),
    "Pyongyang": (39.0392, 125.7625), "Nampo": (38.7375, 125.4078),
    "Osaka": (34.6937, 135.5023), "Naha airport": (26.1958, 127.6459),
    "Ishigaki": (24.3448, 124.1572), "Yonaguni": (24.4672, 123.0045),
    "Ulleungdo": (37.4845, 130.9057), "Hanoi": (21.0285, 105.8542),
    "Ha Long": (20.9599, 107.0425), "Vung Tau": (10.3460, 107.0843),
    "Con Dao": (8.6826, 106.6090), "Ly Son": (15.3800, 109.1190),
    "Ho Chi Minh City": (10.8231, 106.6297), "Vientiane": (17.9757, 102.6331),
    "Bangkok": (13.7563, 100.5018), "Pattaya": (12.9236, 100.8825),
    "Phnom Penh": (11.5564, 104.9282), "Yangon": (16.8409, 96.1735),
    "Sittwe": (20.1460, 92.8983), "Kawthaung": (9.9822, 98.5503),
    "Kuala Lumpur": (3.1390, 101.6869), "Penang": (5.4141, 100.3288),
    "Malacca": (2.1896, 102.2501), "Kota Kinabalu": (5.9804, 116.0735),
    "Singapore": (1.3521, 103.8198), "Manila": (14.5995, 120.9842),
    "Iloilo": (10.7202, 122.5621), "Puerto Princesa": (9.7392, 118.7353),
    "Boracay": (11.9674, 121.9248), "Basco": (20.4487, 121.9702),
    "Koror": (7.3419, 134.4792), "Kochi": (9.9312, 76.2673),
    "Male": (4.1755, 73.5093), "Colombo": (6.9271, 79.8612),
    "Ulaanbaatar": (47.8864, 106.9057), "Vladivostok": (43.1198, 131.8869),
    "Almaty": (43.2220, 76.8512), "Bishkek": (42.8746, 74.5698),
    "Kathmandu": (27.7172, 85.3240), "Thimphu": (27.4728, 89.6390),
    "Dhaka": (23.8103, 90.4125), "New Delhi": (28.6139, 77.2090),
    "Islamabad": (33.6844, 73.0479),
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
    "Yongxing Island": (16.8340, 112.3380), "Dalian": (38.9140, 121.6147),
    "Qingdao": (36.0671, 120.3826), "Xiamen": (24.4798, 118.0894),
    "Zhoushan": (29.9853, 122.2072), "Beihai": (21.4733, 109.1201),
    "Weihai": (37.5131, 122.1204), "Hunchun": (42.8679, 130.3655),
}


def load(name, path):
    url, digest = SOURCES[name]
    if path is None:
        with urllib.request.urlopen(url) as r:
            data = r.read()
    else:
        with open(path, "rb") as fh:
            data = fh.read()
    got = hashlib.sha256(data).hexdigest()
    if got != digest:
        sys.exit(f"{name}: SHA-256 mismatch: {got}")
    return json.loads(data)


def parts(geom):
    return list(getattr(geom, "geoms", [geom]))


def build(countries, minor):
    """The foreign region inside the box, before simplification."""
    frame = box(W - 1.0, S - 1.0, E + 1.0, N + 1.0)
    china = unary_union(
        [shape(f["geometry"]) for f in countries["features"]
         if f["properties"]["ADM0_A3"] in INSIDE]
    ).intersection(frame)
    foreign = unary_union(
        [shape(f["geometry"]).intersection(frame) for f in countries["features"]
         if f["properties"]["ADM0_A3"] not in INSIDE]
    )
    fparts = parts(foreign)
    tree = STRtree(fparts)
    islets = []
    for f in minor["features"]:
        m = shape(f["geometry"])
        if not frame.intersects(m):
            continue
        df = fparts[tree.nearest(m)].distance(m)
        if df < MINOR_REACH and df < china.distance(m):
            islets.append(m)
    region = unary_union([foreign] + islets).buffer(BUFFER, 8).difference(china)
    return region.intersection(box(W, S, E, N)), china, len(islets)


def rings_of(geom):
    out = []
    for p in parts(geom):
        if p.is_empty or p.geom_type != "Polygon":
            continue
        out.append([(round(x, 4), round(y, 4)) for x, y in p.exterior.coords[:-1]])
        for hole in p.interiors:
            out.append([(round(x, 4), round(y, 4)) for x, y in hole.coords[:-1]])
    return [r for r in out if len(r) >= 3]


def even_odd(rings, lon, lat):
    inside = False
    for pts in rings:
        for i in range(len(pts)):
            (x1, y1), (x2, y2) = pts[i - 1], pts[i]
            if (y1 > lat) != (y2 > lat) and lon < x1 + (lat - y1) * (x2 - x1) / (y2 - y1):
                inside = not inside
    return inside


def main():
    countries = load("countries", sys.argv[1] if len(sys.argv) > 1 else None)
    minor = load("minor", sys.argv[2] if len(sys.argv) > 2 else None)
    region, _, n_islets = build(countries, minor)
    region = region.simplify(TOL, preserve_topology=True)

    kinds, cell_start, ring_start, xy = [], [0], [0], []
    cells = {}
    for j in range(NROWS):
        for i in range(NCOLS):
            cell = box(CELL_LON0 + i, CELL_LAT0 + j, CELL_LON0 + i + 1, CELL_LAT0 + j + 1)
            piece = region.intersection(cell)
            if piece.is_empty or piece.area == 0.0:
                kinds.append(0)
            elif abs(piece.area - cell.area) < 1e-12 and len(rings_of(piece)) == 1:
                kinds.append(1)
            else:
                kinds.append(2)
                rings = rings_of(piece)
                cells[(i, j)] = rings
                for pts in rings:
                    for x, y in pts:
                        xy += [x, y]
                    ring_start.append(len(xy) // 2)
            cell_start.append(len(ring_start) - 1)

    def foreign(lat, lon):
        i, j = math.floor(lon - CELL_LON0), math.floor(lat - CELL_LAT0)
        k = kinds[j * NCOLS + i]
        return k == 1 or (k == 2 and even_odd(cells[(i, j)], lon, lat))

    bad = [n for n, (la, lo) in FOREIGN.items() if not foreign(la, lo)]
    bad += [n for n, (la, lo) in CHINA.items() if foreign(la, lo)]
    if bad:
        sys.exit(f"self-check failed for: {bad}")

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

    nfull, npart = kinds.count(1), kinds.count(2)
    out = [
        "// Generated by tools/gen_china_boundary.py. Do not edit.",
        "// Source: Natural Earth 1:10m admin-0 countries, China point of view, and",
        "// 1:10m minor islands (v5.1.2, public domain), https://www.naturalearthdata.com/",
        f"// Foreign land and islets, grown by {BUFFER} deg of sea, minus Chinese land,",
        f"// simplified to {TOL} deg; {n_islets} minor islands added.",
        f"// 1-degree cells from ({CELL_LON0}E, {CELL_LAT0}N): {nfull} full, {npart} partial;",
        f"// {len(ring_start) - 1} rings, {len(xy) // 2} vertices.",
        "",
        "///|",
        "/// Cell grid: first longitude, first latitude, columns, rows (1-degree cells).",
        f"let cell_lon0 : Double = {CELL_LON0}.0",
        "",
        "///|",
        f"let cell_lat0 : Double = {CELL_LAT0}.0",
        "",
        "///|",
        f"let cell_cols : Int = {NCOLS}",
        "",
        "///|",
        f"let cell_rows : Int = {NROWS}",
        "",
        "///|",
        "/// Per cell, row by row: 0 = no foreign region, 1 = all foreign, 2 = test rings.",
        "let cell_kind : FixedArray[Int] = [",
        arr(kinds, str),
        "]",
        "",
        "///|",
        "/// First ring of each cell in `ring_start`; one extra entry at the end.",
        "let cell_ring_start : FixedArray[Int] = [",
        arr(cell_start, str),
        "]",
        "",
        "///|",
        "/// First vertex of each ring in `foreign_xy`; one extra entry at the end.",
        "let ring_start : FixedArray[Int] = [",
        arr(ring_start, str),
        "]",
        "",
        "///|",
        "/// Ring vertices as lon, lat pairs (degrees).",
        "let foreign_xy : FixedArray[Double] = [",
        arr(xy, lambda v: f"{v:.4f}"),
        "]",
        "",
    ]
    with open(OUT, "w") as fh:
        fh.write("\n".join(out))
    print(f"wrote {OUT}: {nfull} full cells, {npart} partial, {len(ring_start) - 1} rings, "
          f"{len(xy) // 2} vertices, {n_islets} islets")


if __name__ == "__main__":
    main()
