"""Reference values used in src/examples/s4_polygon_area/main.mbt.

Prints the area and perimeter of the continental triangle and the area of
the Beijing quad, computed with Python geographiclib (Karney, MIT/X11;
`pip install geographiclib`, checked with 2.1). Paste the printed values
into the example if they ever change.

Usage: python3 tools/gen_s4_reference.py
"""
from geographiclib.geodesic import Geodesic


def area(points):
    p = Geodesic.WGS84.Polygon()
    for lat, lon in points:
        p.AddPoint(lat, lon)
    _, perim, a = p.Compute()
    return a, perim


big = [(53.56, 123.26), (36.0, 73.5), (18.16, 109.5)]
quad = [(39.9990, 116.3800), (39.9990, 116.4030), (40.0170, 116.4030), (40.0170, 116.3800)]
a, perim = area(big)
q, _ = area(quad)
print(f"ref_area = {a!r}")
print(f"ref_perim = {perim!r}")
print(f"q_ref = {q!r}")
