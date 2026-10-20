"""Geodesic-intersection reference values from GeographicLib's IntersectTool.

IntersectTool is the command-line front end of GeographicLib::Intersect
(Charles Karney, GeographicLib, MIT/X11 License). It ships in a
GeographicLib build; this file was checked against version 2.5
(https://github.com/geographiclib/geographiclib/tree/r2.5). The vectors
are numeric results of that tool. They are not a copy of the Intersect
source.

Modes match IntersectTool: -c closest, -o closest with an offset, -n next,
-i segment. Ellipsoids are WGS-84, a sphere, a prolate twin of WGS-84, and
GRS80. The series geodesic is used (no -E).

Usage:
  python3 tools/gen_intersect.py [path/to/IntersectTool]
"""
import subprocess
import sys

TOOL = sys.argv[1] if len(sys.argv) > 1 else "IntersectTool"
OUT = "src/geodesic/intersect_data_test.mbt"

# (code, extra IntersectTool arguments)
ELLIPSOIDS = [
    (0, []),
    (1, ["-e", "6378137", "0"]),
    (2, ["-e", "6378137", "-1/298.257223563"]),
    (3, ["-e", "6378137", "1/298.257222101"]),
]


def lit(tok):
    # Keep IntersectTool's digits. Do not round-trip through binary float.
    t = tok.strip()
    low = t.lower()
    if low == "nan":
        return "@double.not_a_number"
    if low == "inf":
        return "@double.infinity"
    if low == "-inf":
        return "@double.neg_infinity"
    if "e" in low or "E" in t:
        raise SystemExit("unexpected exponent in IntersectTool output: " + t)
    if "." not in t:
        return t + ".0"
    return t


def version():
    proc = subprocess.run([TOOL, "--version"], text=True, capture_output=True, check=False)
    text = (proc.stdout or proc.stderr).strip()
    if "GeographicLib version" not in text:
        raise SystemExit("IntersectTool --version failed:\n" + text)
    return text[text.index("GeographicLib version") :]


def solve(mode, extra, lines):
    proc = subprocess.run(
        [TOOL, mode, "-p", "10", *extra],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit("IntersectTool failed: " + mode + " " + " ".join(extra))
    rows = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("ERROR"):
            raise SystemExit("unexpected IntersectTool line: " + line)
        rows.append(parts)
    if len(rows) != len(lines):
        raise SystemExit(
            "IntersectTool returned %d rows for %d inputs" % (len(rows), len(lines))
        )
    return rows


def fmt_row(nums):
    return "[" + ", ".join(lit(v) if isinstance(v, str) else lit("%.12f" % v) for v in nums) + "]"


# Closest: latX lonX aziX latY lonY aziY -> x y c
CLOSEST = [
    (0.0, 0.0, 45.0, 45.0, 10.0, 135.0),
    (50.0, -4.0, -147.7, 0.0, 0.0, 90.0),
    (50.0, -4.0, -147.7, 0.0, 180.0, 0.0),
    (0.0, 0.0, 90.0, 0.0, 10.0, 90.0),
    (0.0, 0.0, 90.0, 0.0, 0.0, -90.0),
    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (0.0, 179.0, 90.0, 0.0, -179.0, 90.0),
    (90.0, 0.0, 0.0, 90.0, 10.0, 90.0),
    (90.0, 0.0, 0.0, -90.0, 0.0, 0.0),
    (89.9, 0.0, 90.0, 89.9, 180.0, 90.0),
    (0.0, 0.0, 10.0, 0.0, 0.0, 10.0000001),
    (0.0, 0.0, 0.0, 10.0, 0.0, 90.0),
    (39.9087, 116.3975, 152.94, 31.2304, 121.4737, 0.0),
    (-45.0, 20.0, 30.0, 45.0, -20.0, 200.0),
    (0.0, 0.0, 45.0, 0.0, 1.0, 45.0),
    (80.0, 170.0, 10.0, -80.0, -170.0, 350.0),
    (-89.0, 0.0, 180.0, 89.0, 0.0, 0.0),
    (0.0, 0.0, 90.0, 0.0, 90.0, 0.0),
    (10.0, -170.0, 45.0, -10.0, 170.0, 135.0),
    (0.0, 0.0, 1.0, 0.0, 0.001, 179.0),
    (60.0, 30.0, 270.0, -60.0, 30.0, 90.0),
    (0.0, 180.0, 0.0, 0.0, -180.0, 180.0),
    (45.0, 45.0, 45.0, -45.0, -45.0, 225.0),
    (1.0e-8, 1.0e-8, 20.0, -1.0e-8, 2.0e-8, 200.0),
    (-30.0, 10.0, 15.0, -25.0, -20.0, 85.0),
    (0.0, 10.0, 75.0, 5.0, -20.0, 145.0),
    (40.0, 10.0, 120.0, 45.0, -20.0, 190.0),
    (-30.0, 100.0, 15.0, 20.0, -100.0, 250.0),
    (15.0, -45.0, 300.0, -15.0, 45.0, 30.0),
]
CLOSEST_SHORT = CLOSEST[:8]

# Offset closest: latX lonX aziX latY lonY aziY x0 y0 -> x y c
OFFSET = [
    (50.0, -4.0, -147.7, 0.0, 180.0, 0.0, 19529110.0, -5932344.0),
    (0.0, 0.0, 45.0, 45.0, 10.0, 135.0, 0.0, 0.0),
    (0.0, 0.0, 45.0, 45.0, 10.0, 135.0, 1.0e7, 1.0e7),
    (0.0, 0.0, 90.0, 0.0, 10.0, 90.0, 1.0e6, -1.0e6),
    (0.0, 0.0, 0.0, 0.0, 90.0, 90.0, 5.0e6, -2.0e6),
    (20.0, 30.0, 10.0, -20.0, 40.0, 200.0, -8.0e6, 3.0e6),
]
OFFSET_SHORT = OFFSET[:2]

# Next: lat lon aziX aziY -> x y c
NEXT = [
    (0.0, 0.0, 45.0, 135.0),
    (0.0, 0.0, 45.0, 90.0),
    (90.0, 0.0, 0.0, 90.0),
    (0.0, 0.0, 90.0, 90.0),
    (0.0, 0.0, 90.0, -90.0),
    (-90.0, 0.0, 0.0, 45.0),
    (0.0, 179.0, 10.0, 80.0),
    (45.0, -90.0, 0.0, 180.0),
    (10.0, 10.0, 20.0, 20.5),
    (0.0, 0.0, 1.0, 2.0),
    (70.0, 140.0, 300.0, 40.0),
    (-20.0, 179.5, 90.0, 270.0),
]
NEXT_SHORT = NEXT[:4]

# Segment: latX1 lonX1 latX2 lonX2 latY1 lonY1 latY2 lonY2 -> x y c k
# The equatorial pairs are the coincident cases from GeographicLib's
# intersecttest.cpp, stated here in degrees on the real ellipsoid.
SEGMENT = [
    (0.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, 0.0),
    (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0),
    (0.0, 0.0, 0.0, 10.0, 0.0, 20.0, 0.0, 30.0),
    (0.0, 0.0, 10.0, 10.0, 0.0, 0.0, 10.0, 0.0),
    (10.0, 0.0, 20.0, 0.0, 0.0, 10.0, 0.0, 20.0),
    (0.0, 179.0, 0.0, -179.0, 1.0, 179.0, 1.0, -179.0),
    (-10.0, -10.0, 10.0, 10.0, -10.0, 10.0, 10.0, -10.0),
    (39.9, 116.3, 39.95, 116.5, 39.92, 116.2, 39.88, 116.45),
    (0.0, 0.0, 0.0, 40.0, 0.0, -20.0, 0.0, -10.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, -10.0, 0.0, 0.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, -8.0, 0.0, 2.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, -2.0, 0.0, 8.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 0.0, 0.0, 10.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 2.0, 0.0, 12.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 15.0, 0.0, 25.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 30.0, 0.0, 40.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 32.0, 0.0, 42.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 38.0, 0.0, 48.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 40.0, 0.0, 50.0),
    (0.0, 0.0, 0.0, 40.0, 0.0, 50.0, 0.0, 60.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, -20.0, 0.0, -10.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, -8.0, 0.0, 2.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, -2.0, 0.0, 8.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 2.0, 0.0, 12.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 15.0, 0.0, 25.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 30.0, 0.0, 40.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 32.0, 0.0, 42.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 38.0, 0.0, 48.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 40.0, 0.0, 50.0),
    (0.0, 40.0, 0.0, 0.0, 0.0, 50.0, 0.0, 60.0),
]
SEGMENT_SHORT = SEGMENT[:6]


def num(v):
    # Decimal form. Scientific notation is rejected: IntersectTool treats
    # the letter e as an east-longitude marker.
    return "%.12f" % float(v)


def inputs(cases):
    return [" ".join(num(v) for v in row) for row in cases]


def emit(name, rows):
    lines = ["let %s : Array[FixedArray[Double]] = [" % name]
    for row in rows:
        lines.append("  " + fmt_row(row) + ",")
    lines.append("]")
    return "\n".join(lines)


def collect(mode, cases, short, n_out):
    rows = []
    for code, extra in ELLIPSOIDS:
        use = cases if code == 0 else short
        got = solve(mode, extra, inputs(use))
        if any(len(g) != n_out for g in got):
            raise SystemExit("bad %s row width: %s" % (mode, got[0]))
        for src, out in zip(use, got):
            rows.append(["%.1f" % code, *["%.12f" % v for v in src], *out])
    return rows


def main():
    ver = version()
    closest = collect("-c", CLOSEST, CLOSEST_SHORT, 3)
    offset = collect("-o", OFFSET, OFFSET_SHORT, 3)
    nxt = collect("-n", NEXT, NEXT_SHORT, 3)
    segment = collect("-i", SEGMENT, SEGMENT_SHORT, 4)
    body = "\n\n".join(
        [
            emit("ix_closest", closest),
            emit("ix_offset", offset),
            emit("ix_next", nxt),
            emit("ix_segment", segment),
        ]
    )
    text = (
        "// Generated by tools/gen_intersect.py. Do not edit.\n"
        "// Reference: %s\n"
        "// Series geodesic (IntersectTool without -E). Lengths in metres.\n"
        "// Ellipsoid code: 0 WGS-84, 1 sphere a=6378137 f=0,\n"
        "// 2 prolate a=6378137 f=-1/298.257223563, 3 GRS80.\n"
        "// ix_closest: code latX lonX aziX latY lonY aziY x y c\n"
        "// ix_offset: code latX lonX aziX latY lonY aziY x0 y0 x y c\n"
        "// ix_next: code lat lon aziX aziY x y c\n"
        "// ix_segment: code latX1 lonX1 latX2 lonX2 latY1 lonY1 latY2 lonY2 x y c k\n"
        "\n"
        "///|\n%s\n"
    ) % (ver, body)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(text)
    print(
        "wrote %s (%d closest, %d offset, %d next, %d segment) from %s"
        % (OUT, len(closest), len(offset), len(nxt), len(segment), ver)
    )


if __name__ == "__main__":
    main()
