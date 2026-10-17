"""Rhumb-line reference values from GeographicLib's RhumbSolve (WGS-84).

RhumbSolve is an independent implementation of the rhumb inverse, direct,
and area problems (Charles Karney, GeographicLib, MIT/X11 License). It ships
in the Debian/Ubuntu package geographiclib-tools, and it can be built from
https://github.com/geographiclib/geographiclib (this file was checked
against version 2.5). The vectors cover poles, the antimeridian, and due
east/west lines.

Usage:
  python3 tools/gen_rhumb.py [path/to/RhumbSolve]
"""
import subprocess
import sys

RHUMB = sys.argv[1] if len(sys.argv) > 1 else "RhumbSolve"
OUT = "src/geodesic/rhumb_data_test.mbt"


def lit(v):
    if v != v:  # NaN
        return "@double.not_a_number"
    r = repr(float(v))
    if r == "inf":
        return "@double.infinity"
    if r == "-inf":
        return "@double.neg_infinity"
    if "e" in r:
        m, e = r.split("e")
        if "." not in m:
            m += ".0"
        return m + "e" + e
    if "." not in r:
        return r + ".0"
    return r


def solve(args, lines):
    proc = subprocess.run(
        [RHUMB, "-p", "10", *args],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit("RhumbSolve failed: " + " ".join(args))
    rows = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("ERROR"):
            raise SystemExit("unexpected RhumbSolve line: " + line)
        rows.append([float(p) for p in parts])
    if len(rows) != len(lines):
        raise SystemExit("RhumbSolve returned %d rows for %d inputs" % (len(rows), len(lines)))
    return rows


def unique(seq):
    seen = set()
    out = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def inverse_cases():
    lats = [-90.0, -89.0, -45.0, 0.0, 45.0, 89.0, 90.0]
    lons = [-180.0, -179.0, -90.0, 0.0, 90.0, 179.0, 180.0]
    pts = [(la, lo) for la in lats for lo in lons]
    cases = []
    for p in pts[::7]:
        cases.append((p[0], p[1], p[0], p[1]))
    n = len(pts)
    for i in range(n):
        for j in ((i + 1) % n, (i * 3 + 7) % n):
            if i != j:
                a = pts[i]
                b = pts[j]
                cases.append((a[0], a[1], b[0], b[1]))
    cases.extend([
        (0.0, 179.0, 0.0, -179.0),
        (0.0, -179.0, 0.0, 179.0),
        (45.0, 10.0, 45.0, 20.0),
        (45.0, 20.0, 45.0, 10.0),
        (-30.0, 0.0, -30.0, 180.0),
        (10.0, 0.0, 10.0, 180.0),
        (0.0, 0.0, 0.0, 0.000001),
        (0.0, 0.0, 0.000001, 0.0),
        (89.9, 0.0, 89.9, 0.1),
        (-89.9, 10.0, -89.9, -10.0),
        (90.0, 0.0, 0.0, 0.0),
        (90.0, 0.0, 0.0, 50.0),
        (-90.0, 20.0, 0.0, 20.0),
        (90.0, 0.0, -90.0, 0.0),
        (90.0, 0.0, 90.0, 40.0),
        (-90.0, 0.0, -90.0, 10.0),
        (90.0, 20.0, -90.0, 30.0),
        (39.9087, 116.3975, 31.2304, 121.4737),
        (51.4778, -0.0015, 35.6762, 139.6503),
        (0.0, 0.0, 0.0, 90.0),
        (0.0, 0.0, 0.0, -90.0),
        (60.0, 179.9, -60.0, -179.9),
        (-45.0, 179.0, -45.0, -179.0),
        (80.0, -180.0, 80.0, 180.0),
    ])
    return unique(cases)


def direct_cases():
    starts = [
        (0.0, 0.0),
        (0.0, 179.0),
        (45.0, -170.0),
        (-45.0, 10.0),
        (89.0, 0.0),
        (90.0, 0.0),
        (-90.0, 30.0),
        (39.9, 116.4),
    ]
    azis = [0.0, 90.0, -90.0, 180.0, 45.0, 89.9]
    dists = [0.0, 10.0, 100000.0, 2000000.0, 15000000.0, -1000.0]
    cases = []
    for lat, lon in starts:
        for azi in azis:
            for s in dists:
                cases.append((lat, lon, azi, s))
    cases.extend([
        (80.0, 10.0, 0.0, 3000000.0),
        (0.0, 179.0, 90.0, 500000.0),
        (40.6, -73.8, 51.5, 5500000.0),
        (10.0, 20.0, 45.0, -1000.0),
        (0.0, 0.0, 0.0, 20000000.0),
        (-80.0, -20.0, 180.0, 4000000.0),
        (0.0, 179.5, 90.0, 200000.0),
        (0.0, -179.5, -90.0, 200000.0),
    ])
    return unique(cases)


def emit(name, rows):
    lines = ["///|", "let %s : Array[FixedArray[Double]] = [" % name]
    for row in rows:
        lines.append("  [" + ", ".join(lit(v) for v in row) + "],")
    lines.append("]")
    lines.append("")
    return lines


def main():
    inv = inverse_cases()
    dire = direct_cases()
    # Fixed point: RhumbSolve treats the letter e as an east-hemisphere mark,
    # so scientific notation is not safe on its input.
    def line(c):
        return " ".join("%.16f" % v for v in c)

    inv_out = solve(["-i"], [line(c) for c in inv])
    dir_out = solve([], [line(c) for c in dire])
    un_out = solve(["-u"], [line(c) for c in dire])
    inv_rows = [list(c) + inv_out[i] for i, c in enumerate(inv)]
    dir_rows = [list(c) + dir_out[i] for i, c in enumerate(dire)]
    un_rows = [list(c) + un_out[i] for i, c in enumerate(dire)]
    # inverse: azi s12 S12; direct: lat2 lon2 S12
    if any(len(r) != 7 for r in inv_rows + dir_rows + un_rows):
        raise SystemExit("unexpected column count")
    ver = subprocess.run([RHUMB, "--version"], text=True, capture_output=True, check=True)
    # "RhumbSolve: GeographicLib version 2.5" — drop the local executable path.
    ver_txt = ver.stdout.strip().split(": ", 1)[-1]
    header = [
        "// Generated by tools/gen_rhumb.py. Do not edit.",
        "// Reference: GeographicLib RhumbSolve (%s), WGS-84, series (not -E)." % ver_txt,
        "// Inverse columns: lat1 lon1 lat2 lon2 azi12 s12 area12.",
        "// Direct columns: lat1 lon1 azi12 s12 lat2 lon2 area12.",
        "// rhumb_direct_unroll is the direct problem with longitude unrolled.",
        "",
    ]
    body = []
    body.extend(emit("rhumb_inverse", inv_rows))
    body.extend(emit("rhumb_direct", dir_rows))
    body.extend(emit("rhumb_direct_unroll", un_rows))
    open(OUT, "w").write("\n".join(header + body))
    print(
        "wrote %s (%d inverse, %d direct)" % (OUT, len(inv_rows), len(dir_rows)),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
