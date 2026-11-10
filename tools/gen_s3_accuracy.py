"""Generate src/examples/s3_accuracy/data.mbt from GeodTest-short.dat.

Takes every 50th line (200 cases) of Karney's GeodTest-short.dat (CC0,
https://doi.org/10.5281/zenodo.32156; download GeodTest-short.dat.gz from
https://sourceforge.net/projects/geographiclib/files/testdata/) and keeps
the columns lat1 lon1 azi1 lat2 lon2 azi2 s12. The output is formatted
with `moon fmt` afterwards.

Usage: python3 tools/gen_s3_accuracy.py /path/to/GeodTest-short.dat
"""
import hashlib
import sys

SHA256 = "601de5a47f798a3c4067ca9326e8009c0acd60cd2d355b7ba3080bf114734710"
OUT = "src/examples/s3_accuracy/data.mbt"


def lit(x):
    r = repr(float(x))
    if "e" in r:
        m, e = r.split("e")
        if "." not in m:
            m += ".0"
        r = m + "e" + e
    elif "." not in r:
        r += ".0"
    return r


def main():
    raw = open(sys.argv[1], "rb").read()
    if hashlib.sha256(raw).hexdigest() != SHA256:
        sys.exit("GeodTest-short.dat: unexpected sha256")
    lines = raw.decode("ascii").split("\n")
    rows = [l.split()[:7] for i, l in enumerate(lines) if l.strip() and i % 50 == 0]
    out = [
        "///|",
        "/// 200 cases (every 50th line) of GeographicLib's GeodTest-short.dat",
        "/// (C. Karney, CC0, https://doi.org/10.5281/zenodo.32156).",
        "/// Columns: lat1 lon1 azi1 lat2 lon2 azi2 s12.",
        "let cases : Array[FixedArray[Double]] = [",
    ]
    for r in rows:
        out.append("  [" + ", ".join(lit(x) for x in r) + "],")
    out.append("]")
    open(OUT, "w").write("\n".join(out) + "\n")
    print(len(rows), "rows")


if __name__ == "__main__":
    main()
