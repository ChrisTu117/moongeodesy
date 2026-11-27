"""Library line coverage from `moon coverage analyze -- -f summary`.

Reads the summary on stdin, adds up the lines of the library packages
(everything under src/ except src/examples/, whose programs run through
`moon run`, not tests), prints the result, and exits with status 1 when it
is below the minimum percentage given as the first argument (default 90).

Usage: moon coverage analyze -- -f summary | python3 tools/coverage.py 90
"""
import re
import sys

minimum = float(sys.argv[1]) if len(sys.argv) > 1 else 90.0
covered = total = 0
for line in sys.stdin:
    m = re.match(r"(src/\S+\.mbt): (\d+)/(\d+)\s*$", line.strip())
    if m and not m.group(1).startswith("src/examples/"):
        covered += int(m.group(2))
        total += int(m.group(3))
if total == 0:
    sys.exit("no coverage summary on stdin")
pct = 100.0 * covered / total
print(f"library line coverage: {covered}/{total} = {pct:.1f}% (minimum {minimum:g}%)")
sys.exit(0 if pct >= minimum else 1)
