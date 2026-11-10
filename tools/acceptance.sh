#!/usr/bin/env bash
# Runs every automated check in docs/十一月验收清单.md and prints PASS/FAIL
# per checklist ID. Needs the MoonBit toolchain (`moon`) and python3; the
# public-repo check needs network access.
#
#   bash tools/acceptance.sh            # everything
#   bash tools/acceptance.sh --quick    # skip the per-backend example runs
#
# Set MOON_TARGET_DIR to build somewhere other than ./_build.
set -u
cd "$(dirname "$0")/.."

QUICK=0
[ "${1:-}" = "--quick" ] && QUICK=1
TD=()
[ -n "${MOON_TARGET_DIR:-}" ] && TD=(--target-dir "$MOON_TARGET_DIR")
LOG=$(mktemp)
FAILS=0

run() {
  local id="$1" desc="$2"
  shift 2
  if "$@" >"$LOG" 2>&1; then
    printf 'PASS  %-4s %s\n' "$id" "$desc"
  else
    printf 'FAIL  %-4s %s\n' "$id" "$desc"
    tail -n 20 "$LOG" | sed 's/^/        /'
    FAILS=$((FAILS + 1))
  fi
}

moonbit_share() {
  local mbt py
  mbt=$(find src -name '*.mbt' -print0 | xargs -0 cat | wc -l)
  py=$(find tools -name '*.py' -print0 | xargs -0 cat | wc -l)
  echo "MoonBit $mbt lines, Python tools $py lines"
  [ "$mbt" -gt $((9 * py)) ]
}

public_repo() {
  [ "$(curl -s -o /dev/null -w '%{http_code}' https://github.com/ChrisTu117/moongeodesy)" = 200 ]
}

november_history() {
  local n
  n=$(git log --since=2026-11-01 --oneline | wc -l)
  echo "$n commits since 2026-11-01"
  [ "$n" -ge 20 ]
}

license_ok() {
  grep -q 'MIT License' LICENSE && grep -q 'license = "MIT"' moon.mod
}

port_sources() {
  grep -q '## Port source and license' README.md &&
    grep -q 'GeographicLib' src/conic/lcc.mbt &&
    grep -q 'GeographicLib' src/conic/albers.mbt &&
    grep -q 'Chulliat' src/magnetic/magnetic.mbt &&
    grep -q 'GeographicLib' src/utm/tm.mbt &&
    grep -q 'IOGP' src/conic/gigs_data_test.mbt &&
    grep -q '17 U.S.C. 403' src/magnetic/wmm2025_data.mbt
}

ai_docs() {
  [ -s docs/十一月新增功能说明.md ] && grep -q '## AI assistance' README.md
}

readme_sync() {
  python3 tools/sync_readme_tests.py && moon fmt "${TD[@]}" &&
    git diff --exit-code -- src/examples/readme_check
}

mbti_ok() {
  moon info "${TD[@]}" && git diff --exit-code -- '*.mbti'
}

no_ffi() {
  ! grep -rn 'extern ' src/geodesic src/utm src/cn src/internal src/magnetic src/conic
}

examples_all_backends() {
  local t s
  for t in wasm wasm-gc js native; do
    for s in demo s1_cn_distance s2_utm_mgrs s3_accuracy s4_polygon_area s5_compass s6_conic_area; do
      echo "== $t $s"
      moon run "${TD[@]}" --target "$t" "src/examples/$s" || return 1
    done
  done
}

cli_commands() {
  local c=(moon run "${TD[@]}" src/examples/cli --)
  "${c[@]}" geod inverse 39.9087 116.3975 31.2304 121.4737 &&
    "${c[@]}" geod direct 22.5431 114.0579 90 1000000 &&
    "${c[@]}" rhumb inverse 51.4778 -0.0015 35.6762 139.6503 &&
    "${c[@]}" intersect closest 0 0 45 45 10 135 &&
    "${c[@]}" enu forward 39.9087 116.3975 39.9097 116.3985 50 &&
    "${c[@]}" coord wgs2gcj 39.9087 116.3975 &&
    "${c[@]}" coord gcj2wgs 39.91010 116.40374 &&
    "${c[@]}" coord wgs2bd 39.9087 116.3975 &&
    "${c[@]}" coord bd2wgs 39.91645 116.41012 &&
    "${c[@]}" mag 39.9087 116.3975 2026-11-10 &&
    "${c[@]}" albers forward 25 47 0 105 39.9087 116.3975 &&
    "${c[@]}" lcc reverse 25 47 0 105 1000000 4000000 &&
    moon run "${TD[@]}" src/examples/cli | grep -q 'enu reverse' &&
    ! "${c[@]}" nope
}

echo "== Official acceptance rules (A)"
run A1 "MoonBit is the main language" moonbit_share
run A2 "public repository" public_repo
run A3 "continuous commits this period" november_history
run A8 "open-source licence (MIT)" license_ok
run A9 "porting sources and licences stated" port_sources
run A10 "AI work explained (Chinese notes, README section)" ai_docs

echo "== Build and test (C, same as CI)"
run C1 "moon check --deny-warn" moon check --deny-warn "${TD[@]}"
run C2 "moon fmt --check" moon fmt --check "${TD[@]}"
run C3 "README code blocks are the tested ones" readme_sync
run C4 "interface files up to date" mbti_ok
run C5 "moon test --target all (every test, 4 backends)" moon test --target all "${TD[@]}"
run C6 "library packages use no FFI" no_ffi
if [ "$QUICK" = 0 ]; then
  run C7 "examples on every backend (each aborts on a wrong result)" examples_all_backends
fi
run C8 "CLI commands from the README" cli_commands
run C9 "moon package --list" moon package --list "${TD[@]}"

rm -f "$LOG"
echo
if [ "$FAILS" = 0 ]; then
  echo "All automated checks passed. Manual items: see docs/十一月验收清单.md (M1-M4)."
else
  echo "$FAILS check(s) failed."
  exit 1
fi
