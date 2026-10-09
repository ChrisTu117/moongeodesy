# moongeodesy progress

Updated: 2026-10-09 (Asia/Shanghai). Toolchain: moon 0.1.20260920 / moonc v0.10.14.

## Status: all three priorities implemented and tested
`moon test --target all` → 28 tests pass on wasm, wasm-gc, js and native;
`moon check --deny-warn`, `moon fmt --check` and `moon package --list` are clean.

### 1. `cn` — WGS-84 / GCJ-02 / BD-09 ✅
- forward transforms, exact iterative inverses (round trip < 1e-9°), single-step approximate inverses, `out_of_china`
- 200 reference points from `coord-convert` (Krasovsky a=6378245): agreement < 1e-11°

### 2. `geodesic` — Karney / GeographicLib port ✅
- `inverse`, `direct`, `arc_direct`, `distance`; `GeodesicLine` (`position`, `arc_position`, `densify`, `inverse_line`, `direct_line`)
- `PolygonArea` (`add_point`, `add_edge`, `compute`, `test_point`, `test_edge`, polyline mode) + `polygon_area`
- ellipsoids: `wgs84()`, `grs80()`, `cgcs2000()`, `new(a, f)`
- GeodTest-short.dat (CC0), 1000 cases: inverse max Δs12 7.5 nm; direct max position error 8.6 nm
- vs Python geographiclib 2.1 (2020 cases): relative diff < 1e-11 (an accurate `hypot` was needed to match on very short lines)
- polygons (pole-enclosing, antimeridian) vs geographiclib: < 1e-14 × Earth area

### 3. `utm` — TM / UTM / UPS / MGRS ✅
- transverse Mercator (Krüger n^6), polar stereographic, UTM/UPS with Norway/Svalbard zones, convergence + scale
- vs PROJ 9 (pyproj): 600 points, max 5.5 nm
- MGRS encode/decode (UTM + polar UPS grids, precision 0–11): 315/315 strings identical to NGA GeoTrans

### Repo hygiene
- README.md (English + 中文), MIT LICENSE, `src/examples/demo` (run: `moon run src/examples/demo`), README snippets compiled as tests
- `.github/workflows/ci.yml`: check (deny-warn), fmt, mbti diff, tests on all backends, demo, package (Linux + macOS)
- `tools/*.py` regenerate every reference vector (needs `pip install geographiclib pyproj mgrs coord-convert`)
- module metadata in `moon.mod` (the current toolchain deprecates `moon.mod.json` and migrates it to `moon.mod`)

## Left / ideas for further commits
- [ ] CLI (`geod inverse …`, `coord wgs2gcj …`) as a second example
- [ ] rhumb lines; geodesic intersections; `GeodesicExact` for |f| > 0.02
- [ ] optional finer China boundary polygon for `out_of_china` (currently the usual bounding box)
- [ ] benchmarks (`moon bench`)
- [ ] docs: README.mbt.md doc-tests, API docs on mooncakes

## Needs the user (not done by the agent)
- create public GitHub repo `ChrisTu117/moongeodesy`, then push (`git remote add origin … && git push -u origin main`)
- `moon login` (GitHub auth) and `moon publish` — publishing is irreversible, needs approval
- hackathon registration (Feishu form) and WeChat group (nickname = ChrisTu117)
