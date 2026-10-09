# moongeodesy progress

Updated: 2026-10-09 (Asia/Shanghai)

## Done
- [x] Scaffold: moon.mod (mooncakes metadata), MIT LICENSE, .gitignore
- [x] `cn`: WGS-84 ↔ GCJ-02 ↔ BD-09 (forward, exact iterative inverse, approximate inverse), out_of_china; 6 tests vs 200 coord-convert reference points
- [x] `geodesic`: Karney inverse/direct, GeodesicLine (position, arc_position, densify), inverse_line/direct_line
  - GeodTest-short.dat (1000 cases): max |Δs12| 7.5 nm, direct position error 8.6 nm
  - cross-check vs Python geographiclib (1014 inverse + 1006 direct)
## In progress
- [ ] PolygonArea (geodesic polygon area/perimeter)
## Todo
- [ ] UTM/UPS (transverse Mercator, Krüger series) + MGRS
- [ ] README (EN + 中文), examples (cmd/), GitHub Actions CI
- [ ] `moon publish` (needs user login; not done)
