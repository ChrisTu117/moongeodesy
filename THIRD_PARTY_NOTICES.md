# Third-party notices

moongeodesy is MIT licensed (see [LICENSE](LICENSE)). It contains work
derived from, or data taken from, the sources below. Their notices and
terms follow.

## GeographicLib (MIT)

The `geodesic`, `utm`, and `conic` packages and parts of `internal/geomath`
are MoonBit ports of algorithms and code structure from GeographicLib by
Charles F. F. Karney, https://geographiclib.sourceforge.io/ (C++ and Python,
version 2.x). GeographicLib's licence:

```
The MIT License (MIT).

Copyright (c) 2008-2023, Charles Karney

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
```

## World Magnetic Model 2025 (public domain)

`src/magnetic/wmm2025_data.mbt` holds the WMM2025 coefficients and
`src/magnetic/wmm2025_testvalues_data_test.mbt` NOAA's test values, from
NOAA NCEI and the British Geological Survey
(https://www.ncei.noaa.gov/products/world-magnetic-model). Notice required by
17 U.S.C. 403: this is U.S. Government material, not subject to copyright
protection. NOAA NCEI Geomagnetic Modeling Team; British Geological Survey.
2024: World Magnetic Model 2025. https://doi.org/10.25921/aqfd-sd83

## Natural Earth (public domain)

`src/cn/boundary_data.mbt` is derived from Natural Earth 1:10m admin-0
countries (China point of view) and minor islands, v5.1.2,
https://www.naturalearthdata.com/. Natural Earth is in the public domain.

## Test data

These files are used only by tests. They are not part of the library's API.

- `src/geodesic/geodtest_data_test.mbt` and `src/examples/s3_accuracy/data.mbt`:
  rows of GeodTest-short.dat, C. F. F. Karney (2010),
  https://doi.org/10.5281/zenodo.32156, CC0 1.0 (public domain dedication).
- `src/cn/towns_data_test.mbt`: town coordinates from GeoNames,
  https://www.geonames.org/, licensed under Creative Commons Attribution 4.0
  (https://creativecommons.org/licenses/by/4.0/). Attribution: GeoNames.
  The file keeps only latitude and longitude (and three names).
- `src/conic/gigs_data_test.mbt`: points from the IOGP Geospatial Integrity
  of Geoscience Software (GIGS) Test Dataset v2.1, conversion tests 5102,
  5103, and 5109, © International Association of Oil & Gas Producers (IOGP),
  https://github.com/IOGP-GIGS/GIGSTestDataset. The GIGS Software is
  published by IOGP at no charge and may be used, copied, and distributed
  under the GIGS Terms of Use v2.1 (https://gigs.iogp.org/TermsOfUse),
  which require acknowledging IOGP's ownership, forbid distributing the GIGS
  Software itself for profit, and allow it in commercial packages whose
  value lies in what the provider adds. This file is under those terms, not
  under this project's MIT licence; the data is provided "as is" without
  warranty, as stated in the Terms of Use.
- The remaining reference vectors were computed by the scripts in `tools/`
  with Python geographiclib (MIT), PROJ via pyproj (MIT), the `mgrs` package,
  `coord-convert` (MIT), GeographicLib's command-line tools (MIT), and
  mpmath (BSD). They are numbers produced by those programs, not copies of
  their code.

## Formulas

The GCJ-02 and BD-09 transforms use the widely published reverse-engineered
formulas, as in `coord-convert` (MIT); the code here is our own. MGRS letter
tables follow NGA GeoTrans and GeographicLib's `MGRS`.
