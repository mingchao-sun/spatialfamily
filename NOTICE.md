# Notices

All estimators are written from the published formulas (Cliff and Ord, 1981; Anselin, 1988, 1995; Ord and Getis, 1995; Ord, 1975; Anselin et al., 1996;
LeSage and Pace, 2009) in base R and the Matrix package, and are checked against the R packages `spdep`, `spatialreg`, `splm`, `GWmodel`, `spgwr` (GPL) and `fixest`, and against Stata 19 (`spregress`, `spxtregress`, `ppmlhdfe`, `reghdfe`), none of which are included.
The Columbus data come from the `spdep` / `spData` packages (Anselin, 1988) and the North Carolina SIDS data from the `sf` package (Cressie, 1993).

The example data for agglomeration, gravity, the spatial panel and the spatial DID are simulated for this module (`dev/make_data2.R`, `tests/panel_ref.R`).
