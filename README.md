# spatialfamily — spatial statistics and spatial econometrics in jamovi

A jamovi module for the spatial analysis that urban, regional and spatial economists use: build a spatial weights matrix from coordinates, a boundary file or a neighbour list,
test and map spatial clustering, and estimate spatial regression models by maximum likelihood. Chinese summary: [README.zh.md](README.zh.md). Japanese summary: [README.ja.md](README.ja.md).

## Analyses (jamovi menu: Spatial)

| Analysis | What it does | References |
|---|---|---|
| **Spatial Weights** | Builds W and describes it: links, neighbours per unit, isolated units, connected components, connectivity map, histogram of neighbour counts. Sources: coordinate columns (k nearest neighbours, distance band, inverse distance; planar or longitude/latitude), a polygon shapefile or GeoJSON (queen or rook contiguity, or centroid-based), a neighbour list (edge list csv/txt, GAL, GWT). Higher-order neighbours, symmetrising, row-standardised / binary / global style. | Anselin (1988) |
| **Global (Moran, Geary, G)** | Moran's I, Geary's C and the Getis-Ord general G for several variables at once; randomisation, normality or permutation inference; Moran scatterplot. | Cliff and Ord (1981) |
| **Local (LISA, Gi\*)** | Local Moran's I with High-High / Low-Low / High-Low / Low-High classes, Getis-Ord Gi and Gi\*; conditional permutation p-values, FDR or Bonferroni adjustment; cluster map (polygons or points), Moran scatterplot, table of units; saves the class, statistic and p-value to the dataset. | Anselin (1995); Ord and Getis (1995) |
| **Spatial Regression** | OLS, SLX, SAR (lag), SEM (error), SDM (Durbin), SDEM and SAC by maximum likelihood; LM tests (error, lag, robust, SARMA) on the OLS residuals; average direct, indirect and total effects with simulated standard errors; Moran's I of the residuals; comparison of all seven models by log-likelihood, AIC and BIC; saves fitted values and residuals. | Ord (1975); Anselin et al. (1996); LeSage and Pace (2009) |
| **Spatial Panel Regression** | Balanced panels: unit, time or two-way fixed effects (or pooled) with SAR, SEM, SDM, SDEM or SLX by maximum likelihood (within transformation, T ln\|I - rho W\|), direct / indirect / total effects. | Elhorst (2003, 2014) |
| **Geographically Weighted Regression** | Local coefficients that vary over space: bi-square, Gaussian, exponential, tricube or box-car kernels, adaptive or fixed bandwidth selected by AICc or cross-validation, local standard errors and t values, local R2, maps of the coefficients, AICc and effective parameters against OLS. | Brunsdon et al. (1996) |
| **Industrial Agglomeration** | From region x industry employment: locational Gini, Herfindahl, Ellison-Glaeser gamma (with establishment counts), Krugman specialisation index of regions, location quotients (saved to the dataset). | Krugman (1991); Ellison and Glaeser (1997) |
| **Gravity Model (PPML)** | Poisson pseudo-ML for bilateral flows with zeros, high-dimensional fixed effects (origin, destination, pair, origin x period, destination x period, period) absorbed by alternating projections, robust, clustered and two-way clustered standard errors. | Santos Silva and Tenreyro (2006); Correia et al. (2020) |
| **Spatial DID (Spillovers)** | Two-way fixed effects DID with distance rings around the treated units, a spatially lagged treatment, covariates, and standard errors clustered by unit, robust or Conley spatial HAC. | Conley (1999); Butts (2021) |

Every analysis has *Method notes* and a *Data and setup guide* in its Output section.

## Where the weights come from
- **Coordinates**: two numeric columns. Tick *longitude / latitude* for degrees (great-circle distance in km).
- **Boundary file**: a polygon shapefile (`.shp` + `.dbf` in the same folder) or GeoJSON. Type the path, choose the id column in the data and type the name of the id field in the file.
- **Neighbour list file**: one link per line (`from,to` or `from,to,weight`), or a GAL / GWT file. The ids must match the id column in the data.
- Example files installed with the module are addressed as `sample:columbus.geojson`, `sample:columbus_edges.csv`, `sample:nc.shp`, `sample:nc.geojson`.

## Install
- **jamovi library**: in jamovi open the modules button (top right), *jamovi library*, find *Spatial Statistics Family* and install (once the module is listed).
- **Sideload**: download the `.jmo` file and use *Install from file* in the same menu. A `.jmo` works only on the operating system, processor type and jamovi version it was built for.
- **From source** (needs R and the jamovi compiler `jmvtools`; adjust the path to your jamovi folder):
```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

## Use from R (Rj editor, scripts)
Every analysis is also an R function. Pass `NULL` for the arguments you do not use, because they have no default values. For example, after installing the module,
```r
spatialfamily::moran(data = d, vars = "CRIME", idVar = NULL, xVar = "X", yVar = "Y")
```

## Data format
One row per spatial unit. Rows with missing values in the analysis variables are dropped and reported. For file-based weights the data need an id column that matches the ids in the file.

## Example data
Shown in jamovi under **Open > Data Library > Spatial Statistics Family** (also `columbus_panel.csv` for the spatial panel, `region_industry.csv`, `gravity_sim.csv` and `spillover_did.csv`, all simulated). `columbus.csv`: 49 Columbus, Ohio neighbourhoods (Anselin, 1988; `CRIME`, `INC`, `HOVAL`, coordinates `X`, `Y`);
`nc_sids.csv`: 100 North Carolina counties (Cressie, 1993; SIDS rates, centroids `lon`, `lat`, id `fips`). The polygons `columbus.geojson`, `nc.shp` (with `.dbf`, `.shx`, `.prj`),
`nc.geojson` and the neighbour list `columbus_edges.csv` are installed next to them; each dataset entry says which variable goes where.

## Validation
The engines use base R and Matrix only. `tests/reference.R` compares them with `spdep` and `spatialreg`: Moran's I, Geary's C and general G with their variances and p-values, local Moran
(Ii, E, Var, z), Gi and Gi\*, the LM tests, and the maximum likelihood fits of SAR, SEM, SDM, SDEM, SAC and SLX (log-likelihood, coefficients, rho, lambda and observed-information standard errors)
and direct / indirect effects agree to the printed precision. `tests/weights_ref.R` checks shapefile and GeoJSON readers, queen / rook contiguity, higher orders, k nearest neighbours, distance bands
and the GAL reader against `spdep` and `sf`. `tests/panel_ref.R` (splm), `tests/gwr_ref.R` (GWmodel: five kernels, adaptive and fixed bandwidths; coefficients, standard errors, AICc), `tests/agg_gravity_ref.R` (hand-computed indices; PPML against fixest::fepois including two-way clustered standard errors) and `tests/did_ref.R` (fixest::feols; the Conley matrix against its definition by explicit loops) cover the other analyses.
`../stata_validation/07_spatial.do` and `compare_spatial.R` compare with Stata 19: `spregress` (SAR, SEM, SDM, SAC with standard errors and impacts), `estat moran` (the LM error test), `spxtregress, fe`, `ppmlhdfe` and `reghdfe`; all agree.
`tests/smoke.R` runs all analyses headless, including bad input; `tests/perf.R` times larger samples (5,000 units: SAR in about 2 seconds).

## Known limitations
- Spatial panel models need a balanced panel and fixed effects (no random effects, no dynamic panels, no Lee-Yu bias correction). Geographically weighted regression has no Monte Carlo test of non-stationarity and is limited to 4,000 units; GWR maps are point maps.
- Not included: spatial RD (use the geographic RD of the `rdfamily` module), the Duranton-Overman distance-based concentration measure (needs plant coordinates), multiscale GWR, spatial multilevel models.
- The Ellison-Glaeser index approximates the plant Herfindahl from establishment counts; plant-level data give exact values.
- Spatial regression with more than 1,500 units needs row-standardised weights (the log-determinant then comes from a sparse decomposition and the direct effects from a stochastic trace approximation, error below 1%).
- Standard errors are the observed-information (numerical Hessian) ones; spatialreg's default reports the asymptotic ones, which differ slightly in small samples.
- jamovi has no interactive map; maps are static images. Shapefiles must be polygons; point and line shapefiles and other formats (GeoPackage, KML) are not read.
- The total effect with non-row-standardised weights follows LeSage and Pace (the mean row sum of the multiplier matrix); spatialreg reports a different quantity for binary weights.

## Notes for developers
- Layout: `R/engine_*.R` is plain R with no jamovi dependency; `R/*.b.R` and `jamovi/*.yaml` are the jamovi front end. The yaml files are generated by `dev/gen_yaml.py` and its parts `dev/gen_more*.py` (the weights options are shared); edit the scripts, not the yaml. `dev/make_data2.R` rebuilds the simulated example data.
- Tests: `Rscript tests/smoke.R`, `tests/reference.R`, `tests/weights_ref.R` (the last two need `spdep`, `spatialreg`, `sf`); `tests/ui_checklist.md` lists what to check by hand in jamovi. `dev/make_data.R` rebuilds the example data.
- Ribbon icon: `icon/analysis-spatial.svg` (74 x 44): a small map of regions (high values orange, low values blue) with neighbour links. jmvtools 28.17 cannot attach an icon to a module, so `dev/icon_preview.ps1` adds one CSS rule for the menu group to jamovi's client stylesheet (administrator PowerShell, restart jamovi afterwards).

## How to cite
Sun, M. *spatialfamily: Spatial Statistics Family for jamovi* (jamovi module). https://github.com/mingchao-sun/spatialfamily. Please also cite the papers behind the methods you use (listed in the output).

## License
GPL-3, as stated in `DESCRIPTION`; third-party data are listed in `NOTICE.md`. Bug reports and suggestions: https://github.com/mingchao-sun/spatialfamily/issues
