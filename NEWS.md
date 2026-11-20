# spatialfamily 1.0.0

First release. Nine analyses in the jamovi menu Spatial:

- **Spatial Weights**: weights from coordinates (k nearest neighbours, distance band, inverse distance), polygon shapefiles and GeoJSON (queen, rook, centroid based) or neighbour lists (edge list, GAL, GWT); higher orders, symmetrising, styles W, B and C; connectivity map.
- **Global (Moran, Geary, G)** and **Local (LISA, Gi\*)** spatial autocorrelation with analytical and permutation inference, Moran scatterplots, cluster maps, multiple-testing adjustment, results saved to the dataset.
- **Spatial Regression**: OLS, SLX, SAR, SEM, SDM, SDEM and SAC by maximum likelihood, LM diagnostics, direct / indirect / total effects, model comparison.
- **Spatial Panel Regression**: unit, time or two-way fixed effects with SAR, SEM, SDM, SDEM and SLX.
- **Geographically Weighted Regression**: five kernels, adaptive or fixed bandwidth selected by AICc or cross-validation, local coefficients, t values and maps.
- **Industrial Agglomeration**: locational Gini, Herfindahl, Ellison-Glaeser gamma, Krugman specialisation index, location quotients.
- **Gravity Model (PPML)**: Poisson pseudo maximum likelihood with high-dimensional fixed effects and multiway clustering.
- **Spatial DID (Spillovers)**: two-way fixed effects DID with distance rings, spatially lagged treatment and Conley standard errors.

Example data: Columbus (cross-section and a simulated panel), North Carolina SIDS, simulated region-industry employment, bilateral trade and a spillover DID.
Validated against spdep, spatialreg, splm, GWmodel, sf and fixest in R and against Stata 19 (spregress, spxtregress, ppmlhdfe, reghdfe).
