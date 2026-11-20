# spatialfamily: checklist for testing in jamovi

Menu: Spatial. Install: `jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`.
Open sample data from File > Open > Data Library > Spatial Statistics Family: "Columbus, Ohio: neighbourhood crime" (49 neighbourhoods) and "North Carolina: SIDS by county" (100 counties).
The example files are addressed as `sample:columbus_edges.csv`, `sample:columbus.geojson`, `sample:nc.shp`, `sample:nc.geojson`. In the Data tab check that `id` / `fips` are numbers or Nominal; either works.

## 1. Spatial Weights
Columbus: Unit id `id`; Spatial weights from "Neighbour list file", File path `sample:columbus_edges.csv`.
Expect: 49 units, 230 links, mean 4.69 neighbours (min 2, max 10), 0 without neighbours, 1 component, symmetric yes. The connectivity map needs polygons or coordinates, so with a neighbour list a note says so.
Switch to "Boundary file", file `sample:columbus.geojson`, Id field `id`, Neighbours defined by "Queen contiguity": 236 links; the map shows polygons with orange links. "Rook contiguity": 200 links.
Switch to "Coordinates": X `X`, Y `Y`, k nearest neighbours 5 -> 245 links, symmetric no; tick Make links symmetric -> 302 links, symmetric yes.
North Carolina: Boundary file `sample:nc.shp`, Id field `fips`, id variable `fips`, Queen: 100 units, 490 links, mean 4.9 (2 to 9). Order 2 with "Include lower orders" ticked: more links.
Messages: a wrong Id field names the available fields; a missing file says it was not found; k = 60 for 49 units says the number must be between 1 and 48.

## 2. Global autocorrelation
Variables `CRIME`, weights as in 1 (neighbour list). Expect Moran's I 0.4858, expected -0.0208, z 5.343, p < .001. Tick Geary's C: 0.5478, z 4.567. Tick Moran scatterplot: points with a red line.
Inference "Permutation": a note about 999 permutations; p about .001. Alternative "Negative autocorrelation": p close to 1.
With the polygons (queen): I 0.5002, z 5.589.

## 3. Local autocorrelation
Variable `CRIME`, neighbour list. Class counts: High-High 10, Low-Low 3, Not significant 36 (permutation, so counts can differ by one unit between seeds). The cluster map needs polygons or coordinates (use the Boundary file `sample:columbus.geojson`): red and blue clusters, grey elsewhere.
Statistic "Getis-Ord Gi*": Hot spot 10, Cold spot 3. Multiple testing "False discovery rate": fewer significant units. Tick Table of all units: 49 rows.
Save to the dataset: tick the three outputs; new columns appear (cluster class 0 to 4, statistic, p-value) and keep the row alignment if you filter some rows out first.

## 4. Spatial regression
Dependent `CRIME`, Covariates `INC`, `HOVAL`, neighbour list `sample:columbus_edges.csv`, model "Spatial lag (SAR)".
Expect: INC -1.0735, HOVAL -0.2700, rho 0.4039 (se 0.1277, z 3.16); log-likelihood -183.168; AIC 376.34, BIC 385.80; LR test 8.42 (1 df, p .004).
LM tests: LM error 4.611 (p .032), LM lag 7.856 (p .005), robust LM error 0.034 (p .855), robust LM lag 3.278 (p .070), SARMA 7.889 (2 df).
Effects: INC direct -1.1225, indirect -0.6784, total -1.8009; HOVAL -0.2823, -0.1706, -0.4529 (the point estimates equal spatialreg::impacts; the simulated standard errors vary with the seed).
OLS: INC -1.5973, HOVAL -0.2739, intercept 68.619. SEM: lambda 0.5209. SDM: rho 0.3825 and coefficients W:INC, W:HOVAL. Tick "Compare all models": seven rows (OLS, SLX, SAR, SEM, SDM, SDEM, SAC); AIC lowest for SAR (376.3).
Move `HOVAL` to "Covariates not to be lagged" and choose SDM: only W:INC appears.
Save fitted values and residuals: two new columns.
Messages: no covariate -> a hint; a constant covariate -> an error naming it.

## 5. Spatial Panel Regression (new in 0.2.0)
Dataset "Spatial panel: simulated, Columbus layout": dependent `y`, covariates `x1`, `x2`, unit id `id`, time period `time`, Neighbour list file `sample:columbus_edges.csv`, Unit fixed effects, model SAR.
Expect: x1 1.419, x2 -0.606, rho 0.862 (agrees with splm and Stata spxtregress, fe); the notes say 49 units x 6 periods = 294 observations. Two-way fixed effects: rho about 0.29. SDM adds W:x1, W:x2.

## 6. Geographically Weighted Regression
Columbus: dependent `CRIME`, covariates `INC`, `HOVAL`, X `X`, Y `Y`; Bandwidth "Set manually" 25 (adaptive, bi-square): AICc 380.22 (OLS 383.66), 16.46 effective parameters, R2 0.785; "Select by AICc" chooses 24. Maps of INC and HOVAL, the intercept map, and the table of local estimates (tick it). Save local R2 and residuals.

## 7. Industrial Agglomeration
Dataset "Industrial agglomeration": region `region`, industry `industry`, employment `employment`, establishments `establishments`. EG gamma: Autos 0.0866, Textiles 0.0766, Software 0.0727, Mining 0.0102, Food 0.0048, Retail 0.0041. Without establishments gamma is computed with H = 0 (a note says so). Save the location quotient: one new column.

## 8. Gravity Model (PPML)
Dataset "Gravity: simulated bilateral trade": flow `trade`, covariates `ldist`, `sameRegion`, `language`, origin `origin`, destination `destination`, origin and destination fixed effects, robust errors.
Expect: ldist -0.9207, sameRegion 0.5602, language 0.4731, pseudo log-likelihood -472.34. Untick both fixed effects: a constant appears. A covariate that is constant within the origin (e.g. a country code) is rejected with a message.

## 9. Spatial DID (Spillovers)
Dataset "Spatial DID: simulated policy with spillovers": outcome `outcome`, treatment `treat`, covariate `x1`, period `year`, unit id `unit`, X `cx`, Y `cy`, rings `10, 20`.
Expect: Treatment 2.125, Ring 0-10 0.873, Ring 10-20 0.137, x1 0.444 (clustered by unit, 99 df). Standard errors "Conley spatial HAC" needs a bandwidth (try 25). Without rings or the lagged treatment a message asks for one of them.
