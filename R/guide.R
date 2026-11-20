# "Data and setup guide" and "Method notes" for each analysis.

sp_settings <- function(rows) {
    paste0("<table style='border-collapse:collapse;margin:4px 0 8px 0'>",
           "<tr><th style='text-align:left;padding:2px 14px 2px 0;border-bottom:1px solid #999'>Field</th>",
           "<th style='text-align:left;padding:2px 0;border-bottom:1px solid #999'>Set to</th></tr>",
           paste0("<tr><td style='padding:2px 14px 2px 0'>", names(rows), "</td><td style='padding:2px 0'>",
                  unlist(rows), "</td></tr>", collapse = ""), "</table>")
}

sp_weights_help <- paste0(
    "<h4>Where the neighbours come from</h4><ul>",
    "<li><b>Coordinates</b>: two numeric columns (x / longitude and y / latitude). Choose nearest neighbours (k), a distance band, or inverse distance. ",
    "Tick <i>Coordinates are longitude / latitude</i> for degrees (distances are then great-circle kilometres); leave it unticked for projected coordinates.</li>",
    "<li><b>Boundary file</b>: a polygon shapefile (<code>.shp</code>, with its <code>.dbf</code> in the same folder) or a GeoJSON file. Enter the full path, ",
    "choose the data column that holds the unit id, and type the name of the field in the file that holds the same id (for example <code>fips</code>). ",
    "Queen contiguity (shared border or corner) and rook contiguity (shared border) are built from the polygons; k nearest neighbours, a distance band or inverse distance use the polygon centroids.</li>",
    "<li><b>Neighbour list file</b>: a text file with one link per line, <code>from,to</code> or <code>from,to,weight</code> (a header line is allowed), or a GeoDa/spdep <code>.gal</code> or <code>.gwt</code> file. ",
    "The ids in the file must be the same codes as the id column in the data.</li></ul>",
    "<p>Type or paste the full path of the file on your computer. The easiest way: in File Explorer hold Shift, right-click the file, choose <i>Copy as path</i>, and paste it into the File path box (quotation marks are removed; <code>%APPDATA%</code> style variables work). The example files that come with the module can be addressed as <code>sample:columbus.geojson</code>, <code>sample:columbus_edges.csv</code>, <code>sample:nc.shp</code> or <code>sample:nc.geojson</code>. The data must have one row per unit. Rows with missing values in the analysis variables are dropped before the weights are built (for file-based weights the unit is removed from the neighbour structure).</p>",
    "<h4>Style</h4><p><i>Row-standardised</i> (W) makes every row sum to 1, so the spatial lag is a weighted mean of the neighbours. <i>As built</i> (B) keeps binary links or distance weights. ",
    "<i>Globally standardised</i> (C) rescales the weights so they sum to the number of units. <i>Order</i> above 1 links units that are that many steps apart (tick <i>cumulative</i> to include the closer neighbours).</p>")

sp_guide_raw <- function(key) {
    switch(key,
    swm = paste0(
        "<h3>Spatial Weights: data and setup</h3>",
        "<p>Builds the spatial weights matrix and describes it: how many neighbours each unit has, whether any unit is isolated, and whether the links form one connected map. Use it to check the weights before any spatial statistic or model.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>North Carolina SIDS (100 counties)</b>", "Unit id" = "fips", "Spatial weights from" = "Coordinates", "X / Y" = "lon / lat (tick longitude / latitude)",
                         "Method" = "5 nearest neighbours")),
        "<p>Or use the county polygons: weights from <b>Boundary file</b>, file <code>sample:nc.shp</code>, id field <code>fips</code>, Queen contiguity. Expected: 100 units, 490 links (average 4.9 neighbours, 2 to 9), one connected component.</p>",
        sp_weights_help),
    moran = paste0(
        "<h3>Global Spatial Autocorrelation: data and setup</h3>",
        "<p>Tests whether a variable is spatially clustered (positive autocorrelation) or dispersed (negative): Moran's I, Geary's C and the Getis-Ord general G.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Columbus, Ohio (49 neighbourhoods)</b>", "Variables" = "CRIME", "Unit id" = "id", "Spatial weights from" = "Neighbour list file <code>sample:columbus_edges.csv</code> (or Boundary file <code>sample:columbus.geojson</code>, id field <code>id</code>, queen contiguity; or coordinates X, Y)")),
        "<p>Expected with the neighbour list <code>columbus_edges.csv</code> (the textbook weights, 230 links) and row-standardised weights: Moran's I = 0.4858 (expected -0.0208), z = 5.34, p &lt; 0.001 (randomisation). With the polygons in <code>columbus.geojson</code> (queen contiguity, 236 links) I = 0.5002, z = 5.59.</p>",
        "<h4>Options</h4><ul><li><b>Inference</b>: randomisation (the default in most software), normality, or a permutation test (the variable is shuffled over the units; set the number of permutations and a seed).</li>",
        "<li><b>Alternative</b>: two-sided, or one-sided in the direction of clustering (greater) or dispersion (less).</li>",
        "<li><b>Getis-Ord G</b> needs a non-negative variable and tests whether high values or low values cluster; use it with binary or distance weights.</li>",
        "<li><b>Moran scatterplot</b>: each unit's standardised value against the average of its neighbours; the line's slope is Moran's I for row-standardised weights.</li></ul>",
        sp_weights_help),
    lisa = paste0(
        "<h3>Local Spatial Autocorrelation: data and setup</h3>",
        "<p>Finds where the clusters are: local Moran's I (LISA) classifies each unit as High-High, Low-Low, High-Low or Low-High and Getis-Ord Gi* marks hot and cold spots. Significance comes from conditional permutation.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Columbus, Ohio (49 neighbourhoods)</b>", "Variable" = "CRIME", "Unit id" = "id", "Spatial weights from" = "Boundary file <code>sample:columbus.geojson</code>, id field <code>id</code>, Queen contiguity", "Statistic" = "Local Moran's I")),
        "<p>The map needs polygons (boundary file) or coordinates; with a neighbour list file only the table and the Moran scatterplot are available. Tick the Output options to add the cluster class (0 not significant, 1 High-High, 2 Low-Low, 3 Low-High, 4 High-Low), the statistic and the p-value to the dataset.</p>",
        "<p>With many units some will be significant by chance; use the false discovery rate or Bonferroni adjustment for the cluster classes.</p>",
        sp_weights_help),
    sreg = paste0(
        "<h3>Spatial Regression: data and setup</h3>",
        "<p>Linear regression with spatial dependence, estimated by maximum likelihood: spatial lag (SAR), spatial error (SEM), spatial Durbin (SDM), spatial Durbin error (SDEM), spatial lag with spatial error (SAC) and spatially lagged covariates (SLX), with ordinary least squares as the benchmark.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Columbus, Ohio (49 neighbourhoods)</b>", "Dependent" = "CRIME", "Covariates" = "INC, HOVAL", "Unit id" = "id", "Spatial weights from" = "Neighbour list file <code>sample:columbus_edges.csv</code> (or Boundary file <code>sample:columbus.geojson</code>, id field <code>id</code>, queen)", "Model" = "Spatial lag (SAR)")),
        "<p>Expected for SAR with the neighbour list <code>columbus_edges.csv</code> and row-standardised weights: rho = 0.4039, INC = -1.0735, HOVAL = -0.2700, log-likelihood -183.168 (identical to R's spatialreg::lagsarlm). With the polygons (queen contiguity) rho = 0.4233, INC = -1.0487, HOVAL = -0.2663, log-likelihood -182.674.</p>",
        "<h4>Reading the output</h4><ul>",
        "<li><b>Diagnostics for OLS</b>: Lagrange multiplier tests of spatial error and spatial lag dependence in the OLS residuals, and their robust versions. If only one robust test is significant, choose that model; if both are, a Durbin model is the usual next step.</li>",
        "<li><b>Effects</b>: with a spatial lag, a change in a covariate in unit i spills over to the others, so the coefficient is not the marginal effect. The table gives the average direct, indirect (spillover) and total effects with simulated standard errors. Spatial lagged covariates (SLX, SDEM) have direct = own coefficient and indirect = the coefficient of the lagged variable.</li>",
        "<li><b>Compare models</b>: log-likelihood, AIC and BIC of all models on the same weights.</li></ul>",
        "<p>Regressors listed in <i>Covariates</i> are lagged (WX) in the SLX, SDM and SDEM models when they are numeric; regressors in <i>Covariates not to be lagged</i> enter only as themselves. A constant is always included.</p>",
        sp_weights_help),
    span = paste0(
        "<h3>Spatial Panel Regression: data and setup</h3>",
        "<p>Panel regression with spatial dependence between units: unit and/or time fixed effects with a spatial lag (SAR), spatial error (SEM), spatial Durbin (SDM), spatial Durbin error (SDEM) or lagged covariates (SLX), estimated by maximum likelihood.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Spatial panel: simulated, Columbus layout</b>", "Dependent" = "y", "Covariates" = "x1, x2", "Unit id" = "id", "Time period" = "time",
                         "Spatial weights from" = "Neighbour list file <code>sample:columbus_edges.csv</code>", "Fixed effects" = "Unit fixed effects", "Model" = "Spatial lag (SAR)")),
        "<p>The data were simulated with rho = 0.4, coefficients 1.5 and -0.8, unit and time effects. Expected: rho and the coefficients close to those values (the estimates agree with R's splm::spml and Stata's spxtregress).</p>",
        "<h4>Data layout</h4><p>Long format: one row per unit and period, a unit id and a period variable. The spatial weights describe the N units (not the NT rows) and are matched through the unit id; with coordinates the first row of each unit gives its location, so coordinates must not change over time. The panel must be balanced: units that are not observed in every period (after dropping missing values) are removed and counted.</p>",
        "<h4>Options</h4><ul><li><b>Fixed effects</b>: unit effects (the usual choice), time effects, both, or none (pooled with a constant). A covariate that does not vary within units cannot be estimated with unit effects.</li>",
        "<li><b>Lagged covariates</b>: numeric covariates in the first box are lagged in SLX, SDM and SDEM; those in the second box are not.</li>",
        "<li><b>Effects</b>: average direct, indirect and total effects as in the cross-section (SAR, SDM, SLX, SDEM).</li></ul>",
        sp_weights_help),
    gwr = paste0(
        "<h3>Geographically Weighted Regression: data and setup</h3>",
        "<p>Estimates a separate regression around every unit, giving more weight to nearby units, so the coefficients can vary over space. Use it to see where a relationship is stronger or weaker or changes sign.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Columbus, Ohio (49 neighbourhoods)</b>", "Dependent" = "CRIME", "Covariates" = "INC, HOVAL", "X / Y" = "X / Y", "Kernel" = "Bi-square, adaptive, bandwidth set manually to 25")),
        "<p>Expected with a manual adaptive bandwidth of 25 (identical to GWmodel::gwr.basic): AICc 380.22 (OLS 383.66), 16.46 effective parameters, R-squared 0.785 (OLS 0.552); the local coefficient of INC runs from -2.71 to +1.49. With the bandwidth selected by AICc, 24 nearest neighbours are chosen.</p>",
        "<h4>Data layout</h4><p>One row per location with numeric coordinates (x, y or longitude, latitude; tick the longitude/latitude option for degrees, distances are then great-circle kilometres). Rows with missing values are dropped. At most 4000 units.</p>",
        "<h4>Options</h4><ul><li><b>Kernel</b>: how weight declines with distance. Bi-square, tricube and box-car give zero weight beyond the bandwidth, Gaussian and exponential never reach zero.</li>",
        "<li><b>Bandwidth</b>: adaptive = the number of nearest neighbours that enter each local regression (the same number everywhere, so it adapts to sparse and dense areas), fixed = a distance. Select by AICc (the usual choice) or cross-validation, or set it by hand.</li>",
        "<li><b>Maps</b>: the local coefficient (or its t value) at each unit.</li></ul>"),
    agg = paste0(
        "<h3>Industrial Agglomeration: data and setup</h3>",
        "<p>Measures how strongly each industry is concentrated in space (Ellison-Glaeser gamma, locational Gini, Herfindahl) and how specialised each region is (Krugman index, Herfindahl, entropy), from a table of employment by region and industry.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Industrial agglomeration: simulated regions and industries</b>", "Region" = "region", "Industry" = "industry", "Employment" = "employment", "Establishments" = "establishments")),
        "<p>Expected: the Ellison-Glaeser gamma is highest for Autos (0.0866), then Textiles (0.0766) and Software (0.0727); Mining is 0.0102 because its few large plants explain much of its concentration (H = 0.050); Food (0.0048) and Retail (0.0041) follow total employment.</p>",
        "<h4>Data layout</h4><p>Long format: one row per region and industry with the employment (or output, or value added) and, optionally, the number of establishments. Region-industry combinations that are missing count as zero; repeated rows are summed. Rows with missing values are dropped.</p>",
        "<h4>Reading the indices</h4><ul><li><b>Raw concentration G</b>: sum over regions of (industry share - region share of total employment)<sup>2</sup>.</li>",
        "<li><b>EG gamma</b> corrects G for the size of plants (Herfindahl H of the plant sizes) so that it is comparable across industries; above 0.05 an industry is highly localised, below 0.02 not very localised. Without establishment counts H is taken as 0.</li>",
        "<li><b>Locational Gini</b>: 0 when the industry is spread like total employment, close to 1 when it sits in few regions.</li>",
        "<li><b>Krugman index</b> of a region: how different its industry structure is from the rest of the economy.</li>",
        "<li><b>Save</b>: the location quotient of each row (industry share in the region divided by the industry's national share) can be added to the dataset.</li></ul>"),
    gravity = paste0(
        "<h3>Gravity Model (PPML): data and setup</h3>",
        "<p>Estimates a gravity equation for bilateral flows (trade, migration, commuting, capital) by Poisson pseudo maximum likelihood, the estimator recommended by Santos Silva and Tenreyro (2006): it keeps zero flows and is consistent under heteroskedasticity. Fixed effects for origin and destination (or origin-by-period and destination-by-period) absorb multilateral resistance.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Gravity: simulated bilateral trade of 20 countries</b>", "Flow" = "trade", "Covariates" = "ldist, sameRegion, language", "Origin" = "origin", "Destination" = "destination", "Fixed effects" = "Origin and destination")),
        "<p>Expected (robust standard errors): ldist -0.9207, sameRegion 0.5602, language 0.4731, pseudo log-likelihood -472.34; the data were simulated with -0.9, 0.6 and 0.4. The results agree with fixest::fepois.</p>",
        "<h4>Data layout</h4><p>One row per origin-destination pair (and period). The flow is a non-negative number in levels (do not take logs); covariates such as log distance enter as ordinary variables. Variables that are constant within an origin or destination (such as GDP) cannot be estimated when origin and destination fixed effects are included; use them only without those fixed effects or interacted with the partner.</p>",
        "<h4>Options</h4><ul><li><b>Fixed effects</b>: origin, destination, pair (for panels), origin x period, destination x period, period. A covariate that is absorbed by the fixed effects is reported as an error.</li>",
        "<li><b>Standard errors</b>: robust, clustered by origin, destination or pair, or two-way clustered (Cameron, Gelbach and Miller, 2011). Clusters with fewer than about 10 groups are unreliable.</li>",
        "<li>Groups whose flows are all zero (their fixed effect would be minus infinity) are dropped and counted.</li></ul>"),
    spdid = paste0(
        "<h3>Spatial DID (Spillovers): data and setup</h3>",
        "<p>A two-way fixed effects difference-in-differences that also asks whether the policy affected <i>untreated</i> units near the treated ones. Distance rings around the treated units (and/or the spatially lagged treatment) measure the spillovers, so that the control group is not contaminated by them.</p>",
        "<h4>Sample data</h4>",
        sp_settings(list("Dataset" = "<b>Spatial DID: simulated policy with spillovers</b>", "Outcome" = "outcome", "Treatment" = "treat", "Covariates" = "x1", "Time period" = "year", "Unit id" = "unit",
                         "Spatial weights from" = "Coordinates: X cx, Y cy", "Distance rings" = "10, 20")),
        "<p>Expected: treatment effect 2.125 (simulated 2), ring 0-10 effect 0.873 (simulated 1.0), ring 10-20 effect 0.137 (simulated 0.4, not significant with these 75 unit-periods), x1 0.444 (0.5); clustered by unit, 99 degrees of freedom.</p>",
        "<h4>Data layout</h4><p>Long format with a unit id, a period variable, the outcome, and a treatment indicator that is 1 in every period in which the unit is treated (0 before adoption). Locations come from coordinate columns (constant within a unit) or from a boundary file (polygon centroids). Distances are great-circle kilometres when the longitude/latitude option is ticked, otherwise in the units of the coordinates.</p>",
        "<h4>Options</h4><ul><li><b>Distance rings</b>: upper limits, e.g. <code>10, 20</code> gives the rings 0-10 and 10-20. In every period an untreated unit is placed in the ring of its nearest treated unit; units farther than the last limit form the control group.</li>",
        "<li><b>Spatially lagged treatment</b>: the weighted share of a unit's neighbours that are treated (any weights source, including a neighbour list).</li>",
        "<li><b>Standard errors</b>: clustered by unit (the default), robust, or Conley spatial HAC (distance bandwidth in the same units as the distances, Bartlett kernel, optional serial correlation lag; Conley, 1999; Hsiang, 2010).</li></ul>",
        sp_weights_help),
    stop("unknown guide: ", key))
}

sp_notes_raw <- function(key) {
    switch(key,
    swm = paste0(
        "<h3>Spatial Weights: method notes</h3>",
        "<p>The weights matrix W has a row and column for each unit; w<sub>ij</sub> &gt; 0 when j is a neighbour of i. The diagonal is zero. ",
        "Nearest-neighbour weights are not symmetric (i can be among j's neighbours without the reverse); tick <i>Make links symmetric</i> to add the reverse links.</p>",
        "<p><b>Distance band.</b> With the threshold left at 0 the smallest distance is used that leaves no unit without a neighbour. A large threshold in a dense area gives some units very many neighbours.</p>",
        "<p><b>Islands.</b> A unit without neighbours has no spatial lag; statistics ignore it. Spatial regression needs every unit to have at least one neighbour in practice, so add links or use nearest neighbours.</p>",
        "<p><b>Connected components.</b> If the links split the map into several components, spatial statistics describe each component separately; a single large component is usually what you want.</p>"),
    moran = paste0(
        "<h3>Global Spatial Autocorrelation: method notes</h3>",
        "<p><b>Moran's I</b> = (n / S0) &sum;&sum; w<sub>ij</sub> z<sub>i</sub> z<sub>j</sub> / &sum; z<sub>i</sub><sup>2</sup>, with z the deviation from the mean and S0 the sum of the weights. Under no spatial pattern E[I] = -1/(n-1). ",
        "Positive values mean that similar values are near each other.</p>",
        "<p><b>Geary's C</b> compares squared differences between neighbours (E[C] = 1); values below 1 mean positive autocorrelation. The reported z is signed so that positive z means positive autocorrelation for both statistics.</p>",
        "<p><b>Getis-Ord G</b> = &sum;&sum; w<sub>ij</sub> x<sub>i</sub> x<sub>j</sub> / &sum;&sum; x<sub>i</sub> x<sub>j</sub> (i &ne; j). A value above its expectation means that high values cluster.</p>",
        "<p><b>Inference.</b> Randomisation and normality use the analytical variances of Cliff and Ord (1981); the permutation test reshuffles the values over the units (pseudo p-value). ",
        "Results agree with spdep's moran.test, geary.test and globalG.test.</p>"),
    lisa = paste0(
        "<h3>Local Spatial Autocorrelation: method notes</h3>",
        "<p><b>Local Moran's I</b> (Anselin, 1995): I<sub>i</sub> = z<sub>i</sub> &sum;<sub>j</sub> w<sub>ij</sub> z<sub>j</sub> / m<sub>2</sub>. A unit is High-High when its value and its neighbours' average are both above the mean, Low-Low when both are below, and High-Low / Low-High when they differ (spatial outliers). ",
        "Only units whose permutation p-value is below the chosen level (after adjustment) are classified.</p>",
        "<p><b>Getis-Ord Gi*</b> (Ord and Getis, 1995) is the z-score of the sum of a unit's value and its neighbours' values (Gi leaves the unit's own value out). Large positive z marks hot spots, large negative z cold spots.</p>",
        "<p><b>Permutation inference.</b> For each unit its own value is held fixed and the other values are randomly assigned to its neighbours (conditional permutation); the p-value is (rank + 1) / (permutations + 1). ",
        "The analytical z-scores agree with spdep's localmoran (unconditional) and localG.</p>"),
    sreg = paste0(
        "<h3>Spatial Regression: method notes</h3>",
        "<p><b>Models.</b> y = &rho;Wy + X&beta; + WX&theta; + u, with u = &lambda;Wu + &epsilon;. SAR: only &rho;; SEM: only &lambda;; SDM: &rho; and &theta;; SDEM: &lambda; and &theta;; SAC: &rho; and &lambda;; SLX: only &theta; (estimated by OLS).</p>",
        "<p><b>Estimation.</b> Maximum likelihood (Ord, 1975): the likelihood is concentrated over &rho; and &lambda; and the log-determinant ln|I - &rho;W| is computed from the eigenvalues of W (up to 1500 units) or from a sparse LU decomposition on a grid (larger samples, row-standardised weights only). ",
        "&rho; and &lambda; are restricted to the range where I - &rho;W is invertible. Standard errors come from the numerical Hessian of the full log-likelihood (observed information); the estimates and log-likelihoods agree with spatialreg.</p>",
        "<p><b>Effects.</b> In a model with &rho; the derivative &part;y/&part;x<sub>k</sub> is the matrix (I - &rho;W)<sup>-1</sup>(&beta;<sub>k</sub>I + &theta;<sub>k</sub>W). The average direct effect is the mean of its diagonal, the total effect the mean row sum, the indirect effect their difference (LeSage and Pace, 2009). Standard errors come from simulating &rho; and the coefficients from their joint asymptotic normal distribution.</p>",
        "<p><b>Diagnostics.</b> LM tests of Anselin et al. (1996) use the OLS residuals; they test the null of no spatial dependence against error or lag alternatives.</p>"),
    span = paste0(
        "<h3>Spatial Panel Regression: method notes</h3>",
        "<p><b>Model.</b> y<sub>it</sub> = &rho; &sum;<sub>j</sub> w<sub>ij</sub> y<sub>jt</sub> + x<sub>it</sub>&beta; + &sum;<sub>j</sub> w<sub>ij</sub> x<sub>jt</sub>&theta; + &mu;<sub>i</sub> + &tau;<sub>t</sub> + &epsilon;<sub>it</sub> (SAR: &theta; = 0; SEM: spatially autocorrelated errors).</p>",
        "<p><b>Estimation.</b> Fixed effects are removed by the within transformation (unit means, time means or both) after the spatial lags are formed; the concentrated log-likelihood is -(NT/2) ln(&sigma;<sup>2</sup>) + T ln|I - &rho;W| (Elhorst, 2003, 2014), maximised over &rho; (or &lambda;). Standard errors come from the Hessian of the full log-likelihood, multiplied by NT / (NT - number of fixed effects) for the degrees of freedom used by the fixed effects (as Stata's spxtregress does). The bias corrections of Lee and Yu (2010) are not applied.</p>",
        "<p><b>Checks.</b> Results agree with spml in the R package splm for unit and two-way fixed effects (SAR, SEM, SDM) and with spatialreg for the pooled model.</p>"),
    gwr = paste0(
        "<h3>Geographically Weighted Regression: method notes</h3>",
        "<p>At each location i the model is fitted by weighted least squares with weights w<sub>ij</sub> = K(d<sub>ij</sub> / b<sub>i</sub>), giving &beta;(i) = (X'W<sub>i</sub>X)<sup>-1</sup>X'W<sub>i</sub>y. The hat matrix S has rows x<sub>i</sub>'(X'W<sub>i</sub>X)<sup>-1</sup>X'W<sub>i</sub>.</p>",
        "<p><b>Diagnostics.</b> AICc = 2n ln(&sigma;) + n ln(2&pi;) + n (n + tr S) / (n - 2 - tr S) with &sigma;<sup>2</sup> = RSS / n; standard errors use &sigma;<sup>2</sup> = RSS / (n - 2 tr S + tr S'S); the effective number of parameters is 2 tr S - tr S'S (GWmodel conventions; the results agree with GWmodel::gwr.basic).</p>",
        "<p><b>Cautions.</b> Local estimates are correlated and the local t tests are many; read the maps as exploratory evidence of non-stationarity, not as separate hypothesis tests. Collinearity between covariates can produce local coefficients with opposite signs that cancel. The local R-squared is 1 - weighted RSS / weighted TSS around the weighted mean.</p>"),
    agg = paste0(
        "<h3>Industrial Agglomeration: method notes</h3>",
        "<p><b>Ellison and Glaeser (1997).</b> With s<sub>r</sub> the region's share of the industry and x<sub>r</sub> its share of total employment, G = &sum;(s<sub>r</sub> - x<sub>r</sub>)<sup>2</sup> and &gamma; = (G - (1 - &sum;x<sub>r</sub><sup>2</sup>) H) / ((1 - &sum;x<sub>r</sub><sup>2</sup>)(1 - H)), where H is the Herfindahl index of plant sizes. With establishment counts H is computed as &sum;<sub>r</sub> E<sub>r</sub><sup>2</sup> / (n<sub>r</sub> E<sup>2</sup>), which treats plants within a region as equally large; plant-level data give the exact value.</p>",
        "<p><b>Locational Gini (Krugman, 1991).</b> Regions are ordered by location quotient; the Gini is twice the area between the Lorenz curve of the industry's cumulative share and the diagonal.</p>",
        "<p><b>Krugman specialisation index.</b> &sum;<sub>i</sub> |s<sub>ir</sub> - s<sub>i,rest</sub>| compares the industry shares of region r with those of all other regions.</p>",
        "<p>The indices do not account for the distances between regions or for the size of the regions' area; the Duranton-Overman distance-based measure needs plant coordinates and is not included.</p>"),
    gravity = paste0(
        "<h3>Gravity Model (PPML): method notes</h3>",
        "<p>The expected flow is modelled as E[y<sub>ij</sub>] = exp(x<sub>ij</sub>&beta; + &alpha;<sub>i</sub> + &delta;<sub>j</sub>). The pseudo-likelihood is the Poisson likelihood but only the conditional mean needs to be correct. It is maximised by iteratively reweighted least squares, absorbing the fixed effects by alternating projections (Correia, Guimaraes and Zylkin, 2020), so thousands of fixed effects are no problem.</p>",
        "<p><b>Standard errors.</b> Sandwich estimator with scores x&#771;<sub>ij</sub> (y<sub>ij</sub> - &mu;<sub>ij</sub>); clustering sums the scores within clusters (with a G/(G-1) correction) and two-way clustering adds the two one-way matrices and subtracts the intersection. A small-sample correction for the number of fixed effects is not applied.</p>",
        "<p><b>Checks.</b> Coefficients, the pseudo log-likelihood and clustered standard errors agree with fixest::fepois.</p>"),
    spdid = paste0(
        "<h3>Spatial DID (Spillovers): method notes</h3>",
        "<p>y<sub>it</sub> = &alpha;<sub>i</sub> + &gamma;<sub>t</sub> + &tau; D<sub>it</sub> + &sum;<sub>k</sub> &delta;<sub>k</sub> R<sup>k</sup><sub>it</sub> + &theta; (W D)<sub>it</sub> + x<sub>it</sub>&beta; + &epsilon;<sub>it</sub>. R<sup>k</sup><sub>it</sub> = 1 when unit i is untreated in period t and its nearest treated unit in t is in ring k. The coefficients are estimated by least squares after the two-way within transformation.</p>",
        "<p><b>Interpretation.</b> &tau; is the effect on treated units relative to units beyond the last ring; &delta;<sub>k</sub> is the spillover on untreated units in ring k relative to the same control group. If the true spillovers reach beyond the last ring, the control group is still contaminated and &tau; is biased; try larger rings.</p>",
        "<p><b>Standard errors.</b> Clustered: sandwich with the small-sample factor G/(G-1) (N-1)/(N-K) as in fixest and Stata. Conley: sum over periods s, t within the serial lag of &Psi;<sub>t</sub>' K &Psi;<sub>s</sub> with Bartlett weights in distance and time. The coefficients and clustered standard errors agree with fixest::feols; the Conley matrix agrees with its definition computed by explicit loops.</p>",
        "<p><b>Caution.</b> With staggered adoption or heterogeneous effects the two-way estimator mixes already-treated units into the comparison (Goodman-Bacon, 2021). Use it for a single adoption date or check the result with the robust estimators in the DiD module.</p>"),
    stop("unknown notes: ", key))
}

# Wrap the plain text in a styled block (kept simple so it renders the same in every theme).
sp_guide_html <- function(key) sp_guide_raw(key)
sp_notes_html <- function(key) sp_notes_raw(key)
