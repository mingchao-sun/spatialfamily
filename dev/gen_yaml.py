"""Generates jamovi/*.a.yaml, *.r.yaml and *.u.yaml of the four analyses (the spatial weights options are shared).
Run from the module folder:  python dev/gen_yaml.py
"""
import os, textwrap

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "jamovi")

# ------------------------------------------------------------------ shared options
LOC_A = """
    - name: idVar
      title: Unit id
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: xVar
      title: X / longitude
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: yVar
      title: Y / latitude
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric
"""

W_A = """
    - name: wSource
      title: Spatial weights from
      type: List
      options:
          - name: coords
            title: Coordinates (x, y columns)
          - name: boundary
            title: Boundary file (shapefile or GeoJSON)
          - name: adjacency
            title: Neighbour list file (edge list, GAL, GWT)
      default: coords

    - name: wFile
      title: File path
      type: String
      default: ''

    - name: wIdField
      title: Id field in the boundary file
      type: String
      default: ''

    - name: wMethod
      title: Neighbours defined by
      type: List
      options:
          - name: auto
            title: Default for the source (k nearest / queen)
          - name: knn
            title: k nearest neighbours
          - name: dist
            title: Distance band
          - name: idw
            title: Inverse distance
          - name: queen
            title: Queen contiguity (boundary file)
          - name: rook
            title: Rook contiguity (boundary file)
      default: auto

    - name: wK
      title: Number of neighbours (k)
      type: Integer
      default: 5
      min: 1
      max: 500

    - name: wDist
      title: Distance threshold (0 = automatic)
      type: Number
      default: 0
      min: 0

    - name: wPower
      title: Distance decay exponent
      type: Number
      default: 1
      min: 0

    - name: wCutoff
      title: Distance cutoff for inverse distance (0 = none)
      type: Number
      default: 0
      min: 0

    - name: wOrder
      title: Order of neighbours
      type: Integer
      default: 1
      min: 1
      max: 10

    - name: wCumul
      title: Include lower orders
      type: Bool
      default: false

    - name: wSymm
      title: Make links symmetric
      type: Bool
      default: false

    - name: wLonLat
      title: Coordinates are longitude / latitude (distances in km)
      type: Bool
      default: false

    - name: wStyle
      title: Weights style
      type: List
      options:
          - name: W
            title: Row-standardised
          - name: B
            title: As built (binary or distance)
          - name: C
            title: Globally standardised
      default: W

    - name: showGuide
      title: Data and setup guide
      type: Bool
      default: false

    - name: showHelp
      title: Method notes
      type: Bool
      default: false
"""

W_CLEAR = ["idVar", "xVar", "yVar", "wSource", "wFile", "wIdField", "wMethod", "wK", "wDist", "wPower", "wCutoff", "wOrder", "wCumul", "wSymm", "wLonLat", "wStyle"]

LOC_U = """
          - type: TargetLayoutBox
            label: Unit id (needed for files, optional for coordinates)
            children:
                - type: VariablesListBox
                  name: idVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: X / longitude (coordinates)
            children:
                - type: VariablesListBox
                  name: xVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Y / latitude (coordinates)
            children:
                - type: VariablesListBox
                  name: yVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
"""

W_U = """
    - type: CollapseBox
      label: Spatial weights
      collapsed: false
      children:
          - type: ComboBox
            name: wSource
          - type: TextBox
            name: wFile
            label: File path
            format: string
            stretchFactor: 1
            enable: (wSource:boundary || wSource:adjacency)
          - type: Label
            label: 'Tip: Shift + right-click the file in File Explorer,'
            enable: (wSource:boundary || wSource:adjacency)
          - type: Label
            label: 'choose Copy as path, then paste it here (or sample:nc.shp).'
            enable: (wSource:boundary || wSource:adjacency)
          - type: TextBox
            name: wIdField
            label: Id field in the boundary file
            format: string
            stretchFactor: 1
            enable: (wSource:boundary)
          - type: ComboBox
            name: wMethod
          - type: TextBox
            name: wK
            label: Number of neighbours (k)
            format: number
          - type: TextBox
            name: wDist
            label: Distance threshold (0 = automatic)
            format: number
          - type: TextBox
            name: wPower
            label: Distance decay exponent
            format: number
          - type: TextBox
            name: wCutoff
            label: Distance cutoff (0 = none)
            format: number
          - type: TextBox
            name: wOrder
            label: Order of neighbours
            format: number
          - type: CheckBox
            name: wCumul
          - type: CheckBox
            name: wSymm
          - type: CheckBox
            name: wLonLat
          - type: ComboBox
            name: wStyle
"""

OUTPUT_U = """
    - type: CollapseBox
      label: Output
      collapsed: true
      children:
          - type: CheckBox
            name: showGuide
          - type: CheckBox
            name: showHelp
"""

R_HEAD = """---
name: {name}
title: {title}
jrs: '1.1'

items:
    - name: guide
      title: Data and setup guide
      type: Html
      visible: (showGuide)

    - name: helpText
      title: Method notes
      type: Html
      visible: (showHelp)

    - name: notes
      title: Notes
      type: Html
      visible: true
"""


def clear(extra):
    items = W_CLEAR + extra
    return "      clearWith:\n" + "".join("          - %s\n" % i for i in items)


def write(name, ext, text):
    with open(os.path.join(OUT, "%s.%s.yaml" % (name, ext)), "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def a_head(name, title, sub, mt, ms):
    return ("---\nname: %s\ntitle: %s\njas: '1.2'\nversion: '1.0.0'\nmenuGroup: Spatial\nmenuSubgroup: %s\nmenuTitle: %s\nmenuSubtitle: %s\n\noptions:\n    - name: data\n      type: Data\n" %
            (name, title, sub, mt, ms))


# ============================================================================ swm
name, title = "swm", "Spatial Weights"
write(name, "a", a_head(name, title, "Weights", "Spatial Weights", "Neighbours and weights matrix") + LOC_A + W_A + """
    - name: linkPlot
      title: Connectivity map
      type: Bool
      default: true

    - name: histPlot
      title: Histogram of the number of neighbours
      type: Bool
      default: false
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: summary
      title: Spatial weights
      type: Table
""" + clear([]) + """      columns:
          - name: item
            title: ''
            type: text
          - name: value
            title: ''
            type: text

    - name: nbTable
      title: Distribution of the number of neighbours
      type: Table
""" + clear([]) + """      columns:
          - name: k
            title: Neighbours
            type: integer
          - name: n
            title: Units
            type: integer
          - name: pct
            title: '%'
            type: number

    - name: linkMap
      title: Connectivity map
      type: Image
      width: 560
      height: 480
      renderFun: .plotLinks
      visible: (linkPlot)
""" + clear([]) + """
    - name: nbHist
      title: Number of neighbours
      type: Image
      width: 480
      height: 360
      renderFun: .plotHist
      visible: (histPlot)
""" + clear([]))
write(name, "u", """name: swm
title: Spatial Weights
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:""" + LOC_U + W_U + """    - type: CheckBox
      name: linkPlot
    - type: CheckBox
      name: histPlot
""" + OUTPUT_U)

# ========================================================================== moran
name, title = "moran", "Global Spatial Autocorrelation"
write(name, "a", a_head(name, title, "Autocorrelation", "Global (Moran, Geary, G)", "Is the variable spatially clustered?") + """
    - name: vars
      title: Variables
      type: Variables
      suggested:
          - continuous
      permitted:
          - numeric
""" + LOC_A + W_A + """
    - name: inference
      title: Inference
      type: List
      options:
          - name: randomisation
            title: Randomisation
          - name: normality
            title: Normality
          - name: permutation
            title: Permutation
      default: randomisation

    - name: alt
      title: Alternative hypothesis
      type: List
      options:
          - name: two.sided
            title: Two-sided
          - name: greater
            title: Positive autocorrelation
          - name: less
            title: Negative autocorrelation
      default: two.sided

    - name: nperm
      title: Number of permutations
      type: Integer
      default: 999
      min: 99
      max: 99999

    - name: seed
      title: Random seed
      type: Integer
      default: 1
      min: 1

    - name: geary
      title: Geary's C
      type: Bool
      default: false

    - name: getisG
      title: Getis-Ord general G
      type: Bool
      default: false

    - name: scatter
      title: Moran scatterplot
      type: Bool
      default: false
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: global
      title: Global spatial autocorrelation
      type: Table
      refs:
          - cliff1981
          - moran1950
          - spdep
""" + clear(["vars", "inference", "alt", "nperm", "seed", "geary", "getisG"]) + """      columns:
          - name: var
            title: Variable
            type: text
          - name: stat
            title: Statistic
            type: text
          - name: est
            title: Estimate
            type: number
          - name: exp
            title: Expected
            type: number
          - name: sd
            title: Std. deviate
            type: number
          - name: z
            title: z
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue

    - name: scatterPlots
      title: Moran scatterplot
      type: Array
      items: (vars)
      visible: (scatter)
      template:
          title: $key
          type: Image
          width: 480
          height: 420
          renderFun: .plotScatter
""" + clear(["vars", "scatter"]))
write(name, "u", """name: moran
title: Global Spatial Autocorrelation
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Variables
            children:
                - type: VariablesListBox
                  name: vars
                  isTarget: true""" + LOC_U + W_U + """    - type: CollapseBox
      label: Statistics
      collapsed: false
      children:
          - type: CheckBox
            name: geary
          - type: CheckBox
            name: getisG
          - type: CheckBox
            name: scatter
    - type: CollapseBox
      label: Inference
      collapsed: false
      children:
          - type: ComboBox
            name: inference
          - type: ComboBox
            name: alt
          - type: TextBox
            name: nperm
            label: Number of permutations
            format: number
          - type: TextBox
            name: seed
            label: Random seed
            format: number
""" + OUTPUT_U)

# =========================================================================== lisa
name, title = "lisa", "Local Spatial Autocorrelation"
write(name, "a", a_head(name, title, "Autocorrelation", "Local (LISA, Getis-Ord Gi*)", "Where are the clusters and hot spots?") + """
    - name: var
      title: Variable
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric
""" + LOC_A + W_A + """
    - name: stat
      title: Statistic
      type: List
      options:
          - name: localmoran
            title: Local Moran's I (LISA)
          - name: gistar
            title: Getis-Ord Gi* (own value included)
          - name: gi
            title: Getis-Ord Gi (own value excluded)
      default: localmoran

    - name: nperm
      title: Number of permutations
      type: Integer
      default: 999
      min: 99
      max: 99999

    - name: seed
      title: Random seed
      type: Integer
      default: 1
      min: 1

    - name: alpha
      title: Significance level
      type: Number
      default: 0.05
      min: 0.0001
      max: 0.5

    - name: adjust
      title: Multiple testing
      type: List
      options:
          - name: none
            title: None
          - name: fdr
            title: False discovery rate (Benjamini-Hochberg)
          - name: bonferroni
            title: Bonferroni
      default: none

    - name: showUnits
      title: Table of all units
      type: Bool
      default: false

    - name: map
      title: Cluster map
      type: Bool
      default: true

    - name: scatter
      title: Moran scatterplot
      type: Bool
      default: false

    - name: saveCluster
      title: Cluster class
      type: Output

    - name: saveStat
      title: Local statistic
      type: Output

    - name: saveP
      title: Pseudo p-value
      type: Output
""")
LISA_CLEAR = ["var", "stat", "nperm", "seed", "alpha", "adjust"]
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: counts
      title: Clusters
      type: Table
      refs:
          - anselin1995
          - ord1995
""" + clear(LISA_CLEAR) + """      columns:
          - name: cls
            title: Class
            type: text
          - name: n
            title: Units
            type: integer
          - name: pct
            title: '%'
            type: number

    - name: units
      title: Units
      type: Table
      visible: (showUnits)
""" + clear(LISA_CLEAR + ["showUnits"]) + """      columns:
          - name: id
            title: Unit
            type: text
          - name: x
            title: Value
            type: number
          - name: lag
            title: Neighbour mean
            type: number
          - name: stat
            title: Statistic
            type: number
          - name: z
            title: z
            type: number
          - name: p
            title: p (permutation)
            type: number
            format: zto,pvalue
          - name: padj
            title: p (adjusted)
            type: number
            format: zto,pvalue
          - name: cls
            title: Class
            type: text

    - name: clusterMap
      title: Cluster map
      type: Image
      width: 600
      height: 500
      renderFun: .plotMap
      visible: (map)
""" + clear(LISA_CLEAR + ["map"]) + """
    - name: scatterPlot
      title: Moran scatterplot
      type: Image
      width: 480
      height: 420
      renderFun: .plotScatter
      visible: (scatter)
""" + clear(LISA_CLEAR + ["scatter"]) + """
    - name: saveCluster
      title: Cluster class
      type: Output
      varTitle: LISA cluster
      varDescription: Local cluster class (0 not significant, 1 High-High, 2 Low-Low, 3 Low-High, 4 High-Low)
      measureType: nominal
""" + clear(LISA_CLEAR) + """
    - name: saveStat
      title: Local statistic
      type: Output
      varTitle: Local statistic
      varDescription: Local Moran's I or Getis-Ord z-score
      measureType: continuous
""" + clear(LISA_CLEAR) + """
    - name: saveP
      title: Pseudo p-value
      type: Output
      varTitle: Local p-value
      varDescription: Permutation pseudo p-value of the local statistic
      measureType: continuous
""" + clear(LISA_CLEAR))
write(name, "u", """name: lisa
title: Local Spatial Autocorrelation
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Variable
            children:
                - type: VariablesListBox
                  name: var
                  isTarget: true
                  maxItemCount: 1""" + LOC_U + W_U + """    - type: CollapseBox
      label: Statistic and inference
      collapsed: false
      children:
          - type: ComboBox
            name: stat
          - type: TextBox
            name: nperm
            label: Number of permutations
            format: number
          - type: TextBox
            name: seed
            label: Random seed
            format: number
          - type: TextBox
            name: alpha
            label: Significance level
            format: number
          - type: ComboBox
            name: adjust
    - type: CollapseBox
      label: Results
      collapsed: false
      children:
          - type: CheckBox
            name: map
          - type: CheckBox
            name: scatter
          - type: CheckBox
            name: showUnits
    - type: CollapseBox
      label: Save to the dataset
      collapsed: true
      children:
          - type: Output
            name: saveCluster
          - type: Output
            name: saveStat
          - type: Output
            name: saveP
""" + OUTPUT_U)

# =========================================================================== sreg
name, title = "sreg", "Spatial Regression"
SREG_CLEAR = ["dep", "covs", "noLag", "model"]
write(name, "a", a_head(name, title, "Regression", "Spatial Regression", "SAR, SEM, SDM, SDEM, SAC, SLX (maximum likelihood)") + """
    - name: dep
      title: Dependent variable
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: covs
      title: Covariates (lagged in SLX / SDM / SDEM)
      type: Variables
      suggested:
          - continuous
          - nominal
      permitted:
          - numeric
          - factor

    - name: noLag
      title: Covariates not to be lagged
      type: Variables
      suggested:
          - continuous
          - nominal
      permitted:
          - numeric
          - factor
""" + LOC_A + W_A + """
    - name: model
      title: Model
      type: List
      options:
          - name: ols
            title: OLS (no spatial terms)
          - name: slx
            title: Spatially lagged covariates (SLX)
          - name: sar
            title: Spatial lag (SAR)
          - name: sem
            title: Spatial error (SEM)
          - name: sdm
            title: Spatial Durbin (SDM)
          - name: sdem
            title: Spatial Durbin error (SDEM)
          - name: sac
            title: Spatial lag + error (SAC)
      default: sar

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9

    - name: diag
      title: Diagnostics for OLS (LM tests)
      type: Bool
      default: true

    - name: effects
      title: Direct, indirect and total effects
      type: Bool
      default: true

    - name: nsim
      title: Simulations for the effects
      type: Integer
      default: 500
      min: 0
      max: 10000

    - name: seed
      title: Random seed
      type: Integer
      default: 1
      min: 1

    - name: moranRes
      title: Moran's I of the residuals
      type: Bool
      default: true

    - name: compare
      title: Compare all models
      type: Bool
      default: false

    - name: saveFitted
      title: Fitted values
      type: Output

    - name: saveResid
      title: Residuals
      type: Output
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: coef
      title: Coefficients
      type: Table
      refs:
          - anselin1988
          - ord1975
          - spatialreg
""" + clear(SREG_CLEAR + ["level"]) + """      columns:
          - name: term
            title: ''
            type: text
          - name: est
            title: Estimate
            type: number
          - name: se
            title: Std. error
            type: number
          - name: z
            title: z / t
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue
          - name: lo
            title: Lower
            type: number
          - name: hi
            title: Upper
            type: number

    - name: fit
      title: Model fit
      type: Table
""" + clear(SREG_CLEAR) + """      columns:
          - name: item
            title: ''
            type: text
          - name: value
            title: ''
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue

    - name: lmTests
      title: Diagnostics for spatial dependence in the OLS residuals
      type: Table
      visible: (diag)
      refs:
          - anselin1996
""" + clear(["dep", "covs", "noLag", "diag"]) + """      columns:
          - name: test
            title: Test
            type: text
          - name: stat
            title: Statistic
            type: number
          - name: df
            title: df
            type: integer
          - name: p
            title: p
            type: number
            format: zto,pvalue

    - name: effectsTable
      title: Average effects of the covariates
      type: Table
      visible: (effects)
      refs:
          - lesage2009
""" + clear(SREG_CLEAR + ["effects", "nsim", "seed"]) + """      columns:
          - name: term
            title: Covariate
            type: text
          - name: kind
            title: Effect
            type: text
          - name: est
            title: Estimate
            type: number
          - name: se
            title: Std. error
            type: number
          - name: z
            title: z
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue

    - name: moranTable
      title: Moran's I of the residuals
      type: Table
      rows: 1
      visible: (moranRes)
""" + clear(SREG_CLEAR + ["moranRes"]) + """      columns:
          - name: stat
            title: Statistic
            type: number
          - name: exp
            title: Expected
            type: number
          - name: z
            title: z
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue

    - name: compareTable
      title: Comparison of models
      type: Table
      visible: (compare)
""" + clear(["dep", "covs", "noLag", "compare"]) + """      columns:
          - name: model
            title: Model
            type: text
          - name: rho
            title: rho
            type: number
          - name: lambda
            title: lambda
            type: number
          - name: ll
            title: Log-likelihood
            type: number
          - name: k
            title: Parameters
            type: integer
          - name: aic
            title: AIC
            type: number
          - name: bic
            title: BIC
            type: number

    - name: saveFitted
      title: Fitted values
      type: Output
      varTitle: Spatial fitted
      varDescription: Fitted values of the spatial regression (including the spatial lag of the dependent variable)
      measureType: continuous
""" + clear(SREG_CLEAR) + """
    - name: saveResid
      title: Residuals
      type: Output
      varTitle: Spatial residuals
      varDescription: Residuals of the spatial regression
      measureType: continuous
""" + clear(SREG_CLEAR))
write(name, "u", """name: sreg
title: Spatial Regression
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Dependent variable
            children:
                - type: VariablesListBox
                  name: dep
                  isTarget: true
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Covariates (lagged in SLX / SDM / SDEM)
            children:
                - type: VariablesListBox
                  name: covs
                  isTarget: true
          - type: TargetLayoutBox
            label: Covariates not to be lagged
            children:
                - type: VariablesListBox
                  name: noLag
                  isTarget: true
                  height: small""" + LOC_U + W_U + """    - type: CollapseBox
      label: Model
      collapsed: false
      children:
          - type: ComboBox
            name: model
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
          - type: CheckBox
            name: diag
          - type: CheckBox
            name: moranRes
          - type: CheckBox
            name: compare
    - type: CollapseBox
      label: Effects
      collapsed: true
      children:
          - type: CheckBox
            name: effects
          - type: TextBox
            name: nsim
            label: Simulations
            format: number
          - type: TextBox
            name: seed
            label: Random seed
            format: number
    - type: CollapseBox
      label: Save to the dataset
      collapsed: true
      children:
          - type: Output
            name: saveFitted
          - type: Output
            name: saveResid
""" + OUTPUT_U)

# =========================================================================== span
name, title = "span", "Spatial Panel Regression"
SPAN_CLEAR = ["dep", "covs", "noLag", "timeVar", "model", "panelEffects"]
write(name, "a", a_head(name, title, "Regression", "Spatial Panel Regression", "Fixed effects SAR, SEM, SDM, SDEM, SLX (maximum likelihood)") + """
    - name: dep
      title: Dependent variable
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: covs
      title: Covariates (lagged in SLX / SDM / SDEM)
      type: Variables
      suggested:
          - continuous
      permitted:
          - numeric

    - name: noLag
      title: Covariates not to be lagged
      type: Variables
      suggested:
          - continuous
      permitted:
          - numeric

    - name: timeVar
      title: Time period
      type: Variable
      permitted:
          - numeric
          - factor
          - id
""" + LOC_A + W_A + """
    - name: model
      title: Model
      type: List
      options:
          - name: ols
            title: Non-spatial panel (OLS)
          - name: slx
            title: Spatially lagged covariates (SLX)
          - name: sar
            title: Spatial lag (SAR)
          - name: sem
            title: Spatial error (SEM)
          - name: sdm
            title: Spatial Durbin (SDM)
          - name: sdem
            title: Spatial Durbin error (SDEM)
      default: sar

    - name: panelEffects
      title: Fixed effects
      type: List
      options:
          - name: individual
            title: Unit fixed effects
          - name: time
            title: Time fixed effects
          - name: twoway
            title: Unit and time fixed effects
          - name: none
            title: None (pooled, with a constant)
      default: individual

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9

    - name: showEffects
      title: Direct, indirect and total effects
      type: Bool
      default: true

    - name: nsim
      title: Simulations for the effects
      type: Integer
      default: 500
      min: 0
      max: 10000

    - name: seed
      title: Random seed
      type: Integer
      default: 1
      min: 1
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: coef
      title: Coefficients
      type: Table
      refs:
          - elhorst2014
          - lesage2009
""" + clear(SPAN_CLEAR + ["level"]) + """      columns:
          - name: term
            title: ''
            type: text
          - name: est
            title: Estimate
            type: number
          - name: se
            title: Std. error
            type: number
          - name: z
            title: z / t
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue
          - name: lo
            title: Lower
            type: number
          - name: hi
            title: Upper
            type: number

    - name: fit
      title: Model fit
      type: Table
""" + clear(SPAN_CLEAR) + """      columns:
          - name: item
            title: ''
            type: text
          - name: value
            title: ''
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue

    - name: effectsTable
      title: Average effects of the covariates
      type: Table
      visible: (showEffects)
""" + clear(SPAN_CLEAR + ["showEffects", "nsim", "seed"]) + """      columns:
          - name: term
            title: Covariate
            type: text
          - name: kind
            title: Effect
            type: text
          - name: est
            title: Estimate
            type: number
          - name: se
            title: Std. error
            type: number
          - name: z
            title: z
            type: number
          - name: p
            title: p
            type: number
            format: zto,pvalue
""")
write(name, "u", """name: span
title: Spatial Panel Regression
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Dependent variable
            children:
                - type: VariablesListBox
                  name: dep
                  isTarget: true
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Covariates (lagged in SLX / SDM / SDEM)
            children:
                - type: VariablesListBox
                  name: covs
                  isTarget: true
          - type: TargetLayoutBox
            label: Covariates not to be lagged
            children:
                - type: VariablesListBox
                  name: noLag
                  isTarget: true
                  height: small
          - type: TargetLayoutBox
            label: Time period
            children:
                - type: VariablesListBox
                  name: timeVar
                  isTarget: true
                  height: small
                  maxItemCount: 1""" + LOC_U.replace("Unit id (needed for files, optional for coordinates)", "Unit id (one value per unit)") + W_U + """    - type: CollapseBox
      label: Model
      collapsed: false
      children:
          - type: ComboBox
            name: model
          - type: ComboBox
            name: panelEffects
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
    - type: CollapseBox
      label: Effects
      collapsed: true
      children:
          - type: CheckBox
            name: showEffects
          - type: TextBox
            name: nsim
            label: Simulations
            format: number
          - type: TextBox
            name: seed
            label: Random seed
            format: number
""" + OUTPUT_U)

exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_more.py"), encoding="utf-8").read())
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_more2.py"), encoding="utf-8").read())
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_more3.py"), encoding="utf-8").read())
print("written", sorted(os.listdir(OUT)))
