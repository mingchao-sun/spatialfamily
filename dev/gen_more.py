# Second part of gen_yaml.py (executed from it): geographically weighted regression, agglomeration, gravity, spatial DID.

# ============================================================================ gwr
name, title = "gwr", "Geographically Weighted Regression"
GWR_CLEAR = ["dep", "covs", "xVar", "yVar", "wLonLat", "kernel", "bwType", "bwMethod", "bw"]
write(name, "a", a_head(name, title, "Regression", "Geographically Weighted Regression", "Local coefficients that vary over space (GWR)") + """
    - name: dep
      title: Dependent variable
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: covs
      title: Covariates
      type: Variables
      suggested:
          - continuous
      permitted:
          - numeric

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

    - name: idVar
      title: Unit id (optional, labels the local table)
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: wLonLat
      title: Coordinates are longitude / latitude (distances in km)
      type: Bool
      default: false

    - name: kernel
      title: Kernel
      type: List
      options:
          - name: bisquare
            title: Bi-square
          - name: gaussian
            title: Gaussian
          - name: exponential
            title: Exponential
          - name: tricube
            title: Tricube
          - name: boxcar
            title: Box-car
      default: bisquare

    - name: bwType
      title: Bandwidth type
      type: List
      options:
          - name: adaptive
            title: Adaptive (number of nearest neighbours)
          - name: fixed
            title: Fixed (distance)
      default: adaptive

    - name: bwMethod
      title: Bandwidth
      type: List
      options:
          - name: aicc
            title: Select by AICc
          - name: cv
            title: Select by cross-validation
          - name: manual
            title: Set manually
      default: aicc

    - name: bw
      title: Manual bandwidth
      type: Number
      default: 30
      min: 0

    - name: showLocal
      title: Table of local estimates
      type: Bool
      default: false

    - name: maps
      title: Maps of the local coefficients
      type: Bool
      default: true

    - name: mapStat
      title: Map shows
      type: List
      options:
          - name: coef
            title: Local coefficient
          - name: t
            title: Local t value
      default: coef

    - name: showGuide
      title: Data and setup guide
      type: Bool
      default: false

    - name: showHelp
      title: Method notes
      type: Bool
      default: false

    - name: saveR2
      title: Local R-squared
      type: Output

    - name: saveResid
      title: Residuals
      type: Output
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: summary
      title: 'Coefficients: global (OLS) and geographically weighted'
      type: Table
      refs:
          - brunsdon1996
          - gwmodel
""" + clear(GWR_CLEAR) + """      columns:
          - name: term
            title: ''
            type: text
          - name: ols
            title: OLS
            type: number
          - name: olsse
            title: OLS s.e.
            type: number
          - name: min
            title: Minimum
            type: number
          - name: q1
            title: 1st quartile
            type: number
          - name: med
            title: Median
            type: number
          - name: q3
            title: 3rd quartile
            type: number
          - name: max
            title: Maximum
            type: number
          - name: nsig
            title: '% units |t| > 1.96'
            type: number

    - name: diag
      title: Diagnostics
      type: Table
""" + clear(GWR_CLEAR) + """      columns:
          - name: item
            title: ''
            type: text
          - name: gwr
            title: GWR
            type: number
          - name: ols
            title: OLS
            type: number

    - name: local
      title: Local estimates
      type: Table
      visible: (showLocal)
""" + clear(GWR_CLEAR + ["idVar", "showLocal"]) + """      columns:
          - name: id
            title: Unit
            type: text
          - name: x
            title: X
            type: number
          - name: y
            title: Y
            type: number
          - name: fitted
            title: Fitted
            type: number
          - name: r2
            title: Local R2
            type: number

    - name: mapPlots
      title: Map
      type: Array
      items: (covs)
      visible: (maps)
      template:
          title: $key
          type: Image
          width: 520
          height: 440
          renderFun: .plotMap
""" + clear(GWR_CLEAR + ["maps", "mapStat"]) + """
    - name: mapInt
      title: Intercept
      type: Image
      width: 520
      height: 440
      renderFun: .plotMapInt
      visible: (maps)
""" + clear(GWR_CLEAR + ["maps", "mapStat"]) + """
    - name: saveR2
      title: Local R-squared
      type: Output
      varTitle: GWR local R2
      varDescription: Local R-squared of the geographically weighted regression
      measureType: continuous
""" + clear(GWR_CLEAR) + """
    - name: saveResid
      title: Residuals
      type: Output
      varTitle: GWR residuals
      varDescription: Residuals of the geographically weighted regression
      measureType: continuous
""" + clear(GWR_CLEAR))
write(name, "u", """name: gwr
title: Geographically Weighted Regression
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
            label: Covariates
            children:
                - type: VariablesListBox
                  name: covs
                  isTarget: true
          - type: TargetLayoutBox
            label: X / longitude
            children:
                - type: VariablesListBox
                  name: xVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Y / latitude
            children:
                - type: VariablesListBox
                  name: yVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Unit id (optional, labels the local table)
            children:
                - type: VariablesListBox
                  name: idVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
    - type: CheckBox
      name: wLonLat
    - type: CollapseBox
      label: Kernel and bandwidth
      collapsed: false
      children:
          - type: ComboBox
            name: kernel
          - type: ComboBox
            name: bwType
          - type: ComboBox
            name: bwMethod
          - type: TextBox
            name: bw
            label: Manual bandwidth (neighbours, or distance)
            format: number
            enable: (bwMethod:manual)
    - type: CollapseBox
      label: Results
      collapsed: false
      children:
          - type: CheckBox
            name: maps
          - type: ComboBox
            name: mapStat
          - type: CheckBox
            name: showLocal
    - type: CollapseBox
      label: Save to the dataset
      collapsed: true
      children:
          - type: Output
            name: saveR2
          - type: Output
            name: saveResid
""" + OUTPUT_U)
