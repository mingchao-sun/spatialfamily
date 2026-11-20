# Third part of gen_yaml.py (executed from it): agglomeration and gravity.

# ============================================================================ agg
name, title = "agg", "Industrial Agglomeration"
AGG_CLEAR = ["region", "industry", "value", "firms"]
write(name, "a", a_head(name, title, "Agglomeration", "Industrial Agglomeration", "Location quotients, Gini, Ellison-Glaeser, Krugman") + """
    - name: region
      title: Region
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: industry
      title: Industry
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: value
      title: Employment (or output)
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: firms
      title: Number of establishments (optional)
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: showRegion
      title: Regional specialisation (Krugman index, Herfindahl)
      type: Bool
      default: true

    - name: sortBy
      title: Sort industries by
      type: List
      options:
          - name: gamma
            title: Ellison-Glaeser gamma
          - name: gini
            title: Locational Gini
          - name: total
            title: Size
          - name: name
            title: Name
      default: gamma

    - name: showGuide
      title: Data and setup guide
      type: Bool
      default: false

    - name: showHelp
      title: Method notes
      type: Bool
      default: false

    - name: saveLQ
      title: Location quotient
      type: Output
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: industryTable
      title: Concentration of the industries
      type: Table
      refs:
          - ellison1997
          - krugman1991
""" + clear(AGG_CLEAR + ["sortBy"]) + """      columns:
          - name: industry
            title: Industry
            type: text
          - name: total
            title: Total
            type: number
          - name: share
            title: Share of the total
            type: number
          - name: nreg
            title: Regions present
            type: integer
          - name: hhi
            title: Herfindahl (regional shares)
            type: number
          - name: G
            title: Raw concentration G
            type: number
          - name: H
            title: Plant Herfindahl H
            type: number
          - name: gamma
            title: EG gamma
            type: number
          - name: gini
            title: Locational Gini
            type: number

    - name: regionTable
      title: Specialisation of the regions
      type: Table
      visible: (showRegion)
      refs:
          - krugman1991
""" + clear(AGG_CLEAR + ["showRegion"]) + """      columns:
          - name: region
            title: Region
            type: text
          - name: total
            title: Total
            type: number
          - name: share
            title: Share of the total
            type: number
          - name: hhi
            title: Herfindahl (industry shares)
            type: number
          - name: entropy
            title: Entropy
            type: number
          - name: krugman
            title: Krugman index
            type: number

    - name: saveLQ
      title: Location quotient
      type: Output
      varTitle: Location quotient
      varDescription: Location quotient of the region-industry cell
      measureType: continuous
""" + clear(AGG_CLEAR))
write(name, "u", """name: agg
title: Industrial Agglomeration
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Region
            children:
                - type: VariablesListBox
                  name: region
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Industry
            children:
                - type: VariablesListBox
                  name: industry
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Employment (or output)
            children:
                - type: VariablesListBox
                  name: value
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Number of establishments (optional)
            children:
                - type: VariablesListBox
                  name: firms
                  isTarget: true
                  height: small
                  maxItemCount: 1
    - type: CollapseBox
      label: Results
      collapsed: false
      children:
          - type: CheckBox
            name: showRegion
          - type: ComboBox
            name: sortBy
    - type: CollapseBox
      label: Save to the dataset
      collapsed: true
      children:
          - type: Output
            name: saveLQ
""" + OUTPUT_U)

# ======================================================================== gravity
name, title = "gravity", "Gravity Model (PPML)"
GRAV_CLEAR = ["dep", "covs", "origin", "dest", "timeVar", "feOrigin", "feDest", "fePair", "feOrigTime", "feDestTime", "feTime", "cluster"]
write(name, "a", a_head(name, title, "Gravity", "Gravity Model (PPML)", "Poisson pseudo-ML with high-dimensional fixed effects") + """
    - name: dep
      title: Flow
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: covs
      title: Covariates (e.g. log distance, contiguity)
      type: Variables
      suggested:
          - continuous
          - nominal
      permitted:
          - numeric
          - factor

    - name: origin
      title: Origin
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: dest
      title: Destination
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: timeVar
      title: Period (optional)
      type: Variable
      permitted:
          - numeric
          - factor
          - id

    - name: feOrigin
      title: Origin fixed effects
      type: Bool
      default: true

    - name: feDest
      title: Destination fixed effects
      type: Bool
      default: true

    - name: fePair
      title: Origin-destination pair fixed effects
      type: Bool
      default: false

    - name: feOrigTime
      title: Origin x period fixed effects
      type: Bool
      default: false

    - name: feDestTime
      title: Destination x period fixed effects
      type: Bool
      default: false

    - name: feTime
      title: Period fixed effects
      type: Bool
      default: false

    - name: cluster
      title: Standard errors
      type: List
      options:
          - name: robust
            title: Heteroskedasticity-robust
          - name: origin
            title: Clustered by origin
          - name: dest
            title: Clustered by destination
          - name: pair
            title: Clustered by origin-destination pair
          - name: twoway
            title: Two-way clustered (origin and destination)
      default: robust

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9

    - name: showGuide
      title: Data and setup guide
      type: Bool
      default: false

    - name: showHelp
      title: Method notes
      type: Bool
      default: false
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: coef
      title: PPML coefficients
      type: Table
      refs:
          - santossilva2006
          - correia2020
""" + clear(GRAV_CLEAR + ["level"]) + """      columns:
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
            title: z
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
          - name: pct
            title: '% change in flow'
            type: number

    - name: fit
      title: Model fit
      type: Table
""" + clear(GRAV_CLEAR) + """      columns:
          - name: item
            title: ''
            type: text
          - name: value
            title: ''
            type: number
""")
write(name, "u", """name: gravity
title: Gravity Model (PPML)
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Flow
            children:
                - type: VariablesListBox
                  name: dep
                  isTarget: true
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Covariates (e.g. log distance, contiguity)
            children:
                - type: VariablesListBox
                  name: covs
                  isTarget: true
          - type: TargetLayoutBox
            label: Origin
            children:
                - type: VariablesListBox
                  name: origin
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Destination
            children:
                - type: VariablesListBox
                  name: dest
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Period (optional)
            children:
                - type: VariablesListBox
                  name: timeVar
                  isTarget: true
                  height: small
                  maxItemCount: 1
    - type: CollapseBox
      label: Fixed effects
      collapsed: false
      children:
          - type: CheckBox
            name: feOrigin
          - type: CheckBox
            name: feDest
          - type: CheckBox
            name: fePair
          - type: CheckBox
            name: feOrigTime
          - type: CheckBox
            name: feDestTime
          - type: CheckBox
            name: feTime
    - type: CollapseBox
      label: Inference
      collapsed: false
      children:
          - type: ComboBox
            name: cluster
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
""" + OUTPUT_U)
