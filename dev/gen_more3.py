# Fourth part of gen_yaml.py (executed from it): spatial difference-in-differences.

name, title = "spdid", "Spatial DID (Spillovers)"
SPDID_CLEAR = ["dep", "covs", "treat", "timeVar", "rings", "lagTreat", "effects", "se", "conleyBw", "conleyLag"]
write(name, "a", a_head(name, title, "Causal", "Spatial DID (Spillovers)", "TWFE with distance rings and spillover to neighbours") + """
    - name: dep
      title: Outcome
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: treat
      title: Treatment (1 when the unit is treated in the period)
      type: Variable
      permitted:
          - numeric
          - factor

    - name: covs
      title: Covariates (optional)
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
    - name: rings
      title: Distance rings (upper limits, e.g. 5, 10, 20)
      type: String
      default: ''

    - name: lagTreat
      title: Add the spatially lagged treatment (W x treatment)
      type: Bool
      default: false

    - name: effects
      title: Fixed effects
      type: List
      options:
          - name: twoway
            title: Unit and time
          - name: unit
            title: Unit only
          - name: time
            title: Time only
      default: twoway

    - name: se
      title: Standard errors
      type: List
      options:
          - name: cluster
            title: Clustered by unit
          - name: robust
            title: Heteroskedasticity-robust
          - name: conley
            title: Conley spatial HAC
      default: cluster

    - name: conleyBw
      title: Conley distance bandwidth
      type: Number
      default: 0
      min: 0

    - name: conleyLag
      title: Conley serial correlation lag (periods)
      type: Integer
      default: 0
      min: 0
      max: 20

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9
""")
write(name, "r", R_HEAD.format(name=name, title=title) + """
    - name: coef
      title: Treatment and spillover effects
      type: Table
      refs:
          - conley1999
          - butts2021
""" + clear(SPDID_CLEAR + ["level"]) + """      columns:
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
            title: t
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
      title: Sample and fit
      type: Table
""" + clear(SPDID_CLEAR) + """      columns:
          - name: item
            title: ''
            type: text
          - name: value
            title: ''
            type: number
""")
write(name, "u", """name: spdid
title: Spatial DID (Spillovers)
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
          - type: TargetLayoutBox
            label: Outcome
            children:
                - type: VariablesListBox
                  name: dep
                  isTarget: true
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Treatment (1 when the unit is treated in the period)
            children:
                - type: VariablesListBox
                  name: treat
                  isTarget: true
                  height: small
                  maxItemCount: 1
          - type: TargetLayoutBox
            label: Covariates (optional)
            children:
                - type: VariablesListBox
                  name: covs
                  isTarget: true
                  height: small
          - type: TargetLayoutBox
            label: Time period
            children:
                - type: VariablesListBox
                  name: timeVar
                  isTarget: true
                  height: small
                  maxItemCount: 1""" + LOC_U.replace("Unit id (needed for files, optional for coordinates)", "Unit id") + W_U + """    - type: CollapseBox
      label: Spillovers
      collapsed: false
      children:
          - type: TextBox
            name: rings
            label: Distance rings (upper limits, e.g. 5, 10, 20)
            format: string
            stretchFactor: 1
          - type: CheckBox
            name: lagTreat
    - type: CollapseBox
      label: Estimation
      collapsed: false
      children:
          - type: ComboBox
            name: effects
          - type: ComboBox
            name: se
          - type: TextBox
            name: conleyBw
            label: Conley distance bandwidth
            format: number
            enable: (se:conley)
          - type: TextBox
            name: conleyLag
            label: Conley serial correlation lag
            format: number
            enable: (se:conley)
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
""" + OUTPUT_U)
