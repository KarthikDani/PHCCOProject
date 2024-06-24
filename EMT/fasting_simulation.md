### Load both the models and compare the growth rates


```python
import cobra
from cobra.flux_analysis import flux_variability_analysis
```


```python
# Load the original and modified models
original_model = cobra.io.read_sbml_model("""/Users/karthik/Desktop/
                                          PHCCO IISc Internship/EMT/
                                          MODEL1602080000_url.xml""")

integrated_model = cobra.io.read_sbml_model("""/Users/karthik/Desktop/
                                            PHCCO IISc Internship/EMT/
                                            integrated_MODEL1602080000.xml""")
```


```python
# Compare growth rates
original_solution = original_model.optimize()
modified_solution = integrated_model.optimize()

print("Original model growth rate:", original_solution.objective_value)
print("Modified model growth rate:", modified_solution.objective_value)
```

    Original model growth rate: 0.10099999998857717
    Modified model growth rate: 0.10099999998857717


### Perform FVA
- To understand the flexibility of each metabolic reaction by determining the maximum and minimum possible flux values (flux range) it can carry under steady-state conditions and given constraints.
- Before and after adjusting reaction bounds using functions like adjust_bounds, FVA can validate whether the adjusted constraints allow the reactions to operate within feasible physiological bounds.
- Identifies reactions that are critical (essential) or flexible (non-essential), which can further help in prioritising the reactions.


```python
# Perform FVA on both models
original_fva = flux_variability_analysis(original_model)
modified_fva = flux_variability_analysis(integrated_model)

print("FVA results for original model:")
print(original_fva)
print("FVA results for modified model:")
print(modified_fva)
```

    FVA results for original model:
                                     minimum       maximum
    10FTHF7GLUtl                0.000000e+00   1666.666667
    10FTHF7GLUtm                0.000000e+00   1666.666667
    10FTHFtl                   -1.666667e+03      0.000000
    10FTHFtm                   -1.000000e+04   3284.464975
    13DAMPPOX                   6.568857e-14   1197.358346
    ...                                  ...           ...
    NCAMUP                      0.000000e+00  10000.000000
    NCAMDe                      0.000000e+00  10000.000000
    sink_ncam_LSQBKT_c_RSQBKT_ -6.521726e-02      0.019581
    EX_pnto_R_LPAREN_e_RPAREN_ -3.171588e-01     -0.049167
    DM_pnto_R                   4.916697e-02      0.317159
    
    [2853 rows x 2 columns]
    FVA results for modified model:
                                     minimum       maximum
    10FTHF7GLUtl                0.000000e+00   1666.666667
    10FTHF7GLUtm                0.000000e+00   1666.666667
    10FTHFtl                   -1.666667e+03      0.000000
    10FTHFtm                   -1.000000e+04   3284.464975
    13DAMPPOX                   6.568857e-14   1197.358346
    ...                                  ...           ...
    NCAMUP                      0.000000e+00  10000.000000
    NCAMDe                      0.000000e+00  10000.000000
    sink_ncam_LSQBKT_c_RSQBKT_ -6.521726e-02      0.019581
    EX_pnto_R_LPAREN_e_RPAREN_ -3.171588e-01     -0.049167
    DM_pnto_R                   4.916697e-02      0.317159
    
    [2853 rows x 2 columns]


### Affected reactions based on change in Flux


```python
# Identify significantly affected reactions
affected_reactions = []

for rxn in integrated_model.reactions:
    if original_solution.fluxes[rxn.id] != modified_solution.fluxes[rxn.id]:
        affected_reactions.append(rxn.id)

print("Significantly affected reactions due to fasting:")
print(affected_reactions)
```

    Significantly affected reactions due to fasting:
    []

