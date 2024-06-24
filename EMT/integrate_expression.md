```python
import cobra
import pandas as pd
```

### Load EMT Model and the Expression data


```python
# Load the SBML model
model = cobra.io.read_sbml_model(
    "/Users/karthik/Desktop/PHCCO IISc Internship/EMT/MODEL1602080000_url.xml")

# Load the expression data
expression_data = pd.read_excel(
    "/Users/karthik/Desktop/PHCCO IISc Internship/EMT/example.xlsx", 
    sheet_name="Sheet1", header=None)
```

### Specify the column names in a sequence


```python
expression_data.columns = ['uniprot_id', 'expression_data', 'p_value']
```

### Convert to expression data to dict for easy handling


```python
# Convert the expression data to a dictionary
expression_dict = expression_data.set_index('uniprot_id')['expression_data'].to_dict()
expression_dict
```




    {'P60484': 4.6975102,
     'Q92934': 4.477713,
     'P28062': 4.448816,
     'P06241': 4.3806553,
     'O60674': 4.317239,
     'P27361': 4.2146053,
     'Q05655': 4.1526284,
     'Q6IAA8': 4.1215506,
     'P29353': -2.2021408,
     'P18433': -2.4953408,
     'Q5T1C6': -2.6014206,
     'Q8NDV7': -2.614702}



### Integrate expression
Returns all the affected up and down reactions, using a predefined threshold, as tuple of lists of the each.


```python
# Define function to integrate expression data
def integrate_expression(model, expression_dict, up_threshold=2, down_threshold=-2):
    up_genes = [gene for gene, expr in expression_dict.items() if expr >= up_threshold]
    down_genes = [gene for gene, expr in expression_dict.items() if expr <= down_threshold]

    affected_up_rxns = []
    affected_down_rxns = []

    for rxn in model.reactions:
        genes = [g.id for g in rxn.genes]
        if any(gene in up_genes for gene in genes):
            affected_up_rxns.append(rxn.id)
        if any(gene in down_genes for gene in genes):
            affected_down_rxns.append(rxn.id)

    return affected_up_rxns, affected_down_rxns

# Identify affected reactions
up_rxns, down_rxns = integrate_expression(model, expression_dict)
```

### Adjust bounds
Modified upper and lower bounds of reactions based on the up and down affected reactions list above. They are scaled by `factor`


```python
# Define function to adjust bounds
def adjust_bounds(model, up_rxns, down_rxns, factor=0.01):
    for rxn_id in up_rxns:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.upper_bound *= factor
        rxn.lower_bound *= factor

    for rxn_id in down_rxns:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.upper_bound *= factor
        rxn.lower_bound *= factor

# Adjust reaction bounds based on expression data
adjust_bounds(model, up_rxns, down_rxns)
```

### Save the integrated model


```python
# Save the modified model
cobra.io.write_sbml_model(model, 'integrated_MODEL1602080000.xml')
```
