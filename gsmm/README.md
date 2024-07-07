# GSSM 
## An Open Source Tool for Building and Analysing Context Specific Metabolic Models and Generate visualisations.

## Overview

GSMM is a Python-based toolkit designed for analyzing and visualizing metabolic models during my PHCCO Internship, specifically tailored for integrating transcriptomics based expression data and conducting nutrient depletion studies [ongoing]. It utilizes computational biology techniques to analyze flux distributions, reaction deletions, and sink fluxes across context-specific metabolic models.

## Features

- **Model Integration:** Create Context Specific Metabolic Models with easy to use functions
- **Flux Analysis:** Analyze and visualize flux distributions and correlations between metabolic models.
- **Reaction Deletion Analysis:** Conduct reaction deletion studies and visualize results using cluster maps.
- **Sink Flux Analysis:** Analyze and visualize sink fluxes for metabolites of interest.
- **Visualization:** Generate high-quality plots for presentation and analysis.

## File Structure

```
gsmm/
│
├── README.md             # Project overview and instructions
├── license.txt           # License information
├── setup.py              # Setup script for installation
│
└── gsmm/
    ├── __init__.py
    │
    └── csm/
        ├── build_csm.py                  # Scripts for building metabolic models
        ├── analyse_csm.py                # Scripts for analyzing metabolic models
        ├── visualisation.py              # Scripts for visualizing analysis results
        ├── flux_data.pkl                 # Example flux data (pickle format)
        ├── sink_flux_data.pkl            # Example sink flux data (pickle format)
        └── Model_paths_and_metabolites.py# File to store model paths and metabolites of interest
```

### Example Usage

```python
from gsmm.csm.visualisation import plot_fluxes, plot_reaction_deletion_clustermap

# Plot fluxes from default data file
plot_fluxes(flux_filepath='flux_data.pkl')

# Plot reaction deletion clustermap
plot_reaction_deletion_clustermap('ModelE_reaction_deletion_results.pkl', metabolites_of_interest, 'reaction_deletion_clustermap.png')
```

## License

This project is licensed under the [MIT License](LICENSE.txt).

---