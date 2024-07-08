# GSMM 
## An Open Source Tool for Building and Analysing Context Specific Metabolic Models and Generate visualisations.

## Overview

GSMM is a Python-based toolkit for Genome Scale Metabolic Modeling, designed for [currently] analyzing and visualizing metabolic models during my PHCCO Internship, specifically tailored for integrating transcriptomics based expression data and conducting nutrient depletion studies [ongoing]. It utilizes computational biology techniques to analyze flux distributions, reaction deletions, and sink fluxes across context-specific metabolic models.

## Features

- **Model Integration:** Create Context Specific Metabolic Models with easy to use functions
- **Flux Analysis:** Analyze and visualize flux distributions and correlations between metabolic models.
- **Reaction Deletion Analysis:** Conduct reaction deletion studies and visualize results using cluster maps.
- **Sink Flux Analysis:** Analyze and visualize sink fluxes for metabolites of interest.
- **Visualization:** Generate high-quality plots for presentation and analysis.

## Documentation

For detailed documentation of the source code, please visit [GSMM Documentation](https://gsmm.readthedocs.io/en/latest/).

## File Structure

```
gsmm/
├── __init__.py
│
└── csm/
    ├── build_csm.py                  # Building metabolic models
    ├── analyse_csm.py                # Analyzing metabolic models
    ├── visualisation.py              # Visualizing analysis results generated from CSMs
    └── config.py                     # Constants to store model paths, metabolites of interest, etc
```

## License

This project is licensed under the [MIT License](https://github.com/KarthikDani/PHCCOProject/blob/main/gsmm/LICENSE.txt).

---