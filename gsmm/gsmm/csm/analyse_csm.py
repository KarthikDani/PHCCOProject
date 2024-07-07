# flux_analysis.py

from typing import Dict, List, Union
from cobra.io import read_sbml_model
from cobra.core import Model
import pandas as pd
import numpy as np
from cobra.exceptions import OptimizationError, Infeasible
from cobra.flux_analysis import single_reaction_deletion
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to load SBML models
def load_models(model_paths: Dict[str, str]) -> Dict[str, Model]:
    models = {}
    for name, path in model_paths.items():
        print(f"Loading model: {name} from {path}...")
        if not os.path.exists(path):
            logging.error(f"File not found: {path}")
            print(f"Error: File not found: {path}")
            continue
        try:
            models[name] = read_sbml_model(path)
            logging.info(f"Loaded model: {name} from {path}")
            print(f"Successfully loaded model: {name}")
        except Exception as e:
            logging.error(f"Error loading model {name}: {e}")
            print(f"Error loading model {name}: {e}")
    return models


# Function to optimize the model and get flux distributions
def optimize_and_get_fluxes(model: Model) -> Dict[str, float]:
    try:
        print(f"Optimizing model: {model.id}...")
        solution = model.optimize()
        print(f"Optimization completed for model: {model.id}")
        return solution.fluxes.to_dict()
    except (OptimizationError, Infeasible) as e:
        logging.error(f"Optimization failed for model {model.id}: {e}")
        print(f"Optimization failed for model {model.id}: {e}")
        return {}


# Function to extract fluxes from models
def extract_fluxes(models: Dict[str, Model]) -> pd.DataFrame:
    fluxes_data = []
    for model_name, model in models.items():
        print(f"Extracting fluxes for model: {model_name}...")
        fluxes = optimize_and_get_fluxes(model)
        for reaction, flux in fluxes.items():
            fluxes_data.append({
                'Reaction': reaction,
                'Flux': flux,
                'Model': model_name
            })
    print("Flux extraction complete.")
    return pd.DataFrame(fluxes_data)


# Function to filter non-zero fluxes
def filter_non_zero_fluxes(df: pd.DataFrame, threshold: float = 50) -> pd.DataFrame:
    if 'Flux' not in df.columns:
        raise ValueError("DataFrame must contain 'Flux' column")
    print(f"Filtering fluxes with threshold: {threshold}...")
    df['Flux'].fillna(0, inplace=True)
    filtered_df = df[df['Flux'] >= threshold]
    print(f"Filtering complete. Number of reactions above threshold: {len(filtered_df)}")
    return filtered_df


# Function to extract sink reaction fluxes
def get_sink_fluxes(model: Model, model_name: str, metabolites: List[str]) -> List[Dict[str, Union[str, float]]]:
    sink_fluxes = []
    for met_id in metabolites:
        try:
            met = model.metabolites.get_by_id(met_id)
            for reaction in met.reactions:
                if reaction.id.startswith("SK"):
                    with model:
                        solution = model.optimize()
                        flux = solution.fluxes.get(reaction.id, 0)
                        sink_fluxes.append({
                            "Metabolite": met.name,
                            "Flux": flux,
                            "Context_Model": model_name,
                        })
        except KeyError:
            logging.warning(f"Metabolite {met_id} not found in model {model_name}")
            print(f"Warning: Metabolite {met_id} not found in model {model_name}")
        except (OptimizationError, Infeasible) as e:
            logging.error(f"Optimization failed for reaction {reaction.id} in model {model_name}: {e}")
            print(f"Optimization failed for reaction {reaction.id} in model {model_name}: {e}")
    print(f"Extracted sink fluxes for model: {model_name}")
    return sink_fluxes


# Function to collect sink fluxes for all context-specific models
def collect_sink_fluxes(models: Dict[str, Model], metabolites: List[str]) -> pd.DataFrame:
    all_sink_fluxes = []
    for model_name, model in models.items():
        print(f"Collecting sink fluxes for model: {model_name}...")
        sink_fluxes = get_sink_fluxes(model, model_name, metabolites)
        all_sink_fluxes.extend(sink_fluxes)
    print("Sink flux collection complete.")
    return pd.DataFrame(all_sink_fluxes)


# Function to save data for visualization
def save_data_for_visualization(df_fluxes: pd.DataFrame, df_sink_fluxes: pd.DataFrame, filename: str = 'flux_data.pkl'):
    try:
        print(f"Saving data to {filename}...")
        df_fluxes.to_pickle(filename)
        df_sink_fluxes.to_pickle(filename.replace('flux', 'sink_flux'))
        logging.info(f"Data saved to {filename} and {filename.replace('flux', 'sink_flux')}")
        print(f"Data saved successfully to {filename} and {filename.replace('flux', 'sink_flux')}")
    except Exception as e:
        logging.error(f"Failed to save data: {e}")
        print(f"Error: Failed to save data: {e}")


def main():
    # Define model paths
    model_paths = {
        "ModelE": "/Users/karthik/Desktop/PHCCO IISc Internship/Models/Epithelial_csm.xml",
        "ModelM": "/Users/karthik/Desktop/PHCCO IISc Internship/Models/Mesenchymal_csm.xml",
        "ModelMF": "/Users/karthik/Desktop/PHCCO IISc Internship/Models/mesenchymal_fasting_integrated_csm.xml",
    }

    # List of metabolites of interest
    metabolites_of_interest = [
        "ala__L_c", "arg__L_c", "asn__L_c", "asp__L_c", "cys__L_c", "gln__L_c",
        "glu__L_c", "gly_c", "his__L_c", "ile__L_c", "leu__L_c", "lys__L_c",
        "met__L_c", "phe__L_c", "pro__L_c", "ser__L_c", "thr__L_c", "trp__L_c",
        "tyr__L_c", "val__L_c", "damp_c", "dcmp_c", "dgmp_c", "dtmp_c",
        "cmp_c", "gmp_c", "ump_c", "amp_c", "glygn2_c", "sphmyln_hs_c",
        "chsterol_c", "xolest_hs_c", "mag__hs_c", "dag_hs_c", "pail_hs_c",
        "pe_hs_c", "ps_hs_c", "pchol_hs_c", "lpchol_hs_c", "clpn_hs_c",
        "pa_hs_c", "hdcea_c", "hdca_c", "ocdcea_c", "ocdca_c", "ptrc_c",
        "spmd_c", "sprm_c", "gthrd_c", "nad_c", "nadp_c", "Q10_c",
        "paps_c", "thbpt_c", "crn_c", "atp_c", "adp_c", "pi_c",
        "h2o_c", "h_c",
    ]

    # Load models
    models = load_models(model_paths)

    # Extract and filter fluxes
    df_fluxes = extract_fluxes(models)
    df_fluxes_filtered = filter_non_zero_fluxes(df_fluxes)

    # Collect sink fluxes
    df_sink_fluxes = collect_sink_fluxes(models, metabolites_of_interest)

    # Save data for visualization
    save_data_for_visualization(df_fluxes_filtered, df_sink_fluxes)


if __name__ == "__main__":
    main()
