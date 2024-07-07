# flux_analysis.py

from typing import Dict, List, Union, Set, Optional
from cobra.io import read_sbml_model
from cobra.core import Model
import pandas as pd
import numpy as np
from cobra.exceptions import OptimizationError, Infeasible
from cobra.flux_analysis import single_reaction_deletion
import logging
import os
from config import *

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
        
        
# Single reaction deletions
# Function to perform single reaction deletion on a model
def perform_single_reaction_deletion(model: Model) -> pd.DataFrame:
    """Perform single reaction deletion on the given model."""
    deletion_results = single_reaction_deletion(model)
    return deletion_results

# Function to save reaction deletion results for each model
def save_reaction_deletion_results(models_dict: Dict[str, Model], output_dir: str) -> Dict[str, pd.DataFrame]:
    """Perform single reaction deletion for each model in models_dict and save the results."""
    results: Dict[str, pd.DataFrame] = {}
    for model_name, model in models_dict.items():
        logging.info(f"Performing single reaction deletion on {model_name}...")
        deletion_results = perform_single_reaction_deletion(model)
        results[model_name] = deletion_results
        output_path = f"{output_dir}/{model_name}_reaction_deletion_results.pkl"
        deletion_results.to_pickle(output_path)
        logging.info(f"Saved single reaction deletion results for {model_name} to {output_path}")
    return results

# Function to find common reactions across all models
def find_common_reactions(reaction_deletion_results: Dict[str, pd.DataFrame]) -> List[str]:
    """Find common reactions across all models."""
    reaction_sets: List[Set[str]] = [set(results.index) for results in reaction_deletion_results.values()]
    common_reactions: Set[str] = set.intersection(*reaction_sets)
    return list(common_reactions)

# Function to filter deletion results to keep only common reactions
def filter_deletion_results(deletion_results: Dict[str, pd.DataFrame], common_reactions: List[str]) -> Dict[str, pd.DataFrame]:
    """Filter reaction deletion results to keep only common reactions."""
    filtered_results: Dict[str, pd.DataFrame] = {}
    for model_name, results in deletion_results.items():
        filtered_results[model_name] = results.loc[results.index.isin(common_reactions)]
    return filtered_results

# Function to save filtered reaction deletion results
def save_filtered_reaction_deletion_results(filtered_results: Dict[str, pd.DataFrame], output_dir: str) -> None:
    """
    Save filtered reaction deletion results to pickle files.

    Parameters:
    - filtered_results (Dict[str, pd.DataFrame]): Dictionary mapping model names to filtered deletion results DataFrames.
    - output_dir (str): Directory path where to save the results.
    """
    try:
        for model_name, df in filtered_results.items():
            output_path = f"{output_dir}/{model_name}_filtered_reaction_deletion_results.pkl"
            df.to_pickle(output_path)
            print(f"Saved filtered reaction deletion results for {model_name} to {output_path}")
    except Exception as e:
        print(f"Error saving filtered reaction deletion results: {str(e)}")

def run_single_reaction_deletion():

    # Load models
    models = load_models(model_paths)

    # Perform and save reaction deletion results
    reaction_deletion_results = save_reaction_deletion_results(models, output_dir)

    # Find common reactions
    common_reactions = find_common_reactions(reaction_deletion_results)

    # Filter deletion results to keep only common reactions
    filtered_reaction_deletion_results = filter_deletion_results(reaction_deletion_results, common_reactions)

    # Save filtered reaction deletion results
    save_filtered_reaction_deletion_results(filtered_reaction_deletion_results, output_dir__reaction_deletion_results)

def analyse_and_save_fluxes(model_paths: Optional[Dict[str, str]] = None) -> None:
    """
    Analyze fluxes from given model paths and save the results.

    Parameters:
    - model_paths (dict, optional): Dictionary containing model names as keys and file paths as values.
                                    Defaults to None.

    Returns:
    - None
    """
    if not isinstance(model_paths, dict):
        raise TypeError("Expected model_paths to be a dictionary.")
    
    # Load models
    models = load_models(model_paths)

    # Extract and filter fluxes
    df_fluxes = extract_fluxes(models)
    df_fluxes_filtered = filter_non_zero_fluxes(df_fluxes)

    # Collect sink fluxes
    df_sink_fluxes = collect_sink_fluxes(models, metabolites_of_interest)

    # Save data for visualization
    save_data_for_visualization(df_fluxes_filtered, df_sink_fluxes)