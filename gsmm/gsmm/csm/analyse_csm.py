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
from .config import *

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to load SBML models
def load_models(model_paths: Dict[str, str]) -> Dict[str, Model]:
    """
    Load SBML models from specified paths into COBRApy Model objects.

    Args:
        model_paths (Dict[str, str]): Dictionary where keys are model names and values are paths to SBML files.

    Returns:
        Dict[str, Model]: Dictionary mapping model names to corresponding COBRApy Model objects.

    Notes:
        This function iterates through the provided dictionary of model paths, attempts to load each SBML model
        using COBRApy's read_sbml_model function, and stores the loaded models in a dictionary. If any path is
        invalid or model loading fails, an error message is printed and that model is skipped.
    """
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
    """
    Optimize the given COBRApy model and retrieve flux distributions.

    Args:
        model (Model): COBRApy Model object to be optimized.

    Returns:
        Dict[str, float]: Dictionary mapping reaction IDs to their corresponding flux values.

    Notes:
        This function attempts to optimize the input model using COBRApy's default solver.
        If successful, it returns a dictionary containing reaction IDs as keys and their respective flux values.
        If optimization fails due to an OptimizationError or Infeasible solution, an empty dictionary is returned.
    """
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
    """
    Extract fluxes for all reactions from a set of COBRApy models.

    Args:
        models (Dict[str, Model]): Dictionary of model names mapped to COBRApy Model objects.

    Returns:
        pd.DataFrame: DataFrame containing extracted flux data with columns 'Reaction', 'Flux', and 'Model'.

    Notes:
        This function iterates over each model in the input dictionary and optimizes it to obtain fluxes
        for all reactions in the model. The extracted flux data is structured into a DataFrame, where each row
        represents a reaction with its associated flux and model name.
    """
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
def filter_non_zero_fluxes(df: pd.DataFrame, threshold: float = 1) -> pd.DataFrame:
    """
    Filter non-zero fluxes from a DataFrame based on a specified threshold.

    Args:
        df (pd.DataFrame): Input DataFrame containing flux data.
        threshold (float, optional): Threshold value for flux. Default is 1.

    Returns:
        pd.DataFrame: Filtered DataFrame containing reactions with fluxes above or equal to the threshold.

    Raises:
        ValueError: If the 'Flux' column is not present in the input DataFrame.

    Notes:
        This function filters the input DataFrame to retain only rows where the 'Flux' column
        meets or exceeds the specified threshold. Missing flux values are filled with 0.
    """
    if 'Flux' not in df.columns:
        raise ValueError("DataFrame must contain 'Flux' column")
    
    print(f"Filtering fluxes with threshold: {threshold}...")
    df['Flux'].fillna(0, inplace=True)
    filtered_df = df[df['Flux'] >= threshold]
    print(f"Filtering complete. Number of reactions above threshold: {len(filtered_df)}")
    return filtered_df


# Function to extract sink reaction fluxes
def get_sink_fluxes(model: Model, model_name: str, metabolites: List[str]) -> List[Dict[str, Union[str, float]]]:
    """
    Extract sink reaction fluxes for specified metabolites from a given COBRApy model.

    Args:
        model (Model): COBRApy Model object from which sink fluxes will be extracted.
        model_name (str): Name or identifier of the model for logging purposes.
        metabolites (List[str]): List of metabolite IDs or names for which sink fluxes will be extracted.

    Returns:
        List[Dict[str, Union[str, float]]]: List of dictionaries, each containing information about a sink reaction flux.
            Each dictionary typically includes keys 'Metabolite', 'Flux', and 'Context_Model'.

    Notes:
        This function iterates over each metabolite in the provided list of metabolites.
        For each metabolite, it identifies sink reactions (reactions starting with 'SK')
        and calculates their flux under optimal conditions using COBRApy's optimization.
    """
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
    """
    Collect sink fluxes for all context-specific models.

    Args:
        models (Dict[str, Model]): Dictionary containing context-specific models.
            Keys are model names and values are COBRApy Model objects.
        metabolites (List[str]): List of metabolites of interest for sink flux collection.

    Returns:
        pd.DataFrame: DataFrame containing collected sink flux data.
            Columns typically include 'model', 'metabolite', and 'flux'.

    Notes:
        This function iterates over each model in the provided dictionary of models.
        For each model, it collects sink fluxes for specified metabolites using the
        get_sink_fluxes function. The collected sink flux data from all models is
        concatenated into a single DataFrame and returned.
    """
    all_sink_fluxes = []
    for model_name, model in models.items():
        print(f"Collecting sink fluxes for model: {model_name}...")
        sink_fluxes = get_sink_fluxes(model, model_name, metabolites)
        all_sink_fluxes.extend(sink_fluxes)
    print("Sink flux collection complete.")
    return pd.DataFrame(all_sink_fluxes)


# Function to save data for visualization
def save_data_for_visualization(df_fluxes: pd.DataFrame, df_sink_fluxes: pd.DataFrame, filename: str = 'flux_data.pkl'):
    """
    Save flux and sink flux dataframes to pickle files for visualization.

    Args:
        df_fluxes (pd.DataFrame): DataFrame containing flux data to be saved.
        df_sink_fluxes (pd.DataFrame): DataFrame containing sink flux data to be saved.
        filename (str, optional): Name of the pickle file to save. Default is 'flux_data.pkl'.

    Notes:
        This function saves two DataFrames, df_fluxes and df_sink_fluxes, to pickle files.
        df_fluxes is saved directly to the specified filename.
        df_sink_fluxes is saved to a file with 'sink_flux' substituted for 'flux' in the filename.
    """
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
    """
    Perform single reaction deletion analysis on a given COBRApy model.

    Args:
        model (Model): COBRApy Model object on which single reaction deletion will be performed.

    Returns:
        pd.DataFrame: DataFrame containing the results of single reaction deletion.
            Index represents reaction IDs, and columns typically include 'growth' or similar metrics.

    Notes:
        This function utilizes COBRApy's single_reaction_deletion method to analyze the impact
        of deleting each individual reaction in the model on growth or another specified metric.
        It returns a DataFrame where each row corresponds to a reaction and columns represent
        the results of the deletion analysis.
    """
    deletion_results = single_reaction_deletion(model)
    return deletion_results

# Function to save reaction deletion results for each model
def save_reaction_deletion_results(models_dict: Dict[str, Model], output_dir: str) -> Dict[str, pd.DataFrame]:
    """
    Perform single reaction deletion on each model in models_dict and save the deletion results.

    Args:
        models_dict (Dict[str, Model]): Dictionary containing models to analyze.
            Keys are model names and values are COBRApy Model objects.
        output_dir (str): Directory path where the deletion results will be saved.

    Returns:
        Dict[str, pd.DataFrame]: Dictionary containing deletion results for each model.
            Keys are model names and values are DataFrames with deletion results.

    Notes:
        This function performs single reaction deletion for each model provided in models_dict.
        It saves the deletion results as pickle files in the specified output_dir with filenames
        formatted as "{model_name}_reaction_deletion_results.pkl".
    """
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
    """
    Find common reactions across all models based on their deletion results.

    Args:
        reaction_deletion_results (Dict[str, pd.DataFrame]): Dictionary containing deletion results for each model.
            Keys are model names and values are DataFrames with deletion results.

    Returns:
        List[str]: List of reaction IDs that are common across all models.

    Notes:
        This function identifies common reactions present in the deletion results of all models provided.
        It computes the intersection of reaction sets from all model deletion results to find common reactions.
    """
    reaction_sets: List[Set[str]] = [set(results.index) for results in reaction_deletion_results.values()]
    common_reactions: Set[str] = set.intersection(*reaction_sets)
    return list(common_reactions)

# Function to filter deletion results to keep only common reactions
def filter_deletion_results(deletion_results: Dict[str, pd.DataFrame], common_reactions: List[str]) -> Dict[str, pd.DataFrame]:
    """
    Filter deletion results to keep only common reactions across models.

    Args:
        deletion_results (Dict[str, pd.DataFrame]): Dictionary containing deletion results for each model.
            Keys are model names and values are DataFrames with deletion results.
        common_reactions (List[str]): List of reaction IDs representing common reactions to retain.

    Returns:
        Dict[str, pd.DataFrame]: Filtered deletion results dictionary.
            Keys remain the same as deletion_results, and values are DataFrames filtered to include only common reactions.

    Notes:
        This function filters deletion results for each model provided in deletion_results based on the list of common reactions.
        Only reactions present in common_reactions are retained in each model's deletion results DataFrame.
    """
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

# def run_single_reaction_deletion():

#     # Load models
#     models = load_models(model_paths)

#     # Perform and save reaction deletion results
#     reaction_deletion_results = save_reaction_deletion_results(models, output_dir)

#     # Find common reactions
#     common_reactions = find_common_reactions(reaction_deletion_results)

#     # Filter deletion results to keep only common reactions
#     filtered_reaction_deletion_results = filter_deletion_results(reaction_deletion_results, common_reactions)

#     # Save filtered reaction deletion results
#     save_filtered_reaction_deletion_results(filtered_reaction_deletion_results, output_dir__reaction_deletion_results)

def analyse_and_save_fluxes(model_paths: Dict[str, str]) -> None:
    """
    Analyze fluxes from given model paths and save the results.

    Parameters:
    - model_paths (dict, optional): Dictionary containing model names as keys and file paths as values. Defaults to None.

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