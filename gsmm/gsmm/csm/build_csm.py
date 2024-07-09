import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import remove_genes, prune_unused_metabolites, prune_unused_reactions
from corda import CORDA
from typing import Union, Optional
import cobra
from .config import *

def set_glpk_solver() -> None:
    """
    Sets the default solver in COBRApy to GLPK (GNU Linear Programming Kit).

    Raises:
        RuntimeError: If setting the GLPK solver fails for any reason.

    Notes:
        This function updates the global solver configuration for COBRApy.
        GLPK is chosen as it is an open-source linear programming solver suitable
        for solving many types of optimization problems encountered in constraint-based
        modeling of metabolic networks.
    """
    try:
        # Get the default solver configuration from COBRApy
        config = cobra.Configuration()

        # Set the solver to GLPK
        config.solver = "glpk"

        # Update the solver configuration in COBRApy
        cobra.Configuration._config = config
        print("Cobra Configuration set to GLPK solver!")

    except Exception as e:
        raise RuntimeError(f"Failed to set GLPK solver: {str(e)}")

def read_parent_model(model_path: str) -> cobra.Model:
    """
    Reads an SBML model from the specified path using COBRApy.

    Args:
        model_path (str): Path to the SBML model file path to reconstruct from. This is usully base model.

    Returns:
        cobra.Model: The COBRApy model object loaded from the SBML file.

    Raises:
        FileNotFoundError: If the specified model file path does not exist.
        IOError: If there is an error reading the SBML model file.

    Notes:
        This function uses COBRApy's `read_sbml_model` function to read an SBML model
        file and return the corresponding COBRApy model object.
    """
    print(f"Loading SBML model from {model_path}...")
    model = read_sbml_model(model_path)
    print("SBML model loaded.")
    return model

def load_expression_data(data_path: str) -> pd.DataFrame:
    """
    Loads expression data from a CSV file into a pandas DataFrame.

    Args:
        data_path (str): Path to the CSV file containing expression data.

    Returns:
        pd.DataFrame: DataFrame containing the loaded expression data.

    Raises:
        FileNotFoundError: If the specified data file path does not exist.
        IOError: If there is an error reading the CSV file.

    Notes:
        This function uses pandas' `read_csv` function to read expression data from
        a CSV file and returns it as a DataFrame. The CSV file should contain columns
        representing Gene IDs and corresponding expression values.
    """
    print(f"Loading expression data from {data_path}...")
    expression_data = pd.read_csv(data_path)
    print("Expression data loaded.")
    return expression_data

def extract_genes(data: pd.DataFrame, id_column: Union[str, None]) -> list:
    """
    Extracts gene IDs or names from a pandas DataFrame column.

    Args:
        data (pd.DataFrame): DataFrame containing gene data.
        id_column (Union[str, None]): Name of the column containing gene IDs or names.

    Returns:
        list: List of gene IDs or names extracted from the specified column.

    Raises:
        ValueError: If `id_column` is not provided or does not exist in the DataFrame.

    Notes:
        This function extracts gene IDs or names from a specified column in the provided
        DataFrame. It drops missing values (NaNs) and returns a list of extracted gene IDs
        or names.
    """
    # if name_column:
    #     print(f"Extracting genes from {name_column} column...")
    #     genes = data[name_column].dropna().tolist()
    if id_column:
        print(f"Extracting genes from {id_column} column...")
        genes = data[id_column].dropna().tolist()
    else:
        raise ValueError("gene_id_column from the expression data must be provided.")
    
    print(f"Number of Genes extracted: {len(genes)}")
    return genes

def filter_model_by_genes(model: cobra.Model, genes: list) -> cobra.Model:
    """
    Filters a COBRApy model by removing genes not present in a specified list.

    Args:
        model (cobra.Model): The COBRApy model to be filtered.
        genes (list): List of gene IDs or names to retain in the model.

    Returns:
        cobra.Model: Filtered COBRApy model containing only the specified genes.

    Notes:
        This function creates a copy of the input model and removes genes that are not
        present in the provided list from the copied model. It also prunes unused reactions
        and metabolites to ensure consistency after gene removal.
    """
    print("Filtering model by genes...")
    new_model = model.copy()
    genes_in_model = {gene.id for gene in new_model.genes} 
    genes_to_remove = list(genes_in_model - set(genes))
    
    genes_to_remove_by_id = {gene.id for gene in new_model.genes if gene.id in genes_to_remove}
    
    remove_genes(new_model, genes_to_remove_by_id, remove_reactions=True)
    new_model, _ = prune_unused_reactions(new_model)
    new_model, _ = prune_unused_metabolites(new_model)
    
    print("Model filtered.")
    return new_model

def normalize_expression_data(expression_data: pd.DataFrame, scores_column: str) -> pd.DataFrame:
    """
    Normalizes expression data and assigns confidence levels based on quantiles.

    Args:
        expression_data (pd.DataFrame): DataFrame containing expression data.
        scores_column (str): Name of the column to normalize and use for confidence levels.

    Returns:
        pd.DataFrame: DataFrame with normalized scores and assigned confidence levels.

    Notes:
        This function uses MinMaxScaler from scikit-learn to normalize values in the specified
        column of the input DataFrame. It calculates confidence levels based on quantiles of
        the normalized scores.
    """
    print("Normalizing expression data...")
    scaler = MinMaxScaler()
    expression_data[f'Normalized_{scores_column}'] = scaler.fit_transform(expression_data[[scores_column]])
    
    print("Expression data has been normalized.")
    return expression_data

def assign_reaction_confidences(model: cobra.Model, 
                                expression_data: pd.DataFrame, 
                                gene_id_column: str, 
                                scores_column: str) -> pd.DataFrame:
    """
    Assigns confidence levels to reactions based on associated genes' expression data.

    Args:
        model (cobra.Model): The metabolic model to assign reaction confidence levels to.
        expression_data (pd.DataFrame): DataFrame containing expression data and gene confidence levels.
        gene_id_column (str): Column name in `expression_data` containing gene IDs or names.
        scores_column (str): Name of the column in `expression_data` used for gene expression values.

    Returns:
        pd.DataFrame: DataFrame with assigned confidence levels for reactions.

    Notes:
        This function assigns confidence levels (ranging from 0 to 3) to reactions in the provided metabolic model 
        based on the associated genes' expression levels in `expression_data`. Specific biomass reactions are assigned
        the highest confidence level (3) by default.

        The number of different confidence levels assigned depends on the quartiles of `scores_column` in 
        `expression_data`, distributed as follows:
        - (-1): Less than 25th percentile
        - 0: 25th to 49.99th percentile
        - 1: 50th to 74.99th percentile
        - 2: 75th to 89.99th percentile
        - 3: >= 90th percentile

    Example:
        If `scores_column` represents gene expression levels, the function assigns the confidence levels based on
        quartiles of these expression values.
    """
    print("Assigning confidence levels to reactions giving highest confidence levels to Biomass reactions...")
    
    expression_data['Gene_Confidence_Level'] = expression_data[f'Normalized_{scores_column}'].apply(lambda value: (
    3 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.90) else
    2 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.75) else
    1 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.50) else
    0 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.25) else
    -1
    ))
    
    reaction_confidence = {}
    
    for reaction in model.reactions:
        associated_genes = {gene.id for gene in reaction.genes}
        gene_confidences = expression_data.loc[expression_data[gene_id_column].isin(associated_genes), 'Gene_Confidence_Level']
        reaction_confidence[reaction.id] = gene_confidences.min() if not gene_confidences.empty else -1
    
    # Assign highest confidence level to specific biomass reactions
    for biomass_reaction_id in ['BIOMASS_maintenance_noTrTr', 'BIOMASS_maintenance', 'BIOMASS_reaction']:
        if biomass_reaction_id in reaction_confidence:
            reaction_confidence[biomass_reaction_id] = 3
    
    reaction_confidence_df = pd.DataFrame(reaction_confidence.items(), columns=['Reaction_ID', 'Confidence_Level'])
    
    # Get count of each confidence level after assignment
    counts = reaction_confidence_df['Confidence_Level'].value_counts().to_dict()
    for level in range(-1, 4):
        print(f"Number of reactions with confidence level {level}: {counts.get(level, 0)}")
    
    reaction_confidence_df.to_csv(f"reaction_{scores_column}_confidence_levels.csv", index=False)
    print(f"Reaction confidence levels saved as reaction_{scores_column}_confidence_levels.csv.")
    return reaction_confidence_df

def optimize_model(model: cobra.Model, reaction_confidence_dict: dict) -> CORDA:
    """
    Initializes and optimizes a COBRApy model using CORDA with reaction confidence levels.

    Args:
        model (cobra.Model): The COBRApy model to be optimized.
        reaction_confidence_dict (dict): Dictionary mapping reaction IDs to confidence levels.

    Returns:
        CORDA: Optimized CORDA object containing the optimized model and solution.

    Notes:
        This function initializes a CORDA optimization object using the provided COBRApy model
        and reaction confidence levels. CORDA (Cost Optimisation Reaction Dependency Assessment)
        optimizes metabolic models based on context-specific data and reactions' confidence levels.
    """
    print("Building CORDA model, Please be patient as this will take some time...")
    opt = CORDA(model, reaction_confidence_dict)
    opt.build()
    print("CORDA model optimization completed.")
    return opt

def main(model_path: str, data_path: str, gene_id_column: Union[str, None], scores_column: str, base_model_path: str):
    """
    Executes a pipeline to load, filter, normalize data, assign reaction confidences, optimize a metabolic model using CORDA.

    Args:
        model_path (str): Path to the SBML model file (Using big model like Recon3D).
        data_path (str): Path to the CSV file containing expression data.
        gene_id_column (Union[str, None]): Column name in `expression_data` containing gene IDs or names.
        scores_column (str): Name of the column in `expression_data` used for gene confidence levels.
        base_model_path (str): Path to save the filtered SBML model.

    Returns:
        cobra.Model: Optimized COBRApy model (Context Specific Model) object after CORDA optimization.

    Notes:
        This function orchestrates the entire process of setting up the solver, loading the SBML model,
        loading expression data, extracting genes, filtering the model based on genes, normalizing expression
        data, assigning reaction confidence levels based on gene expression, optimizing the model using CORDA,
        and saving the optimized model to a new SBML file.
    """
    set_glpk_solver()
    
    # Load SBML model and expression data
    model = read_parent_model(model_path)
    expression_data = load_expression_data(data_path)
    
    # Extract genes from expression data
    genes = extract_genes(expression_data, gene_id_column)
    
    # Filter model by genes
    filtered_model = filter_model_by_genes(model, genes)
    write_sbml_model(filtered_model, base_model_path)
    
    # Normalize expression data and assign confidence levels
    normalized_expression = normalize_expression_data(expression_data, scores_column)
    
    # Assign reaction confidence levels
    reaction_confidence_df = assign_reaction_confidences(filtered_model, normalized_expression, gene_id_column, 
                                                         scores_column)
    
    # Convert reaction confidence DataFrame to dictionary
    reaction_confidence_dict = reaction_confidence_df.set_index('Reaction_ID')['Confidence_Level'].to_dict()
    
    # Optimize model using CORDA
    optimized_model = optimize_model(filtered_model, reaction_confidence_dict)
    
    # Print optimized model
    print("Optimized CORDA model:")
    print(optimized_model)
    
    return optimized_model.cobra_model(scores_column)

def run_model_reconstruction(model_path: str,
                             base_model_path: str,
                             data_path: str,
                             gene_id_column: str,
                             scores_column: str,
                             ) -> Optional[cobra.Model]:
    """
    Runs the pipeline for reconstructing an optimized metabolic model from expression data.

    Args:
        model_path (str): Path to the Parent SBML model file from which Base model will be derived.
        base_model_path (str): Path to save the filtered SBML Base model which will be used for Reconstruction of Context Specific Models.
        data_path (str): Path to the CSV file containing expression data.
        gene_id_column (Optional[str], optional): Column name in the expression data containing gene IDs or names.
        scores_column (str): Name of the column in the expression data used for gene confidence levels.

    Returns:
        Optional[cobra.Model]: Optimized COBRApy model object after reconstruction, or None if an error occurs.

    Notes:
        This function serves as an interface to run the main pipeline (`main` function) for reconstructing
        an optimized metabolic model from the provided expression data. It handles exceptions and prints
        error messages if reconstruction fails, returning None in case of errors.
    """
    try:
        optimized_model = main(model_path, data_path, gene_id_column, scores_column, base_model_path)
        print(f"Reconstruction completed from the parent model {optimized_model}")
        return optimized_model
    except Exception as e:
        print(f"Error during model reconstruction: {e}")
        return None