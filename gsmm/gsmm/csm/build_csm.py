import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import remove_genes, prune_unused_metabolites, prune_unused_reactions
from corda import CORDA
from typing import Union, Optional
import cobra

def read_parent_model(model_path: str) -> cobra.Model:
    """Load a COBRA model from SBML file."""
    print(f"Loading SBML model from {model_path}...")
    model = read_sbml_model(model_path)
    print("SBML model loaded.")
    return model

def load_expression_data(data_path: str) -> pd.DataFrame:
    """Load expression data from CSV."""
    print(f"Loading expression data from {data_path}...")
    expression_data = pd.read_csv(data_path)
    print("Expression data loaded.")
    return expression_data

def extract_genes(data: pd.DataFrame, name_column: Union[str, None], id_column: Union[str, None]) -> list:
    """Extract genes from expression data based on specified columns."""
    if name_column:
        print(f"Extracting genes from {name_column} column...")
        genes = data[name_column].dropna().tolist()
    elif id_column:
        print(f"Extracting genes from {id_column} column...")
        genes = data[id_column].dropna().tolist()
    else:
        raise ValueError("Either name_column or id_column must be provided.")
    
    print(f"Number of Genes extracted: {len(genes)}")
    return genes

def filter_model_by_genes(model: cobra.Model, genes: list) -> cobra.Model:
    """
    Filter a COBRA model to include only reactions associated with specified genes.
    """
    print("Filtering model by genes...")
    new_model = model.copy()
    genes_in_model = {gene.id if gene.id else gene.name for gene in new_model.genes}  # Handles both ID and name
    genes_to_remove = list(genes_in_model - set(genes))
    
    genes_to_remove_by_id = {gene.id for gene in new_model.genes if gene.id in genes_to_remove or gene.name in genes_to_remove}
    
    remove_genes(new_model, genes_to_remove_by_id, remove_reactions=True)
    new_model, _ = prune_unused_reactions(new_model)
    new_model, _ = prune_unused_metabolites(new_model)
    
    print("Model filtered.")
    return new_model

def normalize_expression_data(expression_data: pd.DataFrame, scores_column: str) -> pd.DataFrame:
    """
    Normalize expression data and calculate confidence levels.
    """
    print("Normalizing expression data...")
    scaler = MinMaxScaler()
    expression_data[f'Normalized_{scores_column}'] = scaler.fit_transform(expression_data[[scores_column]])
    
    expression_data['Gene_Confidence_Level'] = expression_data[f'Normalized_{scores_column}'].apply(lambda value: (
        3 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.90) else
        2 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.75) else
        1 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.50) else
        0 if value >= expression_data[f'Normalized_{scores_column}'].quantile(0.25) else
        -1
    ))
    
    print("Expression data normalized and confidence levels assigned.")
    return expression_data

def assign_reaction_confidences(model: cobra.Model, expression_data: pd.DataFrame, gene_id_column: str, scores_column: str) -> pd.DataFrame:
    """
    Assign confidence levels to reactions based on associated genes.
    """
    print("Assigning confidence levels to reactions...")
    reaction_confidence = {}
    
    for reaction in model.reactions:
        reaction_id = reaction.id
        associated_genes = {gene.name if gene.name else gene.id for gene in reaction.genes}
        gene_confidences = expression_data.loc[expression_data[gene_id_column].isin(associated_genes), 'Gene_Confidence_Level']
        
        reaction_confidence[reaction_id] = gene_confidences.min() if not gene_confidences.empty else -1
    
    # Assign highest confidence level to specific biomass reactions
    for biomass_reaction_id in ['BIOMASS_maintenance_noTrTr', 'BIOMASS_maintenance', 'BIOMASS_reaction']:
        if biomass_reaction_id in reaction_confidence:
            reaction_confidence[biomass_reaction_id] = 3
    
    reaction_confidence_df = pd.DataFrame(reaction_confidence.items(), columns=['Reaction_ID', 'Confidence_Level'])
    print("Reaction Confidences: ", reaction_confidence_df)
    reaction_confidence_df.to_csv(f"reaction_{scores_column}_confidence_levels.csv", index=False)
    print(f"Reaction confidence levels saved as reaction_{scores_column}_confidence_levels.csv.")
    return reaction_confidence_df

def optimize_model(model: cobra.Model, reaction_confidence_dict: dict) -> CORDA:
    """
    Initialize and build CORDA model optimization.
    """
    print("Initializing and building CORDA model...")
    opt = CORDA(model, reaction_confidence_dict)
    opt.build()
    print("CORDA model optimization completed.")
    return opt

def main(model_path: str, data_path: str, gene_name_column: Union[str, None], gene_id_column: Union[str, None], scores_column: str):
    """
    Main function to execute the workflow.
    """
    # Load SBML model and expression data
    model = read_parent_model(model_path)
    expression_data = load_expression_data(data_path)
    
    # Extract genes from expression data
    genes = extract_genes(expression_data, gene_name_column, gene_id_column)
    
    # Filter model by genes
    filtered_model = filter_model_by_genes(model, genes)
    write_sbml_model(filtered_model, "emt_base.xml")
    
    # Normalize expression data and assign confidence levels
    normalized_expression = normalize_expression_data(expression_data, scores_column)
    
    # Assign reaction confidence levels
    reaction_confidence_df = assign_reaction_confidences(filtered_model, normalized_expression, gene_id_column, scores_column)
    
    # Convert reaction confidence DataFrame to dictionary
    reaction_confidence_dict = reaction_confidence_df.set_index('Reaction_ID')['Confidence_Level'].to_dict()
    
    # Optimize model using CORDA
    optimized_model = optimize_model(filtered_model, reaction_confidence_dict)
    
    # Print optimized model
    print("Optimized CORDA model:")
    print(optimized_model)
    
    return optimized_model.cobra_model(scores_column)

def run_model_reconstruction(model_path: str, data_path: str,
                             gene_name_column: Optional[str],
                             gene_id_column: str,
                             scores_column: str) -> Optional[str]:
    """
    Run model reconstruction using the specified paths and columns.

    Parameters:
    - model_path (str): Path to the XML model file.
    - data_path (str): Path to the input data CSV file.
    - gene_name_column (str or None): Column name for gene names (optional).
    - gene_id_column (str): Column name for gene IDs.
    - scores_column (str): Column name for scores.

    Returns:
    - str or None: Optimized model name if reconstruction is successful, else None.
    """
    try:
        optimized_model = main(model_path, data_path, gene_name_column, gene_id_column, scores_column)
        print(f"Reconstruction completed from the parent model {optimized_model}")
        return optimized_model
    except Exception as e:
        print(f"Error during model reconstruction: {e}")
        return None