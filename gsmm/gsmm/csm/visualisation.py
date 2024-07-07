import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from typing import Optional
from config import *
import pickle

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def load_data(filepath: str) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_pickle(filepath)
        logging.info(f"Loaded data from {filepath}")
        return df
    except FileNotFoundError:
        logging.error(f"File not found: {filepath}")
        return None
    except Exception as e:
        logging.error(f"Error loading data from {filepath}: {e}")
        return None

def plot_flux_distribution_clustermap(df_fluxes, save_path):
    """
    Generate a clustermap for flux distribution across different models and reactions.

    Parameters:
    - df_fluxes (pd.DataFrame): DataFrame with 'Model', 'Reaction', and 'Flux' columns.
    - save_path (str): Path to save the clustermap.
    """
    try:
        logging.info("Generating flux distribution clustermap...")

        # Check if the required columns are present
        if not all(col in df_fluxes.columns for col in ['Model', 'Reaction', 'Flux']):
            raise ValueError("Required columns ('Model', 'Reaction', 'Flux') not found in the dataframe.")

        # Pivot the dataframe for clustermap
        pivot_table = df_fluxes.pivot_table(index='Reaction', columns='Model', values='Flux')

        # Drop rows/columns with all NaN values which might cause issues in clustermap
        pivot_table.dropna(axis=0, how='all', inplace=True)
        pivot_table.dropna(axis=1, how='all', inplace=True)

        # Fill NaNs with zeros or appropriate value if required
        pivot_table.fillna(0, inplace=True)

        # Plotting clustermap
        sns.clustermap(pivot_table, cmap="coolwarm", center=0, figsize=(14, 10),
                       method='average', metric='euclidean', standard_scale=1)
        plt.title('Clustermap of Flux Distribution across Models and All Reactions')
        plt.tight_layout()

        # Save the plot
        plt.savefig(save_path)
        logging.info(f"Flux distribution clustermap saved as {save_path}")
        # plt.show()

    except Exception as e:
        logging.error(f"Error generating flux distribution clustermap: {str(e)}")


def plot_sink_fluxes_heatmap(df: pd.DataFrame, output_file: str) -> None:
    """
    Plot a heatmap for sink fluxes across different context models.

    Parameters:
    - df (pd.DataFrame): DataFrame containing sink flux data with columns 'Metabolite', 'Flux', 'Context_Model'.
    - output_file (str): Path to save the sink fluxes heatmap.
    """
    if df is None or df.empty:
        logging.warning("Data frame is empty or None. Skipping heatmap generation.")
        return

    logging.info("Generating sink fluxes heatmap...")
    try:
        pivot_table = df.pivot(index="Metabolite", columns="Context_Model", values="Flux")
        plt.figure(figsize=(12, 10))
        sns.heatmap(pivot_table, annot=True, cmap="coolwarm", center=0)
        plt.title("Sink Fluxes of Reactions associated with Metabolites of interest across all Models")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        logging.info(f"Sink fluxes heatmap saved as {output_file}")
    except Exception as e:
        logging.error(f"Error generating sink fluxes heatmap: {e}")

def plot_flux_correlation_heatmap(df_fluxes, save_path):
    """
    Generate a heatmap of correlation coefficients for fluxes between different models for all reactions.

    Parameters:
    - df_fluxes (pd.DataFrame): DataFrame with 'Model', 'Reaction', and 'Flux' columns.
    - save_path (str): Path to save the heatmap.
    """
    try:
        logging.info("Generating flux correlation heatmap...")

        # Pivot the dataframe for correlation analysis
        pivot_df = df_fluxes.pivot_table(index='Reaction', columns='Model', values='Flux')

        # Compute correlation matrix, but only for different models
        models = pivot_df.columns
        corr_matrix = pd.DataFrame(index=models, columns=models)

        for i in models:
            for j in models:
                if i != j:
                    corr_matrix.loc[i, j] = pivot_df[i].corr(pivot_df[j])

        # Plotting heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix.astype(float), annot=True, cmap="coolwarm", vmin=-1, vmax=1, center=0, square=True)
        plt.title('Correlation Heatmap of All Reaction Fluxes across Models')
        plt.tight_layout()

        # Save the plot
        plt.savefig(save_path)
        logging.info(f"Flux correlation heatmap saved as {save_path}")
        # plt.show()

    except Exception as e:
        logging.error(f"Error generating flux correlation heatmap: {str(e)}")

def plot_sink_flux_correlation_heatmap(df_sink_fluxes, save_path):
    """
    Generate a heatmap of correlation coefficients for sink fluxes between different models.

    Parameters:
    - df_sink_fluxes (pd.DataFrame): DataFrame with 'Metabolite', 'Flux', 'Context_Model' columns.
    - save_path (str): Path to save the heatmap.
    """
    try:
        logging.info("Generating sink flux correlation heatmap...")

        # Pivot the dataframe for correlation analysis
        pivot_df = df_sink_fluxes.pivot_table(index='Metabolite', columns='Context_Model', values='Flux')

        # Compute correlation matrix, but only for different context models
        models = pivot_df.columns
        corr_matrix = pd.DataFrame(index=models, columns=models)

        for i in models:
            for j in models:
                if i != j:
                    corr_matrix.loc[i, j] = pivot_df[i].corr(pivot_df[j])

        # Plotting heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix.astype(float), annot=True, cmap="coolwarm", vmin=-1, vmax=1, center=0, square=True)
        plt.title('Correlation Heatmap of Sink Fluxes across all Models')
        plt.tight_layout()

        # Save the plot
        plt.savefig(save_path)
        logging.info(f"Sink flux correlation heatmap saved as {save_path}")
        # plt.show()

    except Exception as e:
        logging.error(f"Error generating sink flux correlation heatmap: {str(e)}")
        
def filter_and_save_results(filtered_results, file_path):
    """Filter and save valid DataFrames to a pickle file."""

    filtered_results_clean = {}
    for key, df in filtered_results.items():
        if isinstance(df, pd.DataFrame) and set(['ids', 'growth', 'status']).issubset(df.columns):
            filtered_results_clean[key] = df
        else:
            logging.warning(f"Skipping '{key}' as it is not a valid DataFrame for saving.")

    with open(file_path, 'wb') as f:
        pickle.dump(filtered_results_clean, f)


def plot_fluxes(flux_filepath: Optional[str] = flux_filepath,
                sink_flux_filepath: Optional[str] = sink_flux_filepath):
    """
    Plot flux distributions and sink fluxes using default or provided file paths.

    Parameters:
    - flux_filepath (str, optional): File path to flux data (default: 'flux_data.pkl').
    - sink_flux_filepath (str, optional): File path to sink flux data (default: 'sink_flux_data.pkl').

    Returns:
    - None
    """
    df_flux = load_data(flux_filepath)
    df_sink_fluxes = load_data(sink_flux_filepath)

    if df_flux is not None:
        plot_flux_distribution_clustermap(df_flux, 'flux_distribution_clustermap.png')
        plot_flux_correlation_heatmap(df_flux, 'flux_correlation_heatmap.png')

    if df_sink_fluxes is not None:
        plot_sink_fluxes_heatmap(df_sink_fluxes, 'sink_fluxes_heatmap.png')
        plot_sink_flux_correlation_heatmap(df_sink_fluxes, 'sink_flux_correlation_heatmap.png')