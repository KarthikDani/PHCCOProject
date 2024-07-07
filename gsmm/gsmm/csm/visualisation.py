# visualisation.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from typing import Optional

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def load_data(filepath: str) -> Optional[pd.DataFrame]:
    """Load data from a pickle file."""
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

def plot_flux_distribution(df_fluxes, save_path):
    try:
        logging.info(f"Generating flux distribution plot...")
        
        # Check if the required columns are present
        if not all(col in df_fluxes.columns for col in ['Model', 'Reaction', 'Flux']):
            raise ValueError("Required columns ('Model', 'Reaction', 'Flux') not found in the dataframe.")
        
        # Plotting flux distributions across different models and reactions
        plt.figure(figsize=(12, 8))
        sns.violinplot(data=df_fluxes, x='Model', y='Flux', hue='Reaction', split=True, inner='quartile')
        plt.title('Flux Distribution across Models for Different Reactions')
        plt.xlabel('Model')
        plt.ylabel('Flux')
        plt.xticks(rotation=45)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(save_path)
        logging.info(f"Flux distribution plot saved as {save_path}")
        plt.show()
    
    except Exception as e:
        logging.error(f"Error generating flux distribution plot: {str(e)}")

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
        plt.title("Sink Fluxes Heatmap")
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
        logging.info(f"Sink fluxes heatmap saved as {output_file}")
    except Exception as e:
        logging.error(f"Error generating sink fluxes heatmap: {e}")

def main():
    flux_filepath = 'flux_data.pkl'
    sink_flux_filepath = 'sink_flux_data.pkl'

    df_flux = load_data(flux_filepath)
    df_sink_fluxes = load_data(sink_flux_filepath)

    if df_flux is not None:
        plot_flux_distribution(df_flux, 'flux_distribution.png')

    if df_sink_fluxes is not None:
        plot_sink_fluxes_heatmap(df_sink_fluxes, 'sink_fluxes_heatmap.png')

if __name__ == "__main__":
    main()
