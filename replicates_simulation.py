import pandas as pd
import numpy as np
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import os

def simulate_replicates(input_file: str, num_replicates: int, variation_factor: float, output_file: str, random_seed: int):
    """
    Simulates technical replicates for each sample in the input mOTU table.

    Args:
        input_file (str): Path to the input TSV file containing the mOTU table.
        num_replicates (int): The number of replicates to generate for each sample.
        variation_factor (float): The standard deviation for the normal distribution
                                  used to introduce variation. A value of 0.1 means
                                  values will vary by approximately +/- 10%.
        output_file (str): Path to save the new mOTU table with replicates.
        random_seed (int): Seed for the random number generator to ensure reproducibility.
    """
    print(f"🧬 Starting replicate simulation for '{input_file}'...")
    try:
        df = pd.read_csv(input_file, sep='\t', index_col='#OTU_ID')
    except FileNotFoundError:
        print(f"❌ Error: Input file '{input_file}' not found. Please check the path.")
        return
    except Exception as e:
        print(f"❌ Error reading input file: {e}")
        return

    # Identify sample columns dynamically (all columns except the index)
    sample_columns = df.columns.tolist()
    print(f"Identified samples: {sample_columns}")

    # Initialize a reproducible random number generator
    rng = np.random.default_rng(random_seed)

    new_df = df.copy() # Start with a copy of the original data

    for col in sample_columns:
        original_values = df[col].values
        for i in range(1, num_replicates + 1):
            # Generate variation: normal distribution centered at 1 with std dev `variation_factor`
            # Multiplying by this factor ensures variation is proportional to the original value.
            # Clamp values at 0 to avoid negative counts.
            variation = rng.normal(1, variation_factor, size=len(original_values))
            replicate_values = original_values * variation
            replicate_values[replicate_values < 0] = 0 # Ensure no negative counts

            # Round to nearest integer if original data is counts, or keep float if it's relative abundance
            # For simplicity, let's round to integers assuming count data.
            replicate_values = np.round(replicate_values).astype(int)

            new_col_name = f"{col}_rep_{i}"
            new_df[new_col_name] = replicate_values
            print(f"  Generated replicate '{new_col_name}' for sample '{col}'")

    # Save the new mOTU table
    try:
        new_df.to_csv(output_file, sep='\t')
        print(f"✅ New mOTU table with replicates saved to '{output_file}'")
    except Exception as e:
        print(f"❌ Error saving output file: {e}")

def perform_pca(input_file: str, output_plot_file: str):
    """
    Reads the new mOTU table, performs Principal Component Analysis (PCA),
    and creates a scatter plot to visualize the grouping of replicates.

    Args:
        input_file (str): Path to the mOTU table with simulated replicates.
        output_plot_file (str): Path to save the PCA plot (e.g., 'pca_plot.png').
    """
    print(f"\n📊 Performing PCA on '{input_file}'...")
    try:
        df = pd.read_csv(input_file, sep='\t', index_col='#OTU_ID')
    except FileNotFoundError:
        print(f"❌ Error: Input file '{input_file}' not found. Cannot perform PCA.")
        return
    except Exception as e:
        print(f"❌ Error reading input file for PCA: {e}")
        return

    # Transpose the DataFrame so samples are rows and OTUs are columns for PCA
    df_transposed = df.T

    # Handle potential zero variance columns before scaling
    # Remove columns (OTUs) that have zero variance across all samples
    df_transposed = df_transposed.loc[:, (df_transposed != df_transposed.iloc[0]).any()]
    if df_transposed.empty:
        print("⚠️ Warning: No variance found in data after filtering. Cannot perform PCA.")
        return

    # Standardize the data (important for PCA)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_transposed)

    # Perform PCA
    pca = PCA(n_components=2) # Get the first two principal components
    principal_components = pca.fit_transform(scaled_data)

    # Create a DataFrame for plotting
    pca_df = pd.DataFrame(data=principal_components,
                          columns=['Principal Component 1', 'Principal Component 2'],
                          index=df_transposed.index)

    # Extract original sample names for coloring
    # e.g., 'PW_rep_1' -> 'PW', 'IW_rep_2' -> 'IW'
    # Original samples (PW, IW) will also be treated as a group.
    pca_df['SampleGroup'] = pca_df.index.map(lambda x: x.split('_rep_')[0] if '_rep_' in x else x)

    # Plotting
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='Principal Component 1', y='Principal Component 2',
                    hue='SampleGroup', # Color by original sample group
                    data=pca_df,
                    s=100, # Marker size
                    alpha=0.8, # Transparency
                    edgecolor='w', # White edge for markers
                    palette='tab10') # A good palette for distinct colors

    plt.title('PCA of Original Samples and Simulated Replicates', fontsize=16)
    plt.xlabel(f'Principal Component 1 ({pca.explained_variance_ratio_[0]*100:.2f}%)', fontsize=12)
    plt.ylabel(f'Principal Component 2 ({pca.explained_variance_ratio_[1]*100:.2f}%)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Sample Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save the plot
    try:
        plt.savefig(output_plot_file, dpi=300, bbox_inches='tight')
        print(f"✅ PCA plot saved to '{output_plot_file}'")
    except Exception as e:
        print(f"❌ Error saving PCA plot: {e}")
    plt.close() # Close the plot to free memory

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulate technical replicates for microbial samples and perform PCA."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input TSV mOTU table file (e.g., PWIW_motu.tsv)."
    )
    parser.add_argument(
        "--output_data",
        type=str,
        default="simulated_motu_table.tsv",
        help="Path to save the new mOTU table with simulated replicates."
    )
    parser.add_argument(
        "--num_replicates",
        type=int,
        default=3,
        help="Number of technical replicates to generate per sample (default: 3)."
    )
    parser.add_argument(
        "--variation_factor",
        type=float,
        default=0.1,
        help="Degree of variation between replicates (standard deviation for normal distribution, e.g., 0.1 for ~10%% variation). Default: 0.1."
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)."
    )
    parser.add_argument(
        "--output_pca_plot",
        type=str,
        default="pca_replicates_plot.png",
        help="Path to save the PCA plot (e.g., pca_replicates_plot.png)."
    )

    args = parser.parse_args()

    # Run the simulation
    simulate_replicates(
        input_file=args.input,
        num_replicates=args.num_replicates,
        variation_factor=args.variation_factor,
        output_file=args.output_data,
        random_seed=args.random_seed
    )

    # Run the PCA if the output data file was successfully created
    if os.path.exists(args.output_data):
        perform_pca(
            input_file=args.output_data,
            output_plot_file=args.output_pca_plot
        )
    else:
        print("\n🚫 PCA skipped because the simulated data file was not created successfully.")

