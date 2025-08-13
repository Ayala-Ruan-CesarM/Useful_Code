import pandas as pd
import numpy as np
import argparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import MDS # For PCoA
from sklearn.metrics.pairwise import pairwise_distances # For Bray-Curtis
import matplotlib.pyplot as plt
import seaborn as sns
import os
import umap.umap_ as umap # For UMAP, requires umap-learn

def simulate_replicates(input_file: str, num_replicates: int, distribution: str, variation_factor: float, output_file: str, random_seed: int) -> pd.DataFrame:
    """
    Simulates technical replicates for each sample in the input mOTU table
    using a specified distribution and orders them directly next to the original sample.

    Args:
        input_file (str): Path to the input TSV file containing the mOTU table.
        num_replicates (int): The number of replicates to generate for each sample.
        distribution (str): The statistical distribution to use for generating variation
                            ('normal', 'poisson', 'negative_binomial').
        variation_factor (float): Parameter controlling the degree of variation:
                                  - 'normal': Standard deviation of the normal distribution.
                                  - 'poisson': Not directly used to control extra variance beyond mean.
                                  - 'negative_binomial': Dispersion parameter (alpha), higher value means more dispersion.
        output_file (str): Path to save the new mOTU table with replicates.
        random_seed (int): Seed for the random number generator to ensure reproducibility.

    Returns:
        pd.DataFrame: The DataFrame with simulated replicates, with columns ordered
                      as [Original_Sample1, Rep1, Rep2, ..., Original_Sample2, ...].
    """
    print(f"🧬 Starting replicate simulation for '{input_file}' using '{distribution}' distribution...")
    try:
        df = pd.read_csv(input_file, sep='\t', index_col='#OTU_ID')
    except FileNotFoundError:
        print(f"❌ Error: Input file '{input_file}' not found. Please check the path.")
        return pd.DataFrame()
    except Exception as e:
        print(f"❌ Error reading input file: {e}")
        return pd.DataFrame()

    # Identify sample columns dynamically (all columns except the index)
    sample_columns = df.columns.tolist()
    print(f"Identified samples: {sample_columns}")

    # Initialize a reproducible random number generator
    rng = np.random.default_rng(random_seed)

    # Use a list to build the new DataFrame in the correct order
    ordered_columns_data = {}

    for col in sample_columns:
        # Add the original sample column
        ordered_columns_data[col] = df[col]
        original_values = df[col].values

        # Add the simulated replicates for this sample
        for i in range(1, num_replicates + 1):
            replicate_values = np.zeros_like(original_values, dtype=int)

            if distribution == 'normal':
                # Add variation proportional to the original value
                variation = rng.normal(0, variation_factor, size=len(original_values))
                replicate_values = original_values * (1 + variation)
                replicate_values[replicate_values < 0] = 0 # Clamp at 0
                replicate_values = np.round(replicate_values).astype(int)
            elif distribution == 'poisson':
                # For Poisson, the mean (lambda) is the original value.
                for j, ov in enumerate(original_values):
                    replicate_values[j] = rng.poisson(max(1e-9, ov)) if ov > 0 else 0
            elif distribution == 'negative_binomial':
                # Parameterization: mean = mu, dispersion parameter = alpha (variation_factor)
                for j, ov in enumerate(original_values):
                    if ov == 0:
                        replicate_values[j] = 0
                    else:
                        n_param = max(1, int(1 / variation_factor)) if variation_factor > 0 else 1000000
                        p_param = n_param / (n_param + ov)
                        replicate_values[j] = rng.negative_binomial(n_param, p_param)
            else:
                raise ValueError("Invalid distribution type. Choose 'normal', 'poisson', or 'negative_binomial'.")

            new_col_name = f"{col}_rep_{i}"
            ordered_columns_data[new_col_name] = replicate_values
            print(f"  Generated replicate '{new_col_name}' for sample '{col}'")

    # Create the final DataFrame from the ordered dictionary
    new_df = pd.DataFrame(ordered_columns_data, index=df.index)

    # Save the new mOTU table
    try:
        new_df.to_csv(output_file, sep='\t')
        print(f"✅ New mOTU table with replicates saved to '{output_file}'")
    except Exception as e:
        print(f"❌ Error saving output file: {e}")
    return new_df

def filter_low_abundance_otus(data_df: pd.DataFrame, min_total_counts: int, min_prevalence: int, output_file: str) -> pd.DataFrame:
    """
    Filters low-abundance OTUs from the mOTU table.

    Args:
        data_df (pd.DataFrame): The input mOTU table (OTUs as index, samples as columns).
        min_total_counts (int): Minimum total counts for an OTU across all samples to be kept.
        min_prevalence (int): Minimum number of samples/replicates an OTU must be present in (count > 0) to be kept.
        output_file (str): Path to save the filtered mOTU table.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    print(f"\n🧹 Filtering low-abundance OTUs (min_total_counts={min_total_counts}, min_prevalence={min_prevalence})...")
    filtered_df = data_df.copy()

    # Filter by minimum total counts across all samples/replicates for each OTU
    initial_otus = len(filtered_df)
    filtered_df = filtered_df[filtered_df.sum(axis=1) >= min_total_counts]
    print(f"  {initial_otus - len(filtered_df)} OTUs removed based on total counts.")
    initial_otus = len(filtered_df)

    # Filter by minimum prevalence (number of samples/replicates with count > 0)
    filtered_df = filtered_df[(filtered_df > 0).sum(axis=1) >= min_prevalence]
    print(f"  {initial_otus - len(filtered_df)} OTUs removed based on prevalence.")

    if filtered_df.empty:
        print("⚠️ Warning: DataFrame is empty after filtering. No data for analysis.")
        return pd.DataFrame()

    try:
        filtered_df.to_csv(output_file, sep='\t')
        print(f"✅ Filtered mOTU table saved to '{output_file}'")
    except Exception as e:
        print(f"❌ Error saving filtered file: {e}")
    return filtered_df

def perform_pca(data_df: pd.DataFrame, output_plot_file: str, title_suffix: str = ""):
    """
    Performs Principal Component Analysis (PCA) and creates a scatter plot.

    Args:
        data_df (pd.DataFrame): The mOTU table (OTUs as index, samples as columns).
        output_plot_file (str): Path to save the PCA plot.
        title_suffix (str): Suffix to add to the plot title (e.g., " (Filtered Data)").
    """
    print(f"\n📊 Performing PCA on data{title_suffix}...")
    if data_df.empty:
        print("🚫 Cannot perform PCA: Input DataFrame is empty.")
        return

    # Transpose the DataFrame so samples are rows and OTUs are columns for PCA
    df_transposed = data_df.T

    # Handle potential zero variance columns before scaling
    # Remove columns (OTUs) that have zero variance across all samples
    df_transposed = df_transposed.loc[:, (df_transposed != df_transposed.iloc[0]).any()]
    if df_transposed.empty:
        print("⚠️ Warning: No variance found in data after filtering zero-variance OTUs. Cannot perform PCA.")
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

    plt.title(f'PCA of Samples and Replicates {title_suffix}', fontsize=16)
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

def perform_bray_curtis_analysis(data_df: pd.DataFrame, output_heatmap_file: str, output_pcoa_plot_file: str, random_seed: int, title_suffix: str = ""):
    """
    Calculates Bray-Curtis distances, creates a heatmap, and performs PCoA.

    Args:
        data_df (pd.DataFrame): The mOTU table (OTUs as index, samples as columns).
        output_heatmap_file (str): Path to save the Bray-Curtis heatmap.
        output_pcoa_plot_file (str): Path to save the PCoA plot.
        random_seed (int): Seed for MDS (PCoA) for reproducibility.
        title_suffix (str): Suffix to add to the plot titles.
    """
    print(f"\n📏 Performing Bray-Curtis distance analysis{title_suffix}...")
    if data_df.empty:
        print("🚫 Cannot perform Bray-Curtis analysis: Input DataFrame is empty.")
        return

    # Transpose the DataFrame so samples are rows and OTUs are columns
    df_transposed = data_df.T

    # Calculate Bray-Curtis distance matrix
    # Ensure no negative values for Bray-Curtis
    if (df_transposed < 0).any().any():
        print("⚠️ Warning: Negative values found in data. Bray-Curtis distance requires non-negative values. Clamping to 0.")
        df_transposed[df_transposed < 0] = 0

    # Ensure all values are floats for pairwise_distances
    df_transposed = df_transposed.astype(float)

    # Handle cases where a row (sample) sums to zero, which can cause issues with Bray-Curtis
    # Replace NaN with 0 after division, if any occur due to zero sums
    distance_matrix = pairwise_distances(df_transposed, metric='braycurtis')
    distance_matrix = np.nan_to_num(distance_matrix, nan=0.0) # Replace NaN with 0

    distance_df = pd.DataFrame(distance_matrix, index=df_transposed.index, columns=df_transposed.index)

    # --- Plot Heatmap ---
    plt.figure(figsize=(12, 10))
    sns.heatmap(distance_df, cmap='viridis', annot=False, fmt=".2f", linewidths=.5, cbar_kws={'label': 'Bray-Curtis Distance'})
    plt.title(f'Bray-Curtis Distance Heatmap {title_suffix}', fontsize=16)
    plt.xlabel('Samples/Replicates', fontsize=12)
    plt.ylabel('Samples/Replicates', fontsize=12)
    plt.tight_layout()
    try:
        plt.savefig(output_heatmap_file, dpi=300, bbox_inches='tight')
        print(f"✅ Bray-Curtis heatmap saved to '{output_heatmap_file}'")
    except Exception as e:
        print(f"❌ Error saving Bray-Curtis heatmap: {e}")
    plt.close()

    # --- Perform PCoA (MDS) on Bray-Curtis distance matrix ---
    # n_components=2 for 2D plot
    # dissimilarity='precomputed' because we are providing a distance matrix
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=random_seed, normalized_stress=False)
    try:
        pcoa_components = mds.fit_transform(distance_matrix)
    except ValueError as e:
        print(f"❌ Error performing PCoA: {e}. This might happen if the distance matrix contains NaNs or infinities, or if there's insufficient variance.")
        return

    pcoa_df = pd.DataFrame(data=pcoa_components,
                           columns=['PCoA Component 1', 'PCoA Component 2'],
                           index=df_transposed.index)

    pcoa_df['SampleGroup'] = pcoa_df.index.map(lambda x: x.split('_rep_')[0] if '_rep_' in x else x)

    # Plotting PCoA
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PCoA Component 1', y='PCoA Component 2',
                    hue='SampleGroup',
                    data=pcoa_df,
                    s=100,
                    alpha=0.8,
                    edgecolor='w',
                    palette='tab10')

    plt.title(f'PCoA of Samples and Replicates (Bray-Curtis) {title_suffix}', fontsize=16)
    plt.xlabel('PCoA Component 1', fontsize=12)
    plt.ylabel('PCoA Component 2', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Sample Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    try:
        plt.savefig(output_pcoa_plot_file, dpi=300, bbox_inches='tight')
        print(f"✅ Bray-Curtis PCoA plot saved to '{output_pcoa_plot_file}'")
    except Exception as e:
        print(f"❌ Error saving Bray-Curtis PCoA plot: {e}")
    plt.close()

def plot_umap(data_df: pd.DataFrame, output_plot_file: str, random_seed: int, title_suffix: str = ""):
    """
    Performs UMAP dimensionality reduction and creates a scatter plot.

    Args:
        data_df (pd.DataFrame): The mOTU table (OTUs as index, samples as columns).
        output_plot_file (str): Path to save the UMAP plot.
        random_seed (int): Seed for UMAP for reproducibility.
        title_suffix (str): Suffix to add to the plot title.
    """
    print(f"\n🌌 Performing UMAP on data{title_suffix}...")
    if data_df.empty:
        print("🚫 Cannot perform UMAP: Input DataFrame is empty.")
        return

    # Transpose the DataFrame so samples are rows and OTUs are columns for UMAP
    df_transposed = data_df.T

    # Handle potential zero variance columns before scaling
    df_transposed = df_transposed.loc[:, (df_transposed != df_transposed.iloc[0]).any()]
    if df_transposed.empty:
        print("⚠️ Warning: No variance found in data after filtering zero-variance OTUs. Cannot perform UMAP.")
        return

    # Standardize the data (important for UMAP)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_transposed)

    # Perform UMAP
    reducer = umap.UMAP(random_state=random_seed, n_components=2)
    try:
        embedding = reducer.fit_transform(scaled_data)
    except Exception as e:
        print(f"❌ Error performing UMAP: {e}. This might happen if the data is too sparse or has issues.")
        return

    umap_df = pd.DataFrame(data=embedding,
                           columns=['UMAP Component 1', 'UMAP Component 2'],
                           index=df_transposed.index)

    umap_df['SampleGroup'] = umap_df.index.map(lambda x: x.split('_rep_')[0] if '_rep_' in x else x)

    # Plotting
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='UMAP Component 1', y='UMAP Component 2',
                    hue='SampleGroup',
                    data=umap_df,
                    s=100,
                    alpha=0.8,
                    edgecolor='w',
                    palette='tab10')

    plt.title(f'UMAP of Samples and Replicates {title_suffix}', fontsize=16)
    plt.xlabel('UMAP Component 1', fontsize=12)
    plt.ylabel('UMAP Component 2', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Sample Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    try:
        plt.savefig(output_plot_file, dpi=300, bbox_inches='tight')
        print(f"✅ UMAP plot saved to '{output_plot_file}'")
    except Exception as e:
        print(f"❌ Error saving UMAP plot: {e}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulate technical replicates for microbial samples, filter low-abundance OTUs, and perform PCA, PCoA (Bray-Curtis), and UMAP analysis."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input TSV mOTU table file (e.g., PWIW_motu.tsv)."
    )
    parser.add_argument(
        "--output_simulated_data",
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
        "--distribution",
        type=str,
        choices=['normal', 'poisson', 'negative_binomial'],
        default='normal',
        help="Statistical distribution for generating variation ('normal', 'poisson', 'negative_binomial'). Default: 'normal'."
    )
    parser.add_argument(
        "--variation_factor",
        type=float,
        default=0.1,
        help="Degree of variation parameter (float): "
             "  - 'normal': Standard deviation (e.g., 0.1 for ~10%% variation).\n"
             "  - 'poisson': Not applicable (variance tied to mean).\n"
             "  - 'negative_binomial': Dispersion parameter (alpha), higher value means more dispersion. Default: 0.1."
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)."
    )
    parser.add_argument(
        "--output_pca_plot_unfiltered",
        type=str,
        default="pca_replicates_unfiltered.png",
        help="Path to save the PCA plot for unfiltered data."
    )
    parser.add_argument(
        "--output_heatmap_unfiltered",
        type=str,
        default="bray_curtis_heatmap_unfiltered.png",
        help="Path to save the Bray-Curtis heatmap for unfiltered data."
    )
    parser.add_argument(
        "--output_pcoa_plot_unfiltered",
        type=str,
        default="bray_curtis_pcoa_unfiltered.png",
        help="Path to save the Bray-Curtis PCoA plot for unfiltered data."
    )
    parser.add_argument(
        "--output_umap_plot_unfiltered",
        type=str,
        default="umap_replicates_unfiltered.png",
        help="Path to save the UMAP plot for unfiltered data."
    )
    parser.add_argument(
        "--filter_data",
        action='store_true',
        help="Set this flag to enable low-abundance OTU filtering."
    )
    parser.add_argument(
        "--min_total_counts",
        type=int,
        default=10,
        help="Minimum total counts for an OTU across all samples/replicates to be kept (used with --filter_data). Default: 10."
    )
    parser.add_argument(
        "--min_prevalence",
        type=int,
        default=2,
        help="Minimum number of samples/replicates an OTU must be present in (count > 0) to be kept (used with --filter_data). Default: 2."
    )
    parser.add_argument(
        "--output_filtered_data",
        type=str,
        default="filtered_motu_table.tsv",
        help="Path to save the filtered mOTU table (used with --filter_data)."
    )
    parser.add_argument(
        "--output_pca_plot_filtered",
        type=str,
        default="pca_replicates_filtered.png",
        help="Path to save the PCA plot for filtered data (used with --filter_data)."
    )
    parser.add_argument(
        "--output_heatmap_filtered",
        type=str,
        default="bray_curtis_heatmap_filtered.png",
        help="Path to save the Bray-Curtis heatmap for filtered data (used with --filter_data)."
    )
    parser.add_argument(
        "--output_pcoa_plot_filtered",
        type=str,
        default="bray_curtis_pcoa_filtered.png",
        help="Path to save the Bray-Curtis PCoA plot for filtered data (used with --filter_data)."
    )
    parser.add_argument(
        "--output_umap_plot_filtered",
        type=str,
        default="umap_replicates_filtered.png",
        help="Path to save the UMAP plot for filtered data (used with --filter_data)."
    )

    args = parser.parse_args()

    # --- Step 1: Simulate Replicates ---
    simulated_df = simulate_replicates(
        input_file=args.input,
        num_replicates=args.num_replicates,
        distribution=args.distribution,
        variation_factor=args.variation_factor,
        output_file=args.output_simulated_data,
        random_seed=args.random_seed
    )

    if simulated_df.empty:
        print("Exiting due to empty simulated data.")
        exit()

    # --- Step 2: Perform analyses on unfiltered simulated data ---
    print("\n--- Analyzing Unfiltered Simulated Data ---")
    perform_pca(
        data_df=simulated_df,
        output_plot_file=args.output_pca_plot_unfiltered,
        title_suffix="(Unfiltered)"
    )
    perform_bray_curtis_analysis(
        data_df=simulated_df,
        output_heatmap_file=args.output_heatmap_unfiltered,
        output_pcoa_plot_file=args.output_pcoa_plot_unfiltered,
        random_seed=args.random_seed,
        title_suffix="(Unfiltered)"
    )
    plot_umap(
        data_df=simulated_df,
        output_plot_file=args.output_umap_plot_unfiltered,
        random_seed=args.random_seed,
        title_suffix="(Unfiltered)"
    )

    # --- Step 3: Filter Low-Abundance OTUs if requested ---
    if args.filter_data:
        filtered_df = filter_low_abundance_otus(
            data_df=simulated_df,
            min_total_counts=args.min_total_counts,
            min_prevalence=args.min_prevalence,
            output_file=args.output_filtered_data
        )

        if filtered_df.empty:
            print("Exiting analysis of filtered data due to empty DataFrame after filtering.")
        else:
            # --- Step 4: Perform analyses on filtered data ---
            print("\n--- Analyzing Filtered Simulated Data ---")
            perform_pca(
                data_df=filtered_df,
                output_plot_file=args.output_pca_plot_filtered,
                title_suffix="(Filtered)"
            )
            perform_bray_curtis_analysis(
                data_df=filtered_df,
                output_heatmap_file=args.output_heatmap_filtered,
                output_pcoa_plot_file=args.output_pcoa_plot_filtered,
                random_seed=args.random_seed,
                title_suffix="(Filtered)"
            )
            plot_umap(
                data_df=filtered_df,
                output_plot_file=args.output_umap_plot_filtered,
                random_seed=args.random_seed,
                title_suffix="(Filtered)"
            )
