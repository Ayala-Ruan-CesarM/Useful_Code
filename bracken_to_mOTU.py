import argparse
import pandas as pd
import os

def bracken_to_motu(bracken_files, output_file):
    """
    Converts multiple Bracken output files into a single mOTU table.

    Args:
        bracken_files (list): A list of paths to the Bracken output files.
        output_file (str): The path for the output mOTU table (TSV format).
    """
    # A dictionary to hold the data for each sample
    all_samples_data = {}

    # Process each Bracken file
    for file_path in bracken_files:
        try:
            # Extract the sample name from the file name.
            # This assumes a format like 'sample_name.bracken'.
            # You might need to adjust this based on your file naming convention.
            sample_name = os.path.basename(file_path).split('.')[0]

            # Read the Bracken file into a pandas DataFrame
            # The file is tab-separated and has a header.
            df = pd.read_csv(file_path, sep='\t')

            # We are interested in species-level data.
            # Bracken output can contain other taxonomic ranks, so we filter for 'S'.
            #species_df = df[df['taxonomy_level'] == 'S']

            # We need the species name and the estimated read counts.
            # We set the species name as the index to facilitate merging later.
            sample_data = df[['name', 'new_est_reads']].set_index('name')['new_est_reads']

            # Add the data for the current sample to our dictionary
            all_samples_data[sample_name] = sample_data

        except FileNotFoundError:
            print(f"Error: The file {file_path} was not found.")
            continue
        except Exception as e:
            print(f"An error occurred while processing {file_path}: {e}")
            continue

    # Create a single DataFrame from the dictionary of all samples
    # The keys of the dictionary become the columns of the DataFrame.
    motu_table = pd.DataFrame(all_samples_data)

    # The merging process will result in NaN for species not present in a sample.
    # We replace these NaN values with 0, as it means zero reads were detected.
    motu_table = motu_table.fillna(0).astype(int)

    # The index of the DataFrame is the species name. We'll give it a proper name.
    motu_table.index.name = '#OTU_ID'

    # Save the final mOTU table to the specified output file.
    # The output will be a tab-separated file.
    try:
        motu_table.to_csv(output_file, sep='\t')
        print(f"Successfully created mOTU table at: {output_file}")
    except Exception as e:
        print(f"An error occurred while writing the output file: {e}")

if __name__ == '__main__':
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser(description="Convert Bracken output files to a unified mOTU table.")
    parser.add_argument(
        '-i', '--input',
        nargs='+',  # This allows for one or more input files
        required=True,
        help="A list of Bracken output files (space-separated)."
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="The name of the output mOTU table file (e.g., 'motu_table.tsv')."
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Call the main function with the provided arguments
    bracken_to_motu(args.input, args.output)

    # Example of how to run the script from the command line:
    # python your_script_name.py -i sample1.bracken sample2.bracken sample3.bracken -o combined_motu_table.tsv
