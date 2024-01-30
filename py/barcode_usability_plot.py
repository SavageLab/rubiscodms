# == Native Modules ==
# == Installed Modules ==
import pandas as pd
from matplotlib import pyplot as plt
# == Project Modules ==


def main():
    # Snakemake I/O
    # === Inputs
    barcode_path = str(snakemake.input.output_barcode_path)
    barcode_count_path = str(snakemake.input.output_barcode_count_path)
    # === Outputs
    barcode_plot = str(snakemake.output.barcode_plot)

    # Import barcode tables
    barcode_dataframe = pd.read_csv(barcode_path)
    barcode_counts_df = pd.read_csv(barcode_count_path)
    # Define the threshold for the minimum number of reads per barcode to be considered usable
    read_threshold = 2

    # Calculate the number of barcodes that meet or exceed the usability threshold
    usable_barcodes_count = len(barcode_counts_df[barcode_counts_df['Counts'] >= read_threshold])

    # Set the size of the figure for the histogram
    plt.figure(figsize=(20, 5))

    # Merge the reads DataFrame with the unique barcodes DataFrame on the 'Barcode' column
    # This allows us to associate the count of reads with each barcode
    histogram_data = barcode_dataframe.merge(barcode_counts_df, on='Barcode')

    # Define the range and bin size for the histogram
    bin_range = range(0, 100, 1)

    # Plot a histogram of the number of reads per barcode
    plt.hist(histogram_data['Counts'], bins=bin_range, histtype='step', color='black')

    # Set the y-axis to a logarithmic scale to better visualize the distribution
    plt.yscale('log')

    # Label the y-axis and x-axis
    plt.ylabel('Number of barcodes')
    plt.xlabel('Number of reads per barcode')

    # Add a vertical dashed line at the read threshold to indicate the usability cutoff
    plt.axvline(x=read_threshold, linestyle='--', color='red')

    # Annotate the plot with text indicating the usability threshold and the count of potentially usable barcodes
    plt.text(read_threshold + 1, plt.ylim()[1] / 10, f'Usability threshold: {read_threshold} reads', color='red')
    plt.text(read_threshold + 1, plt.ylim()[1] / 20, f'Potentially usable barcodes: {usable_barcodes_count}', color='red')

    # Display the histogram
    plt.savefig(barcode_plot, format="png", dpi=300)


if __name__ == "__main__":
    main()
