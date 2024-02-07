# == Native Modules ==
# == Installed Modules ==
import pandas as pd
# == Project Modules ==


def main():
	# Snakemake I/O
	# === Inputs
	filtered_mutation_table = str(snakemake.input.filtered_mutation_table)
	by_mutation_filtered_table = str(snakemake.input.by_mutation_filtered_table)
	barcode_count_path = str(snakemake.input.barcode_count_path)
	# === Outputs
	mutation_stats_report = str(snakemake.output.mutation_stats_report)

	# Import barcode tables
	barcode_counts_df = pd.read_csv(barcode_count_path)
	filtered_mutation_df = pd.read_csv(filtered_mutation_table)
	by_mutation_filtered_df = pd.read_csv(by_mutation_filtered_table)

	# Calculate statistics
	usable_reads = barcode_counts_df[barcode_counts_df['Barcode'].isin(by_mutation_filtered_df['Barcode'])][
		'count'].sum()
	total_reads = barcode_counts_df['count'].sum()
	usable_reads_perc = str(100 * usable_reads / total_reads)[:5]

	total_mutations_with_barcode = filtered_mutation_df['RbcLCodonMut'].nunique()
	expected_mutations_df = 8778
	perc_mutations_with_barcode = str(100 * total_mutations_with_barcode / expected_mutations_df)[:5]

	# Calculate mutations with at least three barcodes
	totalMutsBCCounts = pd.DataFrame(by_mutation_filtered_df.groupby('RbcLCodonMut').size()).reset_index().rename(
		columns={"0": "HowManyBarcodes"})
	totalMutsWithThreeOrMoreBC = totalMutsBCCounts[totalMutsBCCounts['HowManyBarcodes'] >= 3]
	totalMutsWithThreeOrMoreBCPercentage = str(100 * len(totalMutsWithThreeOrMoreBC) / expected_mutations_df)[:5]

	# Include the second and last two codons for calculation
	totalMutsWithThreeOrMoreBCAdjustedPercentage = str(
		100 * len(totalMutsWithThreeOrMoreBC) / (expected_mutations_df + (19 * 3)))[:5]

	with open(mutation_stats_report, 'w') as report_handle:
		report_handle.write(f"{usable_reads_perc}% of reads should be usable at t0")
		report_handle.write(f"{perc_mutations_with_barcode}% of mutations with at least one barcode")
		report_handle.write(f"{totalMutsWithThreeOrMoreBCPercentage}% of mutations with at least three barcodes")
		report_handle.write(f"{totalMutsWithThreeOrMoreBCAdjustedPercentage}% of mutations with at least three "
							f"barcodes if we include the second and last two codons")


if __name__ == "__main__":
	main()
