# == Native Modules ==
# == Installed Modules ==
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
# == Project Modules ==


def main():
	# Snakemake I/O
	# === Inputs
	filtered_mutation_table = str(snakemake.input.filtered_mutation_table)
	by_mutation_filtered_table = str(snakemake.input.by_mutation_filtered_table)
	barcode_count_path = str(snakemake.input.barcode_count_path)
	# === Outputs
	mutation_stats_report = str(snakemake.output.mutation_stats_report)
	barcode_per_mut_plot_output = str(snakemake.output.barcode_per_mut_plot_output)
	barcode_scatter_output = str(snakemake.output.barcode_scatter_output)
	mutation_heatmap_output = str(snakemake.output.mutation_heatmap_output)

	# Import barcode tables
	bc_counts_df = pd.read_csv(barcode_count_path)
	filtered_mut_df = pd.read_csv(filtered_mutation_table)
	by_mut_filtered_df = pd.read_csv(by_mutation_filtered_table)

	# Calculate statistics
	usable_reads = bc_counts_df[bc_counts_df['Barcode'].isin(by_mut_filtered_df['Barcode'])]['count'].sum()
	total_reads = bc_counts_df['count'].sum()
	usable_reads_perc = str(100 * usable_reads / total_reads)[:5]

	total_mut_with_bc = filtered_mut_df['RbcLCodonMut'].nunique()
	expected_mutations = 8778
	perc_mut_with_bc = str(100 * total_mut_with_bc / expected_mutations)[:5]

	# Calculate mutations with at least three barcodes
	total_mut_bc_counts = pd.DataFrame(by_mut_filtered_df.groupby('RbcLCodonMut').size()).reset_index().rename(
		columns={"0": "HowManyBarcodes"})
	total_mut_with_3_or_more_bc = total_mut_bc_counts[total_mut_bc_counts['HowManyBarcodes'] >= 3]
	total_mut_with_3_or_more_bc_perc = str(100 * len(total_mut_with_3_or_more_bc) / expected_mutations)[:5]

	# Include the second and last two codons for calculation
	total_mut_with_3_or_more_bc_adjusted_perc = str(
		100 * len(total_mut_with_3_or_more_bc) / (expected_mutations + (19 * 3)))[:5]

	# Write statistics to a report file
	with open(mutation_stats_report, 'w') as report_handle:
		report_handle.write(f"{usable_reads_perc}% of reads should be usable at t0\n")
		report_handle.write(f"{perc_mut_with_bc}% of mutations with at least one barcode\n")
		report_handle.write(f"{total_mut_with_3_or_more_bc_perc}% of mutations with at least three barcodes\n")
		report_handle.write(f"{total_mut_with_3_or_more_bc_adjusted_perc}% of mutations with at least three "
							f"barcodes if we include the second and last two codons\n")

	# === HISTOGRAM ===
	# Histogram of barcodes per mutation
	mutBCCounts = by_mut_filtered_df.groupby('RbcLCodonMut').count()['Barcode']

	binSet = range(0, 50, 1)
	plt.hist(mutBCCounts, bins=binSet, color='lightgray')

	plt.title('Internal replicates', size=30)
	plt.ylabel('Mutations', size=20)
	plt.xlabel('Barcodes per Mutation', size=20)

	# Display the histogram
	plt.savefig(barcode_per_mut_plot_output, format="png", dpi=300)
	plt.close('all')
	# === SCATTER PLOT ===
	# Group mutations by amino acid position and count unique mutations
	pos_vs_muts = pd.DataFrame(by_mut_filtered_df.groupby('AApos').apply(lambda x: x.RbcLCodonMut.nunique())) \
		.reset_index().rename(columns={"RbcLCodonMut": "mutsAtPos"})

	# Count barcodes by position
	mutation_bcs_by_position = by_mut_filtered_df.groupby('AApos')['Barcode'].count()

	# Create scatter plot
	fig, ax = plt.subplots(figsize=(20, 5))

	# Scatter plot for mutations at each position
	ax.scatter(x=pos_vs_muts['AApos'], y=pos_vs_muts['mutsAtPos'], s=3, color='black', marker="_")
	ax.set_xlabel("Amino acid position", fontsize=14)
	ax.set_ylabel("Mutations by amino acid position", fontsize=14)
	ax.set_ylim([0, 20])
	ax.set_xlim([-1, 465])
	ax.set_yticks([0, 2, 5, 10, 15, 18, 19])

	# Create a twin object for the second y-axis
	ax2 = ax.twinx()

	# Scatter plot for barcodes of mutations at each position
	ax2.scatter(mutation_bcs_by_position.index, mutation_bcs_by_position, s=25, color='green', marker='_')
	ax2.set_ylabel("Barcodes of mutations by amino acid position", fontsize=14)
	ax2.set_ylim(0, 1000)

	# Save and display the plot
	plt.savefig(barcode_scatter_output, format="png", dpi=300)
	plt.close('all')

	# === MUTATION FREQUENCY HEATMAP ===
	# Calculate barcode frequency for each individual mutation
	position_aa_freq = by_mut_filtered_df.groupby(['AApos', 'mutAA']).size()
	df_position_aa_freq = position_aa_freq.to_frame(name='Barcodes')
	pivot_bcs = df_position_aa_freq.pivot_table(index='mutAA', columns='AApos', values='Barcodes')

	# Create heatmap plot
	fig, ax1 = plt.subplots(1, tight_layout=True, figsize=(10, 5), dpi=100)

	# Plot heatmap
	pcm_bcs = ax1.pcolor(pivot_bcs, cmap='terrain', norm=matplotlib.colors.LogNorm())

	# Set plot labels and titles
	ax1.set_title('Barcodes per mutant', size=15)
	ax1.set_xlabel('Amino acid position in rbcL', size=15)
	ax1.set_ylabel('Amino acid', size=15)

	# Set ticks in the middle
	ax1.set_yticks(np.arange(0.5, 20.5))
	ax1.set_yticklabels(pivot_bcs.index)

	# Add color bar
	fig.colorbar(pcm_bcs, ax=ax1)

	# Save and display the plot
	plt.savefig(mutation_heatmap_output, format="png", dpi=300)
	plt.close('all')


if __name__ == "__main__":
	main()
