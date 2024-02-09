# == Native Modules ==
# == Installed Modules ==
from matplotlib import pyplot as plt
import pandas as pd
# == Project Modules ==


def main():
	# Snakemake I/O
	# === Inputs
	full_mutation_table = str(snakemake.input.full_mutation_table)
	# === Outputs
	mutation_positions_plot = str(snakemake.output.mutation_positions_plot)
	mutation_positions_negative_plot = str(snakemake.output.mutation_positions_negative_plot)

	# Import barcode tables
	mutation_df = pd.read_csv(full_mutation_table)

	mutation_positions = []
	for row in mutation_df.iloc:
		mutation_list = row.Mutation_List
		for mutation in mutation_list:
			mutation_positions.append(mutation[1])

	# Compute mutation positions data
	mutation_position_data = pd.Series(mutation_positions).dropna().value_counts().sort_index()

	# Set plot size and plot the data
	plt.figure(figsize=(20, 5))
	plt.plot(mutation_position_data.index,
			 mutation_position_data,
			 lw=0.2,
			 color='black',
			 drawstyle='steps')

	# Set plot scale and display the plot
	plt.yscale('log')
	plt.xlabel('Mutation Position', fontsize=12)
	plt.ylabel('Frequency', fontsize=12)
	plt.title('Mutation Position Frequency', fontsize=14)
	plt.grid(True, linestyle='--', alpha=0.5)
	plt.savefig(mutation_positions_plot, format="png", dpi=300)
	plt.close('all')

	# Compute mutation positions data
	mutation_position_data = pd.Series(mutation_positions).dropna()

	# Set plot size and plot the histogram
	plt.figure(figsize=(20, 5))
	plt.hist(mutation_position_data,
			 bins=1600,
			 color='black')

	# Set labels, scale, and display the plot
	plt.xlabel('Nucleotide Position in Plasmid', fontsize=12)
	plt.ylabel('Mutations at That Position', fontsize=12)
	plt.title('Mutation Distribution across Nucleotide Positions', fontsize=14)
	plt.yscale('log')
	plt.grid(True, linestyle='--', alpha=0.5)
	plt.savefig(mutation_positions_negative_plot, format="png", dpi=300)


if __name__ == "__main__":
	main()
