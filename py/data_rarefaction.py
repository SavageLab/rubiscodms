# == Native Modules ==
# == Installed Modules ==
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
# == Project Modules ==


def negativeBinomial(x, a, b):
	return a * (1 - np.power(b, x / a))  # a*(1-b^(x/a))


def main():
	# Snakemake I/O
	# === Inputs
	barcode_count_path = str(snakemake.input.barcode_count_path)
	by_mutation_filtered_table = str(snakemake.input.by_mutation_filtered_table)
	# === Outputs
	rarefaction_plot_output = str(snakemake.output.rarefaction_plot_output)

	# Import barcode tables
	bc_counts_df = pd.read_csv(barcode_count_path)
	by_mut_filtered_df = pd.read_csv(by_mutation_filtered_table)

	# Rarefaction curves data initialization
	rarefaction_data = []

	# Sample sizes for rarefaction
	sample_sizes = [10, 30, 100, 300, 1000, 3000, 10000, 30000]
	how_often_to_sample = 100000

	# Add index to the readsBCDF DataFrame
	bc_counts_df['index'] = np.arange(len(bc_counts_df))

	# Filter readsBCDF based on Barcode presence in filtered_df
	filtered_BC_barcode = bc_counts_df['Barcode'].isin(by_mut_filtered_df['Barcode'])
	filteredBCReads = bc_counts_df[filtered_BC_barcode]

	# Iterate through filteredBCReads to compute rarefaction data
	for read in filteredBCReads['index']:
		num = read % how_often_to_sample
		if num == 0 or read in sample_sizes:
			filtered_BCs_up_to_read = filteredBCReads['Barcode'][:read]
			rarefaction_data.append([filtered_BCs_up_to_read.count(), filtered_BCs_up_to_read.nunique()])

	# Transpose rarefaction data
	rarefaction_data_transposed = list(map(list, zip(*rarefaction_data)))

	# Filtered rarefaction data initialization
	rarefaction_data_filtered = []

	# Add index to the filteredBCReads DataFrame
	filteredBCReads['index'] = np.arange(len(filteredBCReads))

	# Iterate through filteredBCReads to compute filtered rarefaction data
	for read in filteredBCReads['index']:
		num = read % how_often_to_sample
		if num == 0 or read in sample_sizes:
			filtered_BCs_up_to_read = filteredBCReads['Barcode'][:read]
			rarefaction_data_filtered.append([filtered_BCs_up_to_read.count(), filtered_BCs_up_to_read.nunique()])

	# Transpose filtered rarefaction data
	rarefaction_filtered_data_transposed = list(map(list, zip(*rarefaction_data_filtered)))

	# Plot rarefaction

	plt.rcParams["figure.figsize"] = (10, 10)

	# Fit to negative binomial distribution curve
	xdata = rarefaction_filtered_data_transposed[0]
	ydata = rarefaction_filtered_data_transposed[1]

	popt, pcov = curve_fit(f=negativeBinomial,
						   xdata=xdata,
						   ydata=ydata,
						   p0=[10000, 0.5],
						   bounds=(0, np.inf)
						   )

	estText = "Negative binomial \ncomplexity estimate: \n" + str(int(popt[0])) + ' barcodes'
	print("Negative binomial skew estimate: " + str(popt[1])[:4])

	# Plot the curve above
	plt.plot(xdata, negativeBinomial(xdata, *popt), color='black', zorder=0)
	plt.tight_layout()
	plt.text(1300000, 150000, estText, fontsize=12, horizontalalignment='center', verticalalignment='bottom')

	plt.scatter(rarefaction_data_transposed[0],
				rarefaction_data_transposed[1],
				color='black',
				label='All Barcodes',
				s=20)
	plt.scatter(rarefaction_filtered_data_transposed[0],
				rarefaction_filtered_data_transposed[1],
				color='teal',
				label='Filtered Barcodes',
				s=100)

	plt.legend(loc='upper left')

	plt.title('Barcode subsampling rarefaction curve', size=20)
	plt.xlabel('Barcodes in subsample', size=15)
	plt.ylabel('Unique barcodes in subsample', size=15)

	plt.savefig(rarefaction_plot_output, format="png", dpi=300)
	plt.close('all')


if __name__ == "__main__":
	main()
