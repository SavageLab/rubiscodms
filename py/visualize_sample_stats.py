# == Native Modules
import pandas as pd
import numpy as np
# == Installed Modules
from matplotlib import pyplot as plt
# == Project Modules


def func(pct, allvals):
    absolute = int(np.round(pct/100.*np.sum(allvals)))
    return f"{pct:.1f}%\n{absolute:d}"


def main():
	# === SNAKEMAKE I/O ===
	# == INPUTS
	full_labeled_barcode_path = str(snakemake.input.full_labeled_barcode_path)
	# == OUTPUTS
	sample_pie_chart_path = str(snakemake.output.sample_pie_chart)
	# == WILDCARDS
	experiment = str(snakemake.wildcards.file_prefix)

	#DEBUG
	full_labeled_barcode_path = "/groups/doudna/projects/daniel_projects/rubiscodms/rubisco_nextseq_processing/result_tables/full_labeled_barcodeCounts.csv"

	# Load dataframe
	full_labeled_barcode = pd.read_csv(full_labeled_barcode_path, low_memory=False)

	experiment_data_dict = {}
	experiment_data_df = full_labeled_barcode[['Barcode', experiment, 'BC_Category']].dropna()

	experiment_data_dict[experiment] = {}
	for category in experiment_data_df['BC_Category'].unique():
		experiment_data_dict[experiment][category] = experiment_data_df[experiment_data_df['BC_Category'] ==
																		category][experiment].sum()

	#
	dataForPie = pd.DataFrame(experiment_data_dict)

	#
	fig, ax = plt.subplots(figsize=(10, 5), subplot_kw=dict(aspect="equal"))
	wedges, texts, autotexts = ax.pie(dataForPie[experiment],
									  autopct=lambda pct: func(pct, dataForPie[experiment]),
									  textprops=dict(color="black"),
									  labels=dataForPie.index,
									  startangle=90)
	plt.title(experiment + ' read breakdown')
	print('Verified reads: ' +
		  str(100 * (dataForPie[experiment]['Verified'] / dataForPie[experiment].sum()))[:4] + '%')
	plt.savefig(sample_pie_chart_path,  format="png", dpi=300)
	plt.show()


if __name__ == "__main__":
	main()
