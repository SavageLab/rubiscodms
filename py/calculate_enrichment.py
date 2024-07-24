# == Native Modules

# == Installed Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, ttest_ind
import itertools
# == Project Modules


def get_t0(t0Lookup, experiment):
	for t0s, experiments in t0Lookup.items():
		if experiment in experiments:
			return t0s
	return None


def enrichment_df(conditions_key, t0_lookup, df0, pseudo_const0, threshold0,
				  include_pval=False,
				  dead_muts = ['K191', 'K166', 'K329', 'D193', 'E194', 'H287'],
				  indexHop_threshold = 20):
	# define names of experimental columns
	expnames = [col for col in df0.columns if ('NP' in col)
				and (col not in conditions_key['t0'] + conditions_key['WT_only'])]

	# get t0s corresponding to expnames
	t0name_match = [get_t0(t0_lookup, en) for en in expnames]

	# make dataframes for the selection and t0 conditions and replace NaNs with zeros
	tF_reads = df0[expnames].fillna(0)
	t0_reads = df0[t0name_match].fillna(0)

	# count the total number of reads across each experiment
	# print(f"Check data before SUM: {df0[t0name_match]}")
	t0read_nums = df0[t0name_match].sum(axis=0)
	expread_nums = df0[expnames].sum(axis=0)

	# get the indices of bad wild-types due to index hopping
	badWT = ((df0['RbcLCodonMut'] == 'WT') & (df0[conditions_key['WT_only'][0]] >= indexHop_threshold))

	# create logical dataframe to exlude the bad wild-types and any barcode for which neither tF or t0 exceed the threshold
	over_threshold = np.logical_and(np.array(~badWT).reshape(-1, 1), np.logical_or(df0[expnames] > threshold0,
																				   np.array(
																					   df0[t0name_match] > threshold0)))

	# compute enrichment values by adding the pseudocount constant multiplied by the number of reads and also normalizing
	# for the read count for each condition.
	enrich_barcode = (np.log10(tF_reads.add(np.array(expread_nums * pseudo_const0))
							   / np.array(t0_reads.add(np.array(t0read_nums * pseudo_const0))))
					  + np.log10(np.array(t0read_nums) / expread_nums)).where(over_threshold)

	# append the mutant label column to the barcode enrichment dataframe
	enrich_barcode['mutant'] = df0['RbcLCodonMut']

	# compute the median enrichment for each variant
	enrich_variant = enrich_barcode.groupby(by='mutant').median()

	# compute the median enrichment for dead mutants by variant
	dead_enrich = enrich_variant.loc[enrich_variant.index.str.startswith(tuple(dead_muts))].median()

	# compute the wild-type enrichment (note that the bad wild-types have been removed by this stage)
	wt_enrich = enrich_variant.loc['WT']

	# compute the normalized enrichment using something akin to a Z-score
	# 0 is dead, 1 is WT
	enrich_norm = (enrich_variant - dead_enrich) / (wt_enrich - dead_enrich)

	if include_pval:
		def ttest_func(x0, wt0):
			# compute two-tailed t-test between wt and mutant while omitting NaNs
			_, pval0 = ttest_ind(wt0, x0, nan_policy='omit')

			# return the geometric mean of the pvalues
			return np.exp(np.mean(np.log(pval0.data)))

		pvalcols = conditions_key['5percent_CO2_20uM_IPTG']
		wt_values = enrich_barcode[enrich_barcode['mutant'] == 'WT']

		# generate a series with the pvalues
		# ttest_pvals = enrich_barcode.groupby(by='mutant').apply(lambda x: ttest_func(x[pvalcols], wt_values[pvalcols]))
		ttest_pvals = enrich_barcode.groupby(by='mutant').apply(
			lambda x: ttest_func(x[pvalcols].mean(axis=1), wt_values[pvalcols].mean(axis=1)))
		enrich_norm['pvalues'] = ttest_pvals

	# return the variant enrichment dataframe.
	return enrich_norm


def rep_pearson_corrs(df0):
	"""
	# Compute the mean pairwise Pearson coefficient for all possible replicate combinations
	:param df0: dataframe with enrichment values across replicates (could be at the barcode or variant level)
	:return: numpy vector of Pearson correlations (n choose 2)
	"""

	# get a list of all pairwise combinations
	pairwise_reps = list(itertools.combinations(df0.columns, 2))
	# create a logical dataframe for determining valid points (pearsonr can't tolerate NaNs)
	nan_df = pd.isna(df0)

	# initialize vector for Pearson corrs.
	pearcorrs = np.zeros(len(pairwise_reps))

	# iterate over all replicate pairs
	for pind, pr in enumerate(pairwise_reps):
		# get valid indices for which neither replicate contains a NaN
		pair_ind = ~nan_df[pr[0]] & ~nan_df[pr[1]]

		# filter the data to ensure no NaNs are passed to pearsonr
		filtered_x = df0[pr[0]][pair_ind]
		filtered_y = df0[pr[1]][pair_ind]

		# Debugging statements to check the data
		print(f"Processing pair: {pr}")
		print(f"Filtered data for {pr[0]}: {filtered_x}")
		print(f"Filtered data for {pr[1]}: {filtered_y}")
		print(f"Are there NaNs in filtered_x? {np.isnan(filtered_x).any()}")
		print(f"Are there NaNs in filtered_y? {np.isnan(filtered_y).any()}")

		# Ensure no NaNs are in the filtered data
		if not np.isnan(filtered_x).any() and not np.isnan(filtered_y).any():
			# compute Pearson correlation
			pearcorrs[pind] = pearsonr(filtered_x, filtered_y)[0]
		else:
			pearcorrs[pind] = np.nan  # Assign NaN if invalid data is found

	return pearcorrs.mean()


def main():
	# Snakemake I/O
	# === Inputs
	# ==
	labeled_barcode_path = str(snakemake.input.labeled_barcode_path)
	# === Outputs
	parameter_sweep_plot = str(snakemake.output.parameter_sweep_plot)
	parameter_sweep_table = str(snakemake.output.parameter_sweep_table)

	# === Constant Variables ===

	# Information about samples
	CONDITIONS_KEY = {
		'WT_only': ['NP_11_64_2'],
		't0': ['NP_11_64_9', 'NP_11_64_11', 'NP_11_64_12', 'NP_11_64_13', 'NP_11_64_15'],
		'0.04percent_CO2': ['NP_11_66_41'],
		'0.2percent_CO2': ['NP_11_66_1', 'NP_11_66_2', 'NP_11_66_3'],
		'0.3percent_CO2': ['NP_11_66_7', 'NP_11_66_8', 'NP_11_66_9'],
		'0.4percent_CO2': ['NP_11_66_4', 'NP_11_66_5', 'NP_11_66_6'],
		'1percent_CO2': ['NP_11_66_10', 'NP_11_66_11', 'NP_11_66_12'],
		'5percent_CO2_20uM_IPTG': ['NP_11_66_13', 'NP_11_66_14', 'NP_11_66_15',
								   'NP_11_66_16', 'NP_11_66_17', 'NP_11_66_18',
								   'NP_11_66_37', 'NP_11_66_38', 'NP_11_66_39'],
		'0uM_IPTG': ['NP_11_66_19', 'NP_11_66_20', 'NP_11_66_21'],
		'5uM_IPTG': ['NP_11_66_22', 'NP_11_66_23', 'NP_11_66_24'],
		'30uM_IPTG': ['NP_11_66_25', 'NP_11_66_26', 'NP_11_66_27'],
		'100uM_IPTG': ['NP_11_66_28', 'NP_11_66_29', 'NP_11_66_30'],
		'300uM_IPTG': ['NP_11_66_31', 'NP_11_66_32', 'NP_11_66_33'],
		'1000uM_IPTG': ['NP_11_66_34', 'NP_11_66_35', 'NP_11_66_36'],
		'Xylose': ['NP_11_66_40']

	}
	T0_LOOKUP = {'NP_11_64_11': ['NP_11_66_1', 'NP_11_66_4', 'NP_11_66_7', 'NP_11_66_10', 'NP_11_66_13',
								 'NP_11_66_16', 'NP_11_66_19', 'NP_11_66_22', 'NP_11_66_25', 'NP_11_66_28',
								 'NP_11_66_31', 'NP_11_66_34', 'NP_11_66_37', 'NP_11_66_40', 'NP_11_66_41'],
				 'NP_11_64_12': ['NP_11_66_2', 'NP_11_66_5', 'NP_11_66_8', 'NP_11_66_11', 'NP_11_66_14',
								 'NP_11_66_17', 'NP_11_66_20', 'NP_11_66_23', 'NP_11_66_26', 'NP_11_66_29',
								 'NP_11_66_32', 'NP_11_66_35', 'NP_11_66_38'],
				 'NP_11_64_13': ['NP_11_66_3', 'NP_11_66_6', 'NP_11_66_9', 'NP_11_66_12', 'NP_11_66_15',
								 'NP_11_66_18', 'NP_11_66_21', 'NP_11_66_24', 'NP_11_66_27', 'NP_11_66_30',
								 'NP_11_66_33', 'NP_11_66_36', 'NP_11_66_39']}

	df_verified_barcodes_trim = pd.read_csv(labeled_barcode_path, index_col=0, low_memory=False)

	# Make list of replicate sets
	rep_names_list = []
	rep_keynames = []
	for key in list(CONDITIONS_KEY.keys()):
		if (key != 't0') & (len(CONDITIONS_KEY[key]) > 1):
			rep_names_list.append(CONDITIONS_KEY[key])
			rep_keynames.append(key)

	# Parameter sweep to find optimal set of
	thresholds = np.arange(0, 20)
	pcounts = np.logspace(-8, -4, 20)

	pearson_mat = np.zeros((len(pcounts), len(thresholds)))
	for tind, th in enumerate(thresholds):
		for pcind, pc in enumerate(pcounts):
			enrich = enrichment_df(CONDITIONS_KEY, T0_LOOKUP, df_verified_barcodes_trim, pc, th, include_pval=False)
			pearson_mat[pcind, tind] = np.mean([rep_pearson_corrs(enrich[reps]) for reps in rep_names_list])

	# Plot the result of the parameter sweep
	thgrid, pcgrid = np.meshgrid(thresholds, pcounts)

	max_index = np.unravel_index(np.argmax(pearson_mat), pearson_mat.shape)

	plt.contourf(pcgrid, thgrid, pearson_mat, 12)
	cbar = plt.colorbar()
	cbar.set_label('mean replicate Pearson correlation', rotation=270, labelpad=20)
	plt.scatter([pcounts[max_index[0]]], [thresholds[max_index[1]]], color='red', marker='x', s=100)
	plt.xscale('log')
	plt.xlabel('pseudocount constant')
	plt.ylabel('minimum count threshold')
	plt.savefig(parameter_sweep_plot, format="png", dpi=300)

	# == Define an array of values for each parameter
	pseudoconsts = np.logspace(-9, -6, 10)
	thresholds = np.linspace(0, 50, 11).astype(int)

	# == Create a meshgrid of parameter values
	pcmat, thmat = np.meshgrid(pseudoconsts, thresholds)

	# == Create an empty list to store the sub-dataframes
	df_list = []

	# == Iterate over the meshgrid using np.nditer
	for pc, th in np.nditer([pcmat, thmat]):
		# print(pc, th)
		#     # run your function with the corresponding parameter values and append the resulting dataframe to the list
		df_list.append(enrichment_df(CONDITIONS_KEY, T0_LOOKUP, df_verified_barcodes_trim, pc, th, include_pval=False))

	# == Create a multi-index from the parameter values
	index = pd.MultiIndex.from_product([pseudoconsts, thresholds], names=['pseudoconst', 'threshold'])

	# == Concatenate the sub-dataframes into a single dataframe with the multi-index
	enrich_combine = pd.concat(df_list, axis=0, keys=index)

	enrich_combine.to_csv(parameter_sweep_table)


if __name__ == "__main__":
	main()
