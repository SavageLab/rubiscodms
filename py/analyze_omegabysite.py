# == Native Modules ==
import natsort
import warnings
import math
# == Installed Modules ==
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
# == Project Modules ==

warnings.filterwarnings('ignore')


def bh_control(pvals, fdr):
	"""Implements Benjamini-Hochberg procedure to control false discovery rate.

	Calling arguments:

	*pvals* : a list of tuples of *(label, p)*  where *label* is some label assigned
	to each data point, and *p* is the corresponding *P-value*.

	*fdr* : the desired false discovery rate

	The return value is the 2-tuple *(pcutoff, significantlabels)*. After applying
	the algorithm, all data points with *p <= pcutoff* are declared significant.
	The labels for these data points are in *significantlabels*. If there are no
	significant sites, *pcutoff* is returned as the maximum P-value that would
	have made a single point significant.
	"""
	num_tests = len(pvals)

	# sort by p-value
	sorted_tests = sorted(pvals, key=lambda tup: tup[1])

	# find maximum rank for which p <= (rank/num_tests)*FDR
	max_rank = 0
	pcutoff = None
	for (rank, (label, p)) in enumerate(sorted_tests):
		rank = rank + 1  # rank beginning with 1 for smallest p-value (there is no rank 0)
		bh_threshold = fdr * float(rank) / num_tests
		if p <= bh_threshold:
			assert rank > max_rank
			max_rank = rank
			pcutoff = bh_threshold

	# pcutoff to have one significant site if there are none
	if pcutoff is None:
		pcutoff = 1.0 / num_tests * fdr

	# collect significant ranks:
	significantlabels = []
	for (rank, (label, p)) in enumerate(sorted_tests):
		rank = rank + 1  # rank beginning with 1 for site with smallest p-vaalue
		if rank <= max_rank:
			assert p <= pcutoff
			significantlabels.append(label)

	return pcutoff, significantlabels


def get_pvalues(omegas):
	# omegas = pandas.read_csv(df, comment='#', sep='\t')
	sites = natsort.realsorted(map(str, omegas.site.values))
	assert set(sites) == set(map(str, omegas['site'].values)), \
		("sites in {0} don't match those in prefs or diffprefs"
		 .format(args['omegabysite']))
	shortname = '$\omega_r$'
	longname = ('{0} $<1 \; \longleftarrow$ $\log_{{10}} P$ '
				'$\longrightarrow \;$ {0} $>1$'.format(shortname))
	prop_d = {}
	for r in sites:
		omegar = float(omegas[omegas['site'].astype(str) == r]['omega'])
		p = float(omegas[omegas['site'].astype(str) == r]['P'])
		if omegar < 1:
			prop_d[r] = max(math.log10(1e-4), math.log10(p))
		else:
			prop_d[r] = -max(math.log10(1e-4), math.log10(p))
	return prop_d.values()


def transform_p(x, p, maxval):
	"""If *x* >= 1, returns *min(-log10(p), maxval)*; if *x < 1* returns *max(-maxval, log10(p))*.

	This transforms *P* values so large values indicate strong evidence for *x > 1*,
	and small (negative) values indicate strong evidence for *x < 1*.
	"""
	assert maxval > 0
	assert 0 <= p <= 1
	if x >= 1:
		if p == 0:
			return maxval
		else:
			return min(-math.log10(p), maxval)
	else:
		if p == 0:
			return -maxval
		else:
			return max(-maxval, math.log10(p))


def main():
	# Snakemake I/O
	# === Inputs
	phydms_root_path = str(snakemake.params.phydms_root_path)
	# === Outputs
	models_summary_plot = str(snakemake.output.models_summary_plot)

	# Use phydms standard names to fetch the 'omegabysite' files
	model_list = [('expcm', f"{phydms_root_path}_ExpCM_aa_preference_omegabysite.txt"),
				  ('yngkp_m5', f"{phydms_root_path}_YNGKP_M5_omegabysite.txt"),
				  ('avg_expcm', f"{phydms_root_path}_averaged_ExpCM_aa_preference_omegabysite.txt"),
				  ('yngkp_m0', f"{phydms_root_path}_YNGKP_M0_omegabysite.txt")
				  ]

	fig, ax = plt.subplots(1, len(model_list), figsize=(20, 5), sharey=True)
	symbols = ['o', 's', '^', 'D', '>', '<', 'v', '*', 'x', '+']  # Different symbols

	for idx, model_entry in enumerate(model_list[0:], start=0):
		df = pd.read_csv(model_entry[1], sep='\t', skiprows=list(range(0, 5)))

		pcutoff, significantsites = bh_control(list(zip(list(df.site), list(df.P))), 0.05)
		print(pcutoff, len(significantsites))

		yvalues = get_pvalues(df)
		points = get_pvalues(df.head(10))
		hlines = [-1 * transform_p(10, pcutoff, 4), transform_p(10, pcutoff, 4)]
		hlines.sort()

		data = pd.DataFrame({'Model': [model_list[idx][0]] * len(yvalues), 'Values': yvalues})
		sns.set_palette("colorblind")

		sns.violinplot(ax=ax[idx], x='Model', y='Values', data=data)
		ax[idx].set_ylabel('$\omega_r < 1 \leftarrow \log_{10} P \\rightarrow \omega_r > 1$')
		ax[idx].hlines(hlines, *ax[idx].get_xlim(), colors='r', linestyles='dashed')

		for point, symbol, label in zip(points, symbols, df.site):
			ax[idx].scatter(x=0, y=point, marker=symbol, label=str(label))
		ax[idx].legend()
	# break
	# Adjust legend outside of the for loop
	# plt.legend(loc='upper right')
	plt.show()

	plt.savefig(models_summary_plot)


if __name__ == "__main__":
	main()
