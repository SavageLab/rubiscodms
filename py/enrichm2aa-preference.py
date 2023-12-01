# Installed modules
import pandas as pd
from matplotlib import pyplot as plt
from Bio import SeqIO


def main():
	# Snakemake I/O
	# Inputs
	dms_in = str(snakemake.input.dms_data)
	# Outputs
	dms_aapref_out = str(snakemake.output.dms_prefs)
	# Params
	position_col = str(snakemake.params.position_col)
	enrichment_col = str(snakemake.params.enrichment_col)
	aminoacid_col = str(snakemake.params.aminoacid_col)

	# DEBUG INPUT
	# dms_in = "/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/scratch/betalactamase/betaLactamase_enrichmentScores.txt" # "/groups/doudna/projects/daniel_projects/prywes_n/input_data/dms_data.csv"
	# position_col = "Site"
	# enrichment_col = "Trial_1_AmpConc_2500" # "5percent_CO2_20uM_IPTG"
	# aminoacid_col = "AminoAcid"

	minenrichment = 1.0e-4  # minimum allowed enrichment
	df = pd.read_csv(dms_in)

	df["preference"] = [max(minenrichment, 10**df[enrichment_col][x]) for x in range(len(df))]

	# Format the AA preference table
	df_aapref = df.pivot(index=position_col, columns=aminoacid_col, values="preference")
	df_aapref.fillna(1, inplace=True)
	df_aapref = df_aapref.div(df_aapref.sum(axis=1), axis=0)
	df_aapref.insert(0, "site", range(1, len(df_aapref) + 1))

	# Export output
	df_aapref.to_csv(dms_aapref_out, index=False)

if __name__ == "__main__":
	main()
