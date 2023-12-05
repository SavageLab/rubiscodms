# == Native Modules
from Bio import AlignIO, SeqIO
# == Installed Modules
import pandas as pd
# == Project Modules


def main():
	# Snakemake I/O
	# Inputs
	dms_prefs = str(snakemake.input.dms_prefs)
	msa = str(snakemake.input.msa)
	refseq_path = str(snakemake.input.refseq)
	# Outputs
	matched_alignments = str(snakemake.output.matched_alignments)
	# Params
	position_col = str(snakemake.params.position_col)

	# DEBUG INPUT
	# dms_prefs = "/groups/doudna/projects/daniel_projects/prywes_n/pgym_input_data/YAP1_HUMAN_Araya_2012.csv"
	# refseq_path = "/groups/doudna/projects/daniel_projects/prywes_n/pgym_input_data/fasta/YAP1_HUMAN.fasta"
	# position_col = 'Site'
	# enrichment_col = 'DMS_score'
	# aminoacid_col = 'Mt_aminoacid'
	# wt_aminoacid_col = "Wt_aminoacid"

	# === Import Substitution Tables ===
	df = pd.read_csv(dms_prefs)
	(site_min, site_max) = (df[position_col].astype(int).min(), df[position_col].astype(int).max())

	# === Import Reference
	refseq = SeqIO.read(open(refseq_path), "fasta")

	# === Import Alignment
	alignment = AlignIO.read(open(msa), "fasta")

	for aligned_record in alignment:
		if len(aligned_record.seq) != len(refseq.seq):
			print("#### WARNING ####\nThe MSA provided has a different length than that of the original reference sequence"
			      "\nThis rule will assume the MSA has been previously cut to the length of the DMS."
			      "\nIf this is not the case, phydms_comprehensive will issue an error."
			      "\nPlease provide the alignment that includes the original length so this rule can make the necessary"
			      "adjustments.\n")
			break
		aligned_record.seq = aligned_record.seq[(site_min - 1):(site_max - 1)]

	with open(matched_alignments, 'w') as f:
		AlignIO.write(alignment, f, 'fasta')


if __name__ == "__main__":
	main()
