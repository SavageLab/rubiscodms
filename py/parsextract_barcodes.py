# == Native Modules ==
# == Installed Modules ==
import pandas as pd
from Bio import SeqIO
# == Project Modules ==


def parseAndExtractBC(file, flanking_sequence):
	"""
	# Parse and extract barcode, make csv files for each one
	# Check for flank even if there's one mistake
	:param file: 
	:param flanking_sequence: 
	:return: 
	"""
	flankSeqsOBO = []
	for position, letter in enumerate(flanking_sequence):
		for nuc in ['A', 'T', 'G', 'C']:
			if flanking_sequence[position] != nuc:
				flankSeqsOBO.append(flanking_sequence[:position] + nuc + flanking_sequence[position + 1:])
	
	fastq_records = SeqIO.parse(file, "fastq")
	barcode_list = []
	for readCount, rec in enumerate(fastq_records):
		if readCount % 3000000 == 0:
			print(readCount)
		readSeq = str(rec.seq)
		barcode_location = readSeq.find(flanking_sequence)
		if barcode_location == -1:
			for OBOseq in flankSeqsOBO:
				barcode_location = readSeq.find(OBOseq)
				if barcode_location > -1:
					break
		if barcode_location == -1:
			barcode_list.append('No flank found')
		else:
			BC = readSeq[barcode_location + len(flanking_sequence): barcode_location + len(flanking_sequence) + 30]
			if 'N' in BC:
				barcode_list.append('BC has Ns')
			else:
				barcode_list.append(BC)
	return barcode_list


def main():
	# === SNAKEMAKE I/O ===
	# == INPUTS
	alignment_file_path = str(snakemake.input.fastq_reads)
	pacbio_barcode_path = str(snakemake.input.pacbio_barcode_path)
	# == PARAMS
	flanking_sequence = str(snakemake.params.flanking_sequence)
	# == OUTPUTS
	barcode_path = str(snakemake.output.barcode_path)
	barcode_count_path = str(snakemake.output.barcode_count_path)
	pacbio_merged_counts_path = str(snakemake.output.pacbio_merged_counts)
	# == WILDCARDS
	file_prefix = str(snakemake.wildcards.file_prefix)

	# DEBUG INPUT
	# alignment_file_path = "/groups/doudna/projects/daniel_projects/rubiscodms/input_data_nextseq/NP_11_66_1.fastq"
	# flanking_sequence = "TTCCGTACGA"
	# file_prefix = "NP_11_66_1"
	# pacbio_barcode_path = "/groups/doudna/projects/daniel_projects/rubiscodms/rubisco_reads_processing/barcodes/np_11_64_10_ScaI_firstPassAllBarcodes1.csv"

	df_pacbio_barcode = pd.read_csv(pacbio_barcode_path).drop(columns = 'Unnamed: 0')
	#
	df_barcode_list = pd.DataFrame(parseAndExtractBC(alignment_file_path, flanking_sequence), columns=[file_prefix])
	#
	df_barcode_counts = pd.DataFrame(df_barcode_list[file_prefix].value_counts())

	df_barcode_counts.index.name = 'Barcode'
	df_merged_counts = df_pacbio_barcode.merge(df_barcode_counts.rename(columns={'count': file_prefix}), on='Barcode', how='outer')
	df_merged_counts = df_merged_counts.convert_dtypes()

	#
	df_barcode_counts.to_csv(barcode_count_path, index=True)
	#
	df_barcode_list.to_csv(barcode_path)
	#
	df_merged_counts.to_csv(pacbio_merged_counts_path)


if __name__ == "__main__":
	main()
