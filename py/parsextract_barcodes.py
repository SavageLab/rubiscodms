# == Native Modules ==
# == Installed Modules ==
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
	# == PARAMS
	# == OUTPUTS

	# These files can be big, save them on a drive if necessary
	for i, file in enumerate(listOfFiles):
		print(listOfFileNames[i])
		dfbarcode_list = pd.DataFrame(parseAndExtractBC(file), columns=[listOfFileNames[i]])
		dfbarcode_list.to_csv('/Users/noamprywes/NP_11_66_42/barcode_lists/' + listOfFileNames[i] + '_barcode_list.csv')


if __name__ == "__main__":
	main()
