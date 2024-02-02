# == Native Modules ==
import pickle
# == Installed Modules ==
from Bio import SeqIO
import pysam
import pandas as pd
# == Project Modules ==


def main():
	# Snakemake I/O
	# === Inputs
	# == Define the path to the GenBank file
	genbank_feature_path = str(snakemake.input.genbank_feature_path)
	# == Define the path to the SAM file containing alignment data
	alignment_file_path = str(snakemake.input.alignment_file_path)
	# Outputs
	output_barcode_path = str(snakemake.output.barcode_path)
	output_barcode_count_path = str(snakemake.output.barcode_count_path)
	feature_location_path = str(snakemake.output.feature_location_path)

	# == DEBUG ==
	# Define the path to the GenBank file
	# genbank_feature_path = 'np-11-64-1_ScaI_index.gb'
	# Define the path to the SAM file containing alignment data
	# alignment_file_path = 'NP_11_64_10_PacBio_ScaI_index.bam'

	# Initialize a dictionary to store read names and their associated barcodes
	barcode_dictionary = {'Read_Name': [], 'Barcode': []}

	# Read the GenBank file using SeqIO
	genbank_record = SeqIO.read(genbank_feature_path, "genbank")

	# Initialize a dictionary to store feature types and their corresponding read names
	features_dictionary = {'Read_name': []}

	# Initialize a dictionary to store the locations of each feature
	feature_locations = {}

	# Iterate over the features in the GenBank record
	for feature in genbank_record.features:
		# Add the feature type to the features dictionary if not already present
		if feature.type not in features_dictionary:
			features_dictionary[feature.type] = []

		# Store the feature locations in the dictionary
		feature_locations[feature.type] = [feature.location.start, feature.location.end]

	# Open the SAM file using pysam.AlignmentFile for reading alignment data
	samfile = pysam.AlignmentFile(alignment_file_path)

	# Fetch all alignments from the SAM file, setting until_eof to True to include unmapped reads
	alignment_iterator = samfile.fetch(until_eof=True)

	# Process a subset of alignments for demonstration purposes
	# Note: The full processing was previously completed and took approximately 30 minutes
	for i, read in enumerate(alignment_iterator):
		# if i < 4:  # Limit to the first 4 reads for demonstration
		# Check if the read is mapped in a proper pair and is not a secondary alignment
		# Flags 0 and 16 correspond to primary alignments on the forward and reverse strands, respectively
		if read.flag in (0, 16):
			# Initialize an empty string to accumulate the barcode sequence
			barcode = ''

			# Iterate over aligned pairs to extract the barcode sequence
			for aligned_pair in read.get_aligned_pairs(with_seq=True):
				reference_position = aligned_pair[1]

				try:
					# Check if the reference position falls within the barcode feature location
					# Adjust for zero-based indexing used in Python
					if feature_locations['Barcode'][0] <= reference_position < feature_locations['Barcode'][1]:
						try:
							# Append the nucleotide from the read sequence to the barcode string
							barcode += read.seq[aligned_pair[0]]
						except IndexError:
							# Handle cases where the aligned pair does not have a corresponding read sequence index
							pass
				except TypeError:
					continue

			# Append the read name and extracted barcode to the dictionary
			barcode_dictionary['Read_Name'].append(read.query_name)
			barcode_dictionary['Barcode'].append(barcode)

	# Convert the barcode dictionary to a pandas DataFrame
	barcode_dataframe = pd.DataFrame(barcode_dictionary)
	# Save the DataFrame to a CSV file for further analysis	
	barcode_dataframe.to_csv(output_barcode_path, index=False)

	# Create a new DataFrame that counts the occurrences of each unique barcode
	barcode_counts_df = barcode_dataframe['Barcode'].value_counts().to_frame(name='Counts')

	# Reset the index of the DataFrame to convert the barcodes from the index to a column
	barcode_counts_df.reset_index(inplace=True)
	barcode_counts_df.rename(columns={'index': 'Barcode'}, inplace=True)

	# Filter out barcodes that are shorter than the minimum required length
	# In this case, we are only interested in barcodes that are longer than 10 bases
	minimum_barcode_length = 10
	barcode_length_mask = barcode_counts_df['Barcode'].str.len() > minimum_barcode_length
	filtered_barcode_counts_df = barcode_counts_df.loc[barcode_length_mask]

	# Display the filtered DataFrame
	filtered_barcode_counts_df.to_csv(output_barcode_count_path, index=False)
	with open(feature_location_path, 'ab') as f:
		pickle.dump(feature_locations, f)


if __name__ == "__main__":
	main()
