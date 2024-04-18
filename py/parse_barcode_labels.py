# == Native Modules

# == Installed Modules
import pandas as pd


# == Project Modules


def main():
	# === SNAKEMAKE I/O ===
	# == INPUTS
	biochem_reference_path = str(snakemake.input.biochem_reference_path)
	concat_counts_path = str(snakemake.input.concat_counts_path)
	# == OUTPUTS
	labeled_barcode_path = str(snakemake.output.labeled_barcode_path)
	full_labeled_barcode_path = str(snakemake.output.full_labeled_barcode_path)
	# Load reference biochemistry data

	biochem_reference_data = pd.read_csv(biochem_reference_path)
	df_concat_barcode = pd.read_csv(concat_counts_path, low_memory=False)
	df_concat_barcode = df_concat_barcode.merge(biochem_reference_data, left_on='RbcLCodonMut', right_on='Mutation',
												how='outer')
	df_concat_barcode = df_concat_barcode.convert_dtypes()

	# Label barcodes
	df_concat_barcode.loc[(df_concat_barcode['BackboneMut'] == False) &
						  (df_concat_barcode['InsertionsFound'] == False) &
						  (df_concat_barcode['DeletionsFound'] == False) &
						  (df_concat_barcode['PRKmut'] == False), 'BC_Category'] = 'Verified'

	df_concat_barcode.loc[(df_concat_barcode['BackboneMut'] == True), 'BC_Category'] = 'BackboneMut'
	df_concat_barcode.loc[(df_concat_barcode['InsertionsFound'] == True), 'BC_Category'] = 'InsertionsFound'
	df_concat_barcode.loc[(df_concat_barcode['DeletionsFound'] == True), 'BC_Category'] = 'DeletionsFound'
	df_concat_barcode.loc[(df_concat_barcode['PRKmut'] == True), 'BC_Category'] = 'PRKmut'
	df_concat_barcode.loc[
		(df_concat_barcode['RbcLCodonMut'] == 'Multiple_mutations'), 'BC_Category'] = 'Multiple_mutations'
	df_concat_barcode.loc[(df_concat_barcode['RbcLCodonMut'] == 'Silent mutation'), 'BC_Category'] = 'Silent mutation'
	df_concat_barcode.loc[
		(df_concat_barcode['RbcLCodonMut'] == 'Ambiguous_mutation'), 'BC_Category'] = 'Ambiguous_mutation'
	df_concat_barcode.loc[(df_concat_barcode['RbcLCodonMut'] == 'Illegal_mutation'), 'BC_Category'] = 'Illegal_mutation'

	# Find primer dimer barcodes

	df_concat_barcode.loc[df_concat_barcode['Barcode'].str.contains('CGCGGGGATT'), 'BC_Category'] = 'Primer_Dimer'
	# readsAndRefs.loc[readsAndRefs['Barcode'].str.contains('CGGGCGCGGGG'), 'BC_Category'] = 'Primer_Dimer'
	df_concat_barcode.loc[df_concat_barcode['Barcode'].str.contains('AGGACGCGCGGGGAT'), 'BC_Category'] = 'Primer_Dimer'

	# Label no BC reads

	df_concat_barcode.loc[df_concat_barcode['Barcode'].str.contains('No flank found'), 'BC_Category'] = 'No_BC_in_read'

	# Change remaining barcodes to "Unknown"

	df_concat_barcode[['BC_Category']] = df_concat_barcode[['BC_Category']].fillna(value='Unknown')

	# Only verified barcodes will be used for further analysis
	df_verified_barcodes = df_concat_barcode[df_concat_barcode['BC_Category'] == 'Verified']

	df_verified_barcodes_trim = df_verified_barcodes.drop(columns=['Mutation_List',
														  'BackboneMut',
														  'InsertionsFound',
														  'DeletionsFound',
														  'PRKmut',
														  'Mutation',
														  'BC_Category'], axis=1)


	df_verified_barcodes_trim.to_csv(labeled_barcode_path)
	df_verified_barcodes.to_csv(full_labeled_barcode_path)


if __name__ == "__main__":
	main()
