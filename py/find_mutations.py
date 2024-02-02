# == Native Modules ==
import pickle
import math
# == Installed Modules ==
import pysam
import pandas as pd
# == Project Modules ==


# Function to find mutation locations
def find_mut_locations(mut_list):
    return [mutation[1] for mutation in mut_list]


# Function to check if backbone mutation exists
def find_backbone_mut(row, feature_locs):
    mut_locs = find_mut_locations(row['Mutation_List'])
    bb_mut_found = False
    for position in mut_locs:
        try:
            if position <= feature_locs['rbcL'][0] or position >= feature_locs['PRK'][1]:
                bb_mut_found = True
            if feature_locs['BCleft'][0] <= position <= feature_locs['BCleft'][1] or \
                    feature_locs['BCright'][0] <= position <= feature_locs['BCright'][1]:
                bb_mut_found = True
        except Exception as e:
            pass
    return bb_mut_found


# Function to check if insertions are found
def find_insertions(row):
    return any(mutation[2] == None for mutation in row['Mutation_List'])


# Function to check if deletions are found
def find_deletions(row):
    return any(mutation[0] == None for mutation in row['Mutation_List'])


# Function to check if PRK mutation is found
def find_prk_mut(row, feature_locs):
    mut_locs = find_mut_locations(row['Mutation_List'])
    prk_mut_found = False
    for position in mut_locs:
        try:
            if position >= feature_locs['PRK'][0] and position <= feature_locs['PRK'][1]:
                prk_mut_found = True
        except Exception as e:
            pass
    return prk_mut_found


# Function to translate DNA sequence to amino acids
def translate(seq):
    table = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N',
              'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L',
              'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q',
              'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
              'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G',
              'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F',
              'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_',
              'TGG':'W' }
    protein = ''
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon] if codon in table else ''
    return protein


# Function to find rbcL mutations
def rbcL_mut(row, feature_locs):
    programmed_mutation_codons = ['CGT', 'CAT', 'AAA', 'GAT', 'GAA',
                                  'AGC', 'ACC', 'AAT', 'CAG', 'TGC',
                                  'GGC', 'CCG', 'GCG', 'GTG', 'ATT',
                                  'CTG', 'ATG', 'TTT', 'TAT', 'TGG']
    rbcLDNA = 'ATGGACCAGAGTTCTCGCTATGTCAATCTTGCTTTGAAAGAAGAGGACTTAATTGCCGGAGGCGAACACGTATTGTGCGCCTATATTATGAAACCAAAAGCTGGGTACGGTTACGTTGCCACCGCAGCGCACTTTGCTGCAGAGAGTTCGACTGGAACAAATGTAGAGGTGTGTACCACTGATGATTTTACTCGTGGCGTGGATGCCTTAGTTTACGAAGTAGATGAGGCCCGCGAGCTTACTAAGATTGCCTACCCAGTAGCACTGTTCGATCGTAATATCACGGATGGAAAAGCTATGATCGCATCATTTCTTACATTGACCATGGGCAATAATCAAGGCATGGGTGACGTAGAATATGCCAAAATGCATGACTTCTATGTGCCTGAAGCGTATCGCGCTCTGTTTGATGGCCCTTCTGTAAATATTTCCGCCTTGTGGAAGGTACTGGGTCGTCCGGAGGTTGATGGAGGCTTAGTTGTCGGTACTATCATTAAACCAAAGCTTGGGCTGCGCCCTAAGCCTTTCGCTGAGGCTTGCCACGCATTTTGGTTAGGTGGTGATTTTATTAAAAATGATGAACCGCAAGGGAACCAACCCTTCGCGCCACTTCGCGACACGATCGCCTTGGTTGCAGATGCGATGCGCCGTGCACAGGACGAAACGGGTGAAGCAAAGTTGTTTAGTGCTAATATCACCGCTGATGACCCGTTTGAAATCATCGCCCGCGGGGAATACGTTCTTGAAACATTCGGTGAAAATGCAAGTCACGTTGCCTTATTAGTAGACGGATACGTAGCCGGCGCAGCTGCTATTACGACCGCTCGTCGCCGTTTCCCCGATAATTTCTTGCATTATCATCGCGCCGGGCATGGAGCGGTCACTAGCCCTCAGTCTAAGCGCGGCTATACTGCTTTTGTGCATTGTAAGATGGCTCGTTTGCAAGGTGCGTCGGGTATCCATACCGGAACTATGGGCTTTGGAAAGATGGAAGGGGAGTCCTCAGATCGTGCCATCGCGTATATGTTAACTCAAGATGAAGCACAAGGCCCCTTTTATCGCCAGTCGTGGGGAGGGATGAAGGCCTGCACCCCTATTATTTCTGGCGGAATGAATGCCTTACGCATGCCCGGCTTCTTCGAAAATCTTGGTAATGCGAATGTCATCCTGACTGCTGGTGGCGGTGCTTTCGGCCATATTGACGGCCCTGTGGCCGGCGCTCGTTCTCTTCGCCAAGCTTGGCAGGCGTGGCGCGATGGCGTCCCTGTGCTTGACTACGCCCGCGAGCATAAGGAGTTGGCACGTGCTTTTGAGTCCTTCCCTGGCGACGCTGACCAGATTTATCCAGGTTGGCGTAAGGCGCTTGGAGTGGAGGATACTCGTTCAGCTTTACCTGCG'
    mut_locs = find_mut_locations(row['Mutation_List'])
    mut_locs_in_frame = []
    affected_codons = []
    rbcL_mut_list = []
    original_codon = 'XXX'
    mutant_codon = 'YYY'
    original_aa = 'X'
    mut_position = -1
    mutant_aa = 'Y'
    aa_mut = 'X0Y'
    for i, position in enumerate(mut_locs):
        try:
            if feature_locs['rbcL'][0] <= position <= feature_locs['rbcL'][1]:
                rbcL_mut_list.append(row['Mutation_List'][i])
                mut_position = math.floor((position - feature_locs['rbcL'][0]) / 3)
                if mut_position not in affected_codons:
                    affected_codons.append(mut_position)
        except Exception as e:
            pass
    if len(affected_codons) == 0:
        aa_mut = 'WT'
        return aa_mut, original_aa, (mut_position + 1), mutant_aa
    try:
        if len(affected_codons) == 1:
            original_codon = rbcLDNA[affected_codons[0] * 3: affected_codons[0] * 3 + 3]
            mut_codon_seq_list = list(original_codon)
            for mutationNumber, mutation in enumerate(rbcL_mut_list):
                within_codon_position = (mutation[1] - feature_locs['rbcL'][0]) % 3
                mut_codon_seq_list[within_codon_position] = mutation[2]
                mutant_codon = ''.join(mut_codon_seq_list)
        else:
            aa_mut = 'Multiple_mutations'
            return aa_mut, original_aa, (mut_position + 1), mutant_aa
    except Exception as e:
        aa_mut = 'Ambiguous_mutation'
        return aa_mut, original_aa, (mut_position + 1), mutant_aa
    if mutant_codon not in programmed_mutation_codons:
        aa_mut = 'Illegal_mutation'
        return aa_mut, original_aa, (mut_position + 1), mutant_aa
    if original_codon != 'XXX' and mutant_codon != 'YYY':
        original_aa = translate(original_codon)
        mutant_aa = translate(mutant_codon)
        if original_aa == mutant_aa:
            aa_mut = mutant_aa = 'Silent mutation'
        elif mutant_aa == '_':
            aa_mut = mutant_aa = 'Nonsense mutation'
        else:
            aa_mut = (original_aa + str(mut_position + 1) + mutant_aa)
    return aa_mut, original_aa, (mut_position + 1), mutant_aa


def main():
    # Snakemake I/O
    # === Inputs
    aligned_consensus_path = str(snakemake.input.aligned_consensus_path)
    feature_location_path = str(snakemake.input.feature_location_path)
    # === Outputs
    mutation_table = str(snakemake.output.mutation_table)
    filtered_mutation_table = str(snakemake.output.filtered_mutation_table)

    # Look through consensus sequences and find mutations
    # This code might handle mutations at the edge of the read oddly, though those are less likely to be real anyway
    bam_file = pysam.AlignmentFile(aligned_consensus_path)
    with open(feature_location_path, 'rb') as feature_location_handle:
        feature_location = pickle.load(feature_location_handle)
    mutation_dict = {}

    for a, BCConsensus in enumerate(bam_file):
        aligned_pairs = BCConsensus.get_aligned_pairs(with_seq=True)
        mutations = []
        refPosNotNone = 0
        for i, (read_pos, ref_pos, ref_base) in enumerate(aligned_pairs):
            # Prevent reference position from being interpreted as None
            if ref_pos != None:
                refPosNotNone = ref_pos
            if ref_base is not None:
                ref_base = ref_base.upper()
            if read_pos is None:
                read_base = None
            else:
                read_base = BCConsensus.query_sequence[read_pos].upper()
            if (read_base != ref_base) and (ref_base != 'N'):
                mutations = mutations + [[ref_base, refPosNotNone, read_base]]

        # This makes a dataframe where each barcode is paired with a list of lists of mutations
        mutation_dict[BCConsensus.qname] = [mutations]
    mutation_df = pd.DataFrame(mutation_dict).transpose().reset_index().rename(columns={'index': 'Barcode',
                                                                                        0: 'Mutation_List'})

    # Dataframe operations
    mutation_df['BackboneMut'] = mutation_df.apply(lambda row: find_backbone_mut(row, feature_location), axis=1)
    mutation_df['InsertionsFound'] = mutation_df.apply(lambda row: find_insertions(row), axis=1)
    mutation_df['DeletionsFound'] = mutation_df.apply(lambda row: find_deletions(row), axis=1)
    mutation_df['PRKmut'] = mutation_df.apply(lambda row: find_prk_mut(row, feature_location), axis=1)
    mutation_df['RbcLCodonMut'], mutation_df['originalAA'], mutation_df['AApos'], mutation_df['mutAA'] = \
        zip(*mutation_df.apply(rbcL_mut, axis=1))

    # Writing dataframe to CSV
    mutation_df.to_csv(mutation_table)

    # Filtering away barcodes that break the rules
    better_muts_df = mutation_df[~(mutation_df['BackboneMut'] |
                                   mutation_df['DeletionsFound'] |
                                   mutation_df['PRKmut'])]

    better_muts_df = better_muts_df[~(better_muts_df['RbcLCodonMut'].isin(['Multiple_mutations', 'Silent mutation',
                                                                           'Ambiguous_mutation', 'Illegal_mutation']))]

    filtered_mutation_df = better_muts_df[better_muts_df.groupby('RbcLCodonMut')['RbcLCodonMut'].transform('size') >= 3]
    filtered_mutation_df[['Barcode', 'RbcLCodonMut', 'originalAA', 'AApos', 'mutAA']].to_csv('NP_11_64_10_lookupTable061623.csv')

    filtered_mutation_df.to_csv(filtered_mutation_table)


if __name__ == "__main__":
    main()
