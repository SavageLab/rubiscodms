# == Native Modules ==
import os
import time
# == Installed Modules ==
import pysam
import pandas as pd
# == Project Modules ==


def main():
    # Snakemake I/O
    # === Inputs
    barcode_path = str(snakemake.input.barcode_path)
    barcode_count_path = str(snakemake.input.barcode_count_path)
    sorted_bam_path = str(snakemake.input.sorted_bam_path)
    # === Params
    sorted_barcode_reads_path = str(snakemake.input.sorted_barcode_reads_path)
    # === Outputs
    bam_barcode_reads_path = str(snakemake.output.bam_barcode_reads_path)
    consensus_fasta_path = str(snakemake.output.consensus_fasta_path)

    # DEBUG
    # barcode_count_path = "/groups/doudna/projects/daniel_projects/prywes_n/rubisco_reads_processing/barcodes/np_11_64_10_ScaI_barcodeCounts.csv"
    # barcode_path = "/groups/doudna/projects/daniel_projects/prywes_n/rubisco_reads_processing/barcodes/np_11_64_10_ScaI_firstPassAllBarcodes1.csv"
    # sorted_bam_path = "/groups/doudna/projects/daniel_projects/prywes_n/rubisco_reads_processing/samtools/toy_sorted.bam"
    # consensus_fasta_path = "/groups/doudna/projects/daniel_projects/prywes_n/rubisco_reads_processing/samtools/test.fasta"
    # bam_barcode_reads_path = "/groups/doudna/projects/daniel_projects/prywes_n/rubisco_reads_processing/samtools/test_barcodes.fasta"
    # bam_mapping_path = 'NP_11_64_10_PacBio_ScaI_index.bam'

    # Sorting the BAM file is a prerequisite for indexing, which facilitates efficient data retrieval
    # pysam.sort("-o", sorted_bam_path, bam_mapping_path)

    # Index the BAM file to enable fast random access
    # This step creates an index file (.bai) associated with the BAM file
    # The index file is crucial for many downstream analyses that require querying specific genomic regions
    # pysam.index(sorted_bam_path)

    # Import barcode tables
    barcode_dataframe = pd.read_csv(barcode_path)
    barcode_count_dataframe = pd.read_csv(barcode_count_path)

    # Open the indexed BAM file using pysam.AlignmentFile
    # This allows us to work with the alignments in the BAM file programmatically
    alignment_bam_file = pysam.AlignmentFile(sorted_bam_path, "rb")

    # Create an IndexedReads object from the BAM file
    # This object allows for random access to reads based on their alignment position
    indexed_bam_file = pysam.IndexedReads(alignment_bam_file)
    indexed_bam_file.build()

    # List of barcodes to process
    barcodesToConsensus = barcode_count_dataframe[barcode_count_dataframe['Counts'] > 1]['Barcode']

    # Start clock
    startTime = time.time()

    with open(consensus_fasta_path, 'w') as fasta_handle:
        # Process each barcode
        for idx, barcode in enumerate(barcodesToConsensus):
            if (idx % 10000 == 0) or (idx == 200) or (idx == 20) or (idx == 1000):
                print(str(100 * idx / len(barcodesToConsensus)) + '% finished in ' + str(
                    time.time() - startTime) + ' seconds')
            # Retrieve reads associated with the current barcode
            barcode_reads = barcode_dataframe[barcode_dataframe['Barcode'] == barcode]
            # Create a temporary BAM file to store reads for the current barcode
            temp_barcode_reads_align = pysam.AlignmentFile(bam_barcode_reads_path, 'wb', template=alignment_bam_file)
            for read_name in barcode_reads['Read_Name']:
                # Find the read in the indexed BAM file and write to the temporary BAM
                try:
                    for alignment in indexed_bam_file.find(read_name):
                        temp_barcode_reads_align.write(alignment)
                except KeyError:
                    continue
            temp_barcode_reads_align.close()
            pysam.sort("-o", sorted_barcode_reads_path, bam_barcode_reads_path)
            consensus = ''

            try:
                consensus = os.popen(
                    f"samtools consensus -l 100000 --config hifi {sorted_barcode_reads_path}").read()
            except Exception as e:
                print(e)
            try:
                if consensus.split()[1]:
                    fasta_handle.write('>' + barcode + '\n')
                    fasta_handle.write(consensus.split()[1] + '\n')
            except IndexError:
                continue


if __name__ == "__main__":
    main()
