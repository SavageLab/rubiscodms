# == Native Modules ==
import os
# == Installed Modules ==
from Bio import SeqIO
import pysam
import pandas as pd
# == Project Modules ==


def main():
    # Snakemake I/O
    # === Inputs
    barcode_path = str(snakemake.input.barcode_path)
    sorted_bam_path = str(snakemake.input.sorted_bam_path)
    # === Outputs
    bam_barcode_reads_path = str(snakemake.output.bam_barcode_reads_path)
    consensus_fasta_path = str(snakemake.output.consensus_fasta_path)

    # Import barcode tables
    barcode_dataframe = pd.read_csv(barcode_path)

    # DEBUG
    # bam_mapping_path = 'NP_11_64_10_PacBio_ScaI_index.bam'

    # Sorting the BAM file is a prerequisite for indexing, which facilitates efficient data retrieval
    # pysam.sort("-o", sorted_bam_path, bam_mapping_path)

    # Index the BAM file to enable fast random access
    # This step creates an index file (.bai) associated with the BAM file
    # The index file is crucial for many downstream analyses that require querying specific genomic regions
    # pysam.index(sorted_bam_path)

    # Open the indexed BAM file using pysam.AlignmentFile
    # This allows us to work with the alignments in the BAM file programmatically
    alignment_bam_file = pysam.AlignmentFile(sorted_bam_path, "rb")

    # Create an IndexedReads object from the BAM file
    # This object allows for random access to reads based on their alignment position
    indexed_bam_file = pysam.IndexedReads(alignment_bam_file)
    indexed_bam_file.build()

    # List of barcodes to process
    barcode_list = barcode_dataframe['Barcode'].tolist()

    # Open the template BAM file for reading
    with pysam.AlignmentFile(sorted_bam_path, 'rb') as template_bam:
        # Open the output FASTA file for writing
        with open(consensus_fasta_path, 'w') as output_fasta:
            # Process each barcode
            for barcode in barcode_list:
                # Retrieve reads associated with the current barcode
                barcode_reads = barcode_dataframe[barcode_dataframe['Barcode'] == barcode]
                # Create a temporary BAM file to store reads for the current barcode
                with pysam.AlignmentFile(bam_barcode_reads_path, 'wb', template=template_bam) as temp_bam:
                    for read_name in barcode_reads['Read_Name']:
                        try:
                            # Find the read in the indexed BAM file and write to the temporary BAM
                            for alignment in indexed_bam_file.find(read_name):
                                temp_bam.write(alignment)
                        except Exception as e:
                            raise Exception(f"Error processing read {read_name}: {e}")

                # Write the barcode as the header in the FASTA file
                # seq_record = SeqIO.SeqRecord(id=barcode, seq=)
                # output_fasta.write(f'>{barcode}\n')

                # Generate the consensus sequence using samtools
                try:
                    consensus_cmd = f'samtools consensus -l 100000 --config hifi {sorted_bam_path}'
                    consensus_result = os.popen(consensus_cmd).read()
                    consensus_sequence = consensus_result.split()[1]
                    # Write the barcode as the header in the FASTA file
                    seq_record = SeqIO.SeqRecord(id=barcode, seq=consensus_sequence)
                    SeqIO.write(seq_record, output_fasta, 'fasta')
                except Exception as e:
                    raise Exception(f"Error generating consensus for barcode {barcode}: {e}")


if __name__ == "__main__":
    main()
