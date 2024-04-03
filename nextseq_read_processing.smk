# **** Variables ****
configfile: "config/nextseq_read_processing.yaml"
configfile: "config/pacbio_read_processing.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile pacbio_read_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --rerun-incomplete --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#
		expand("{run}/barcodes/{file_prefix}_parsed_barcodes.csv",
			run=config["run_nextseq"],file_prefix=config["file_prefix"]),

rule extract_barcodes:
	input:
		fastq_reads = lambda wildcards: glob.glob("{in_dir}/{file_prefix}.fastq".format(
			in_dir=config["input_dir_nextseq"], file_prefix=wildcards.file_prefix)),
		pacbio_barcode_path = lambda wildcards: glob.glob("{run_pacbio}/barcodes/{experiment_id}_barcodeCounts.csv".format(
		run_pacbio=config['run_pacbio'],experiment_id=config["experiment_id"]))
	output:
		barcode_path = "{run}/barcodes/{file_prefix}_parsed_barcodes.csv",
		barcode_count_path = "{run}/barcodes/{file_prefix}_barcodeCounts.csv",
		pacbio_merged_counts = "{run}/barcodes/{file_prefix}_pacbioMerged_barcodeCounts.csv",
	params:
		flanking_sequence=config["flankSeq"]
	conda:
		"envs/pyplot.yaml"
	message:
		"""
Extracting Barcodes from reads:
FASTQ: {input.fastq_reads}
Flanking Sequence: {params.flanking_sequence}
Bracode Report:
{output.barcode_path}
        """
	script:
		"py/parsextract_barcodes.py"

rule merge_barcode_reports:
	input:
		barcode_report_list = expand("{run}/barcodes/{file_prefix}_pacbioMerged_barcodeCounts.csv".format(
			run=config["run_nextseq"],file_prefix=config["file_prefix"])),
	output:
		concat_counts_path = "{run}/barcodes/master_concat_barcodes.csv",
	conda:
		"envs/pyplot.yaml"
	script:
		"py/concat_barcodes.py"

rule add_biochemistry:
	input:
		biochem_reference_data = lambda wildcards: glob.glob("{input_dir}/biochemData061823.csv".format(
			input_dir=config['input_dir_nextseq'])),
		concat_counts_path = "{run}/barcodes/master_concat_barcodes.csv"
	output:
		labeled_barcode_count_path="{run}/barcodes/{file_prefix}_labeled_barcodeCounts.csv",
	conda:
		"envs/pyplot.yaml"
	message:
		"""
Extracting Barcodes from reads:
FASTQ: {input.fastq_reads}
Flanking Sequence: {params.flanking_sequence}
Bracode Report:
{output.barcode_path}
		"""
	script:
		"py/parsextract_barcodes.py"