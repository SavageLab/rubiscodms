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
		expand("{run}/barcodes/{file_prefix}_pacbioMerged_barcodeCounts.csv",
			run=config["run_nextseq"],file_prefix=config["file_prefix"]),
		expand("{run}/barcodes/master_concat_barcodes.csv",
			run=config["run_nextseq"]),
		expand("{run}/result_tables/labeled_barcodeCounts.csv",
			run=config["run_nextseq"]),
		expand("{run}/figures/{file_prefix}_pie_chart.png",
			run=config["run_nextseq"],file_prefix=config["file_prefix"]),
		expand("{run}/result_tables/enrich_parameter_sweep",
			run=config["run_nextseq"])

rule extract_barcodes:
	input:
		fastq_reads = lambda wildcards: glob.glob("{in_dir}/{file_prefix}.fastq".format(
			in_dir=config["input_dir_nextseq"], file_prefix=wildcards.file_prefix)),
		pacbio_barcode_path = lambda wildcards: glob.glob("{run_pacbio}/mutations/{experiment_id}_full_mutation_table.csv".format(
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
Barcode Report:
{output.pacbio_merged_counts}
        """
	script:
		"py/parsextract_barcodes.py"

rule merge_barcode_reports:
	input:
		barcode_report_list = expand("{{run}}/barcodes/{file_prefix}_pacbioMerged_barcodeCounts.csv",
			file_prefix=config["file_prefix"])
	output:
		concat_counts_path = "{run}/barcodes/master_concat_barcodes.csv",
	conda:
		"envs/pyplot.yaml"
	message:
		"""
		Generate the concatenated report: {output.concat_counts_path}
		By merging the following reports:
		{input.barcode_report_list} 
		"""
	script:
		"py/concat_barcodes.py"

rule add_biochemistry:
	input:
		biochem_reference_path = lambda wildcards: glob.glob("{input_dir}/biochemData061823.csv".format(
			input_dir=config['input_dir_nextseq'])),
		concat_counts_path = "{run}/barcodes/master_concat_barcodes.csv"
	output:
		labeled_barcode_path="{run}/result_tables/labeled_barcodeCounts.csv",
		full_labeled_barcode_path="{run}/result_tables/full_labeled_barcodeCounts.csv",
	conda:
		"envs/pyplot.yaml"
	message:
		"""
Add biochemistry data: {input.biochem_reference_path}
Verified and labeled barcode data: {input.concat_counts_path}
Generate files labelled with bichemical data:
{output.labeled_barcode_path}
{output.full_labeled_barcode_path}
		"""
	script:
		"py/parse_barcode_labels.py"

rule visualize_sample_stats:
	input:
		full_labeled_barcode_path="{run}/result_tables/full_labeled_barcodeCounts.csv"
	output:
		sample_pie_chart = "{run}/figures/{file_prefix}_pie_chart.png"
	conda:
		"envs/pyplot.yaml"
	message:
		"""
Visualize Barcode data from sample: {wildcards.file_prefix}
Import data from: {input.full_labeled_barcode_path}
Export figure to : {output.sample_pie_chart}
		"""
	script:
		"py/visualize_sample_stats.py"

rule calculate_enrichment:
	input:
		labeled_barcode_path = "{run}/result_tables/labeled_barcodeCounts.csv"
	output:
		parameter_sweep_plot = "{run}/figures/parameter_sweep_plot.png",
		parameter_sweep_table = "{run}/result_tables/enrich_parameter_sweep",
	conda:
		"envs/pyplot.yaml"
	message:
		"""
Import data from: {input.labeled_barcode_path}
Export table and figue respectively to : 
{output.parameter_sweep_table}
{output.parameter_sweep_plot}
		"""
	script:
		"py/calculate_enrichment.py"

rule bootstrap_annotate:
	input:
		parameter_sweep_table = "{run}/result_tables/enrich_parameter_sweep",
	output:
		bootstrap_data = "{run}/result_tables/bootstrap_data"
	conda:
		"envs/pyplot.yaml"
	message:
		"""
Import data from: {input.parameter_sweep_table}
Export table and figue to : 
{output.bootstrap_data}
		"""
	script:
		"py/bootstrap_annotate.py"
