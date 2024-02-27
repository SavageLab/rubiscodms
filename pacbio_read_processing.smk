# **** Variables ****
configfile: "config/pacbio_read_processing.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile pacbio_read_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --rerun-incomplete --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#   Minimap2 alignment
		expand("{run}/minimap2/{experiment_id}.bam",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Extract barcodes
		expand("{run}/barcodes/{experiment_id}_firstPassAllBarcodes1.csv",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Plot barcodes
		expand("{run}/figures/{experiment_id}_barcodePlot.png",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Recover barcode reads
		expand("{run}/samtools/{experiment_id}_sorted.bam",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Align consensus sequences
		expand("{run}/samtools/{experiment_id}_aligned_consensus.bam",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Annotate mutations
		expand("{run}/mutations/{experiment_id}_full_mutation_table.csv",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Plot mutation-associated statistics
		expand("{run}/mutations/{experiment_id}_mutation_statistics_table.csv",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Conduct data rarefaction
		expand("{run}/figures/{experiment_id}_rarefaction_plot.png",
			run=config["run"],experiment_id=config["experiment_id"]),
		#   Plot mutation positions
		expand("{run}/figures/{experiment_id}_mutation_positions_negative.png",
			run=config["run"],experiment_id=config["experiment_id"]),

rule minimap2_map:
	input:
		fasta_feature_path = lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.fasta".format(
			in_dir=config['input_dir'], experiment_id=wildcards.experiment_id)),
		p1 = lambda wildcards: glob.glob("{in_dir}/{experiment_id}_pacbio.fastq".format(
			in_dir=config['input_dir'], experiment_id=wildcards.experiment_id))
	output:
		bam_map_path = "{run}/minimap2/{experiment_id}.bam",
		sorted_bam_path= "{run}/samtools/{experiment_id}_sorted.bam",
	conda:
		"envs/mapping.yaml"
	params:
		sam_map = "{run}/minimap2/{experiment_id}.sam"
	message:
		"""
Mapping read pairs:
R1: {input.p1}
Reference:
{input.fasta_feature_path}
Mapped Output:
{params.sam_map}
        """
	threads:
		config["threads"]
	shell:
		"""
		minimap2 --MD -Lax map-pb {input.fasta_feature_path} {input.p1} > {params.sam_map} 
		samtools view -@ {threads} -bT {input.fasta_feature_path} {params.sam_map} > {output.bam_map_path}
		samtools sort -@ {threads} {output.bam_map_path} -o {output.sorted_bam_path}
		samtools index -@ {threads} -b {output.sorted_bam_path}
		rm {params.sam_map}
		"""

rule extract_barcodes:
	input:
		genbank_feature_path = lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.gb".format(
			in_dir=config['input_dir'],
			experiment_id=wildcards.experiment_id)),
		bam_map_path = "{run}/minimap2/{experiment_id}.bam"
	output:
		barcode_path = "{run}/barcodes/{experiment_id}_firstPassAllBarcodes1.csv",
		barcode_count_path = "{run}/barcodes/{experiment_id}_barcodeCounts.csv",
		feature_location_path = "{run}/barcodes/{experiment_id}_feature_location.pkl",
	conda:
		"envs/barcodes_processing.yaml"
	message:
		"""
rule: extract_barcodes
Import feature in GenBank format:\n {input.genbank_feature_path}
Export barcode tables to:\n {output.barcode_path}\n {output.barcode_count_path}
		"""
	script:
		"py/extract_barcodes.py"

rule barcode_plot:
	input:
		barcode_path="{run}/barcodes/{experiment_id}_firstPassAllBarcodes1.csv",
		barcode_count_path="{run}/barcodes/{experiment_id}_barcodeCounts.csv"
	output:
		barcode_plot = "{run}/figures/{experiment_id}_barcodePlot.png",
	conda:
		"envs/data_handling.yaml"
	message:
		"""
rule: barcode_plot 		
Import  barcode tables:\n {input.barcode_path}\n {input.barcode_count_path}
Export barcode usability plot to:\n {output.barcode_plot}
		"""
	script:
		"py/barcode_usability_plot.py"

rule pick_barcode_reads:
	input:
		barcode_path= "{run}/barcodes/{experiment_id}_firstPassAllBarcodes1.csv",
		sorted_bam_path= "{run}/samtools/{experiment_id}_sorted.bam"
	output:
		bam_barcode_reads_path = "{run}/samtools/{experiment_id}_barcode_reads.bam",
		consensus_fasta_path = "{run}/samtools/{experiment_id}_barcode_reads.fasta"
	conda:
		"envs/barcodes_processing.yaml"
	message:
		"""
rule: pick_barcode_reads
Import  BAM alignment:\n {input.sorted_bam_path}\n
Export consensus FASTA sequences to:\n {output.consensus_fasta_path}
		"""
	script:
		"py/barcode_consensus.py"

rule align_consensus:
	input:
		fasta_feature_path=lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.fasta".format(
			in_dir=config['input_dir'],experiment_id=wildcards.experiment_id)),
			consensus_fasta_path = "{run}/samtools/{experiment_id}_barcode_reads.fasta"
	output:
		aligned_consensus_path="{run}/samtools/{experiment_id}_aligned_consensus.bam"
	params:
		sam_consensus = "{run}/samtools/{experiment_id}_aligned_consensus.sam"
	conda:
		"envs/mapsam.yaml"
	message:
		"""
rule: align_consensus 		
Import consensus barcode FASTA:\n {input.consensus_fasta_path}
Export aligned consensus BAM to:\n {output.aligned_consensus_path}
		"""
	shell:
		"""
		minimap2 --MD -Lax map-hifi {input.fasta_feature_path} {input.consensus_fasta_path} > {params.sam_consensus}
		samtools view -bT {input.fasta_feature_path} {params.sam_consensus} > {output.aligned_consensus_path}
		rm {params.sam_consensus}
		"""

rule find_mutations:
	input:
		aligned_consensus_path="{run}/samtools/{experiment_id}_aligned_consensus.bam",
		feature_location_path = "{run}/barcodes/{experiment_id}_feature_location.pkl"
	output:
		full_mutation_table = "{run}/mutations/{experiment_id}_full_mutation_table.csv",
		filtered_mutation_table = "{run}/mutations/{experiment_id}_filtered_mutation_table.csv",
		by_mutation_filtered_table = "{run}/mutations/{experiment_id}_by_mutation_filtered_table.csv"
	conda:
		"envs/barcodes_processing.yaml"
	message:
		"""
rule: find_mutations		
Import aligned consensus BAM :\n {input.feature_location_path}\n {input.aligned_consensus_path}
Export mutation table to:\n {output.full_mutation_table}
		"""
	script:
		"py/find_mutations.py"

rule mutation_statistics:
	input:
		filtered_mutation_table = "{run}/mutations/{experiment_id}_filtered_mutation_table.csv",
		barcode_count_path = "{run}/barcodes/{experiment_id}_barcodeCounts.csv"
	output:
		mutation_stats_report = "{run}/mutations/{experiment_id}_mutation_statistics_table.csv",
		barcode_per_mut_plot_output = "{run}/figures/{experiment_id}_barcode_per_mut_plot.png",
		mutation_heatmap_output = "{run}/figures/{experiment_id}_mutation_heatmap.png"
	conda:
		"envs/data_handling.yaml"
	message:
		"""
rule: mutation_statistics
Import barcode and mutation data from:\n {input.filtered_mutation_table}\n {input.barcode_count_path}
Calculate statistics and report on:\n {output.mutation_stats_report}\n
		"""
	script:
		"py/mutation_statistics.py"

rule data_rarefaction:
	input:
		barcode_count_path = "{run}/barcodes/{experiment_id}_barcodeCounts.csv",
	output:
		rarefaction_plot_output = "{run}/figures/{experiment_id}_rarefaction_plot.png"
	conda:
		"envs/scipy.yaml"
		# from matplotlib import pyplot as plt
		# import numpy as np
		# import pandas as pd
		# from scipy.optimize import curve_fit
	message:
		"""
rule: data_rarefaction 		
Import barcode counts:\n {input.barcode_count_path}
Plot rarefaction to:\n {output.rarefaction_plot_output}
		"""
	script:
		"py/data_rarefaction.py"

rule plot_mutation_positions:
	input:
		full_mutation_table = "{run}/mutations/{experiment_id}_full_mutation_table.csv",
	output:
		mutation_positions_plot = "{run}/figures/{experiment_id}_mutation_positions.png",
		mutation_positions_negative_plot = "{run}/figures/{experiment_id}_mutation_positions_negative.png"
	conda:
		"envs/data_handling.yaml"
	message:
		"""
rule: plot_mutation_positions		
Import mutations table: {input.full_mutation_table}
Generate mutation position plot at: {output.mutation_positions_plot}
		"""
	script:
		"py/plot_mutation_positions.py"
