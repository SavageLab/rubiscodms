# **** Variables ****
configfile: "config/pacbio_read_processing.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile dms_workflow.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}/{experiment_id}_firstPassAllBarcodes1.csv",
			run=["run"], experiment_id=["experiment_id"]),
		expand("{run}/{experiment_id}_barcodePlot.png",
		       run=["run"], experiment_id=["experiment_id"]),
		expand("{run}/minimap2/{experiment_id}_index.bam",
		       run=["run"], experiment_id=["experiment_id"]),

rule hisat2_map:
	input:
		index = "{run}/{experiment_id}_index.fasta",
		p1 = "{run}/{experiment_id}_PacBio.fastq",
	output:
		bam_map_path ="{run}/minimap2/{experiment_id}.bam"
	conda:
		"envs/mapping.yaml"
	params:
		sam_map = "{run}/minimap2/{experiment_id}.sam",
	message:
		"""
Mapping read pairs:
R1: {input.p1}
Reference:
{input.index}
Mapped Output:
{params.sam_map}
        """
	threads:
		config["threads"]
	shell:
		"""
		minmap2 
		samtools view -bT {input.index} {params.sam_map} > {output.bam_map_path}
		rm {params.sam_map}
		"""

rule extract_barcodes:
	input:
		genbank_feature_path = lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.gb".format(
			in_dir=config['input_dir'],
			experiment_id=wildcards.experiment_id))
	output:
		barcode_path = "{run}/{experiment_id}_firstPassAllBarcodes1.csv",
		barcode_count_path = "{run}/{experiment_id}_barcodeCounts.csv"
	conda:
		"envs/samtools.yaml"
	message:
		"""
		Import feature in GenBank format:\n {input.genbank_feature_path}
		Export barcode tables to:\n {output.barcode_path}\n {output.barcode_count_path}		
		"""
	script:
		"py/extract_barcodes.py"

rule barcode_plot:
	input:
		barcode_path="{run}/{experiment_id}_firstPassAllBarcodes1.csv",
		barcode_count_path="{run}/{experiment_id}_barcodeCounts.csv"
	output:
		barcode_plot = "{run}/{experiment_id}_barcodePlot.png",
	conda:
		"envs/samtools.yaml"
	message:
		"""
		Import  barcode tables:\n {input.barcode_path}\n {input.barcode_count_path}
		Export barcode usability plot to:\n {output.barcode_plot}		
		"""
	script:
		"py/barcode_usability_plot.py"

rule pick_barcode_reads:
	input:
		bam_map_path = "{run}/minimap2/{experiment_id}.bam",
	output:
		sorted_bam_path = "{run}/minimap2/{experiment_id}_sorted.bam",
		bam_barcode_reads_path = "{run}/minimap2/{experiment_id}_barcode_reads.bam",
		sorted_bam_barcode_reads_path = "{run}/minimap2/{experiment_id}_sorted_barcode_reads.bam",
		consensus_fsata_path = "{run}/minimap2/{experiment_id}_barcode_reads.fasta"
	conda:
		"envs/samtools.yaml"
	message:
		"""
		Import  BAM alignment:\n {input.bam_map_path}\n
		Export consensus FASTA sequences to:\n {output.consensus_fsata_path}		
		"""
