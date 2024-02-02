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

rule minimap2_map:
	input:
		fasta_feature_path = lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.fasta".format(
			in_dir=config['input_dir'], experiment_id=wildcards.experiment_id)),
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
{input.fasta_feature_path}
Mapped Output:
{params.sam_map}
        """
	threads:
		config["threads"]
	shell:
		"""
		minmap2 
		samtools view -bT {input.fasta_feature_path} {params.sam_map} > {output.bam_map_path}
		rm {params.sam_map}
		"""

rule extract_barcodes:
	input:
		genbank_feature_path = lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.gb".format(
			in_dir=config['input_dir'],
			experiment_id=wildcards.experiment_id))
	output:
		barcode_path = "{run}/{experiment_id}_firstPassAllBarcodes1.csv",
		barcode_count_path = "{run}/{experiment_id}_barcodeCounts.csv",
		feature_location_path = "{run}/{experiment_id}_feature_location.pkl",
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
	script:
		"py/barcode_consensus.py"

rule align_consensus:
	input:
		fasta_feature_path=lambda wildcards: glob.glob("{in_dir}/{experiment_id}_index.fasta".format(
			in_dir=config['input_dir'],experiment_id=wildcards.experiment_id)),
			consensus_fasta_path = "{run}/minimap2/{experiment_id}_barcode_reads.fasta"
	output:
		aligned_consensus_path="{run}/minimap2/{experiment_id}_aligned_consensus.bam"
	params:
		sam_consensus = "{run}/minimap2/{experiment_id}_aligned_consensus.sam"
	conda:
		"envs/samtools.yaml"
	message:
		"""
		Import consensus barcode FASTA:\n {input.consensus_fasta_path}\n
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
		aligned_consensus_path="{run}/minimap2/{experiment_id}_aligned_consensus.bam",
		feature_location_path = "{run}/{experiment_id}_feature_location.pkl"
	output:
		mutation_table = "{run}/minimap2/{experiment_id}_mutation_table.csv",
		filtered_mutation_table = "{run}/minimap2/{experiment_id}_filtered_mutation_table.csv"
	conda:
		"envs/samtools.yaml"
	message:
		"""
		Import aligned consensus BAM :\n {input.aligned_consensus_path}\n
		Export mutation table to:\n {output.mutation_table}		
		"""
	script:
		"py/find_mutations.py"

