# **** Variables ****
configfile: "config/nextseq_read_processing.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile pacbio_read_processing.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --rerun-incomplete --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		#
		expand("{run}/minimap2/{experiment_id}.bam",
			run=config["run"],experiment_id=config["experiment_id"]),


rule extract_barcodes:
	input:
		fasta_feature_path = lambda wildcards: glob.glob("{in_dir}/{file_prefix}.fastq".format(
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
	script:
		"py/parsextract_barcodes.py"
