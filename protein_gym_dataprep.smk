# **** Variables ****
configfile: "config/prywes_pgym_dms.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile protein_gym_dataprep.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}_{min_ident}/{experiment_id}/sequence_ids/id_list.txt", run=config["run"], experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/fasta/multi_fasta.fna", run=config["run"], experiment_id=config["reference_seq"],min_ident=config["min_ident"])

# noinspection SmkAvoidTabWhitespace
rule fetch_msa_ids:
	input:
		msa = lambda wildcards: glob.glob("{input_dir}/{experiment_prefix}/processed_inputs/nt_alignment_msa.fna".format(
			input_dir=config["input_dir"], experiment_prefix=wildcards.experiment_id))
	output:
		id_list = "{run}_{min_ident}/{experiment_id}/sequence_ids/id_list.txt"
	conda:
		"envs/dms.yaml"
	message:
		"""
Fetch sequence IDs from {input.msa}
Export ID list to {output.id_list}	
		"""
	shell:
		"""
		python3 py/id_from_msa.py {input.msa} {output.id_list}
		"""
rule protein_to_locus:
	input:
		id_list = "{run}_{min_ident}/{experiment_id}/sequence_ids/id_list.txt"
	output:
		multi_fna = "{run}_{min_ident}/{experiment_id}/fasta/multi_fasta.fna"
	conda:
		"envs/ptn2locus.yaml"
	shell:
		"""
		touch {output.multi_fna}
		python3 py/ptn2locus.py {input.id_list} 'thedoudnalab@gmail.com' -o {output.multi_fna}
		"""
