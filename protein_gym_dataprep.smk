# **** Variables ****
configfile: "config/prywes_pgym_dms.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile protein_gym_dataprep.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}/{experiment_id}/sequence_ids/id_list.txt", run=config["run"], experiment_id=config["experiment_id"]),
		expand("{run}/{experiment_id}/fasta/multi_fasta.fna", run=config["run"], experiment_id=config["experiment_id"])

# noinspection SmkAvoidTabWhitespace
rule fetch_msa_ids:
	input:
		msa = lambda wildcards: glob.glob("{input_dir}/{run}/msa/{experiment_prefix}.a2m".format(
			input_dir="/groups/doudna/projects/daniel_projects/prywes_n/input_data",
			run=wildcards.run, experiment_prefix=config["experiment_id"][wildcards.experiment_id]))
	output:
		id_list = "{run}/{experiment_id}/sequence_ids/id_list.txt"
	shell:
		"""
		python3 id_from_msa.py {input.msa} {output.id_list}
		"""
rule protein_to_locus:
	input:
		id_list = "{run}/{experiment_id}/sequence_ids/id_list.txt"
	output:
		multi_fna = "{run}/{experiment_id}/fasta/multi_fasta.fna"
	conda:
		"envs/ptn2locus.yaml"
	shell:
		"""
		python3 py/ptn2locus.py {input.id_list} 'thedoudnalab@gmail.com' -o {output.multi_fna}
		"""
