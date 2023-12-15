# **** Variables ****
configfile: "config/prywes_pgym_dms.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile protein_gym_dataprep.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# PSI-Blast the MSA input using a list of sequence databases
		expand("{run}_{min_ident}/{experiment_id}/psiblast_out/{experiment_id}.blastout",
			run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		# Export FASTA sequences of PSIBLAST Hits containing taxid labels
		expand("{run}_{min_ident}/{experiment_id}/psiblast_out/{experiment_id}_hits.fasta",
			run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		# Generate Krona plot
		expand("{run}_{min_ident}/{experiment_id}/krona/{experiment_id}_krona-plot.html",
			run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		# Re-align the initial MSA with the PSI-BLAST hit sequences
		expand("{run}_{min_ident}/{experiment_id}/clustalo/{experiment_id}.msa.fasta",
			run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/sequence_ids/id_list.txt",
			run=config["run"], experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/fasta/multi_fasta.fna",
			run=config["run"], experiment_id=config["reference_seq"],min_ident=config["min_ident"])

# noinspection SmkAvoidTabWhitespace
rule iterative_search:
	input:
		refseq = lambda wildcards: glob.glob("{in_dir}/fasta/{reference_seq}.fasta".format(
			in_dir=config['input_dir'],reference_seq=config["reference_seq"][wildcards.experiment_id]))
	output:
		psiblast_out = "{run}_{min_ident}/{experiment_id}/psiblast_out/{experiment_id}.blastout"
	params:
		db = config["shared_db_path"],
		custom_cols = config["blast_custom_cols"],
	threads: config["threads"]
	message:
		"""
Blasting query :\n {input.refseq}\n
Against database:\n {params.db}\nGenerating:\n {output.psiblast_out}
Wildcards: {wildcards}
		"""
	shell:
		"""
        psiblast \
        -query {input.refseq} \
        -db {params.db} \
        -outfmt "10 {params.custom_cols}" \
        -num_threads {threads} \
        -num_iterations 10 \
        -max_hsps 1 \
        -subject_besthit \
        -gapopen 9 \
        -inclusion_ethresh 1e-10 \
        -evalue 1e-5 \
        -qcov_hsp_perc 70 \
        -out {output.psiblast_out}
        """

# noinspection SmkAvoidTabWhitespace
rule taxid_parse:
	input:
		psiblast_out = "{run}_{min_ident}/{experiment_id}/psiblast_out/{experiment_id}.blastout"
	output:
		taxid_counts = "{run}_{min_ident}/{experiment_id}/krona/{experiment_id}_taxid_counts.tsv",
		hits_fasta= "{run}_{min_ident}/{experiment_id}/psiblast_out/{experiment_id}_hits.fasta"
	params:
		blast_col_names = config["blast_custom_cols"]
	conda:
		"envs/bio.yaml"
	script:
		"py/blastout_taxid_count.py"


# noinspection SmkAvoidTabWhitespace
rule krona:
	input:
		taxid_counts = "{run}_{min_ident}/{experiment_id}/krona/{experiment_id}_taxid_counts.tsv"
	output:
		krona_chart = "{run}_{min_ident}/{experiment_id}/krona/{experiment_id}_krona-plot.html"
	# params:
	# 	taxdump_path = config["taxdump_path"]
	conda:
		"envs/krona.yaml"
	message:
		"""
This rule implements Krona for interactive visualization of taxonomic distribution of the sequences.
The first five lines of the shell section are intended to correct Krona's issues with handling of taxonomic
background data.
Input data: {input.taxid_counts}
Output: {output.krona_chart}
Wildcards used in this rule: {wildcards}
		"""
	shell:
		"""		
		mkdir $CONDA_PREFIX/bin/scripts || true
		mkdir $CONDA_PREFIX/bin/taxonomy || true
		ln -s $CONDA_PREFIX/opt/krona/scripts/extractTaxonomy.pl $CONDA_PREFIX/bin/scripts || true 
		ln -s $CONDA_PREFIX/opt/krona/scripts/taxonomy.make $CONDA_PREFIX/bin/scripts || true
		ktUpdateTaxonomy.sh
		ktImportTaxonomy -m 2 -t 1 -tax $CONDA_PREFIX/bin/taxonomy -o {output.krona_chart} {input.taxid_counts}
		"""


# noinspection SmkAvoidTabWhitespace
rule realignment:
	input:
		hits_fasta= "{run}_{min_ident}/{experiment_id}/psiblast_out/{experiment_id}_hits.fasta",
		refseq = lambda wildcards: glob.glob("{in_dir}/fasta/{reference_seq}.fasta".format(
			in_dir=config['input_dir'],reference_seq=config["reference_seq"][wildcards.experiment_id]))
	output:
		post_search_msa = "{run}_{min_ident}/{experiment_id}/clustalo/{experiment_id}.msa.fasta"
	params:
		merged_input = "{run}_{min_ident}/{experiment_id}/clustalo/{experiment_id}_merged-input.fasta",
	threads: config["threads"]
	shell:
		"""
		cat {input.refseq} {input.hits_fasta} > {params.merged_input}
		clustalo --iter 10 --threads {threads} -i {params.merged_input} -o {output.post_search_msa} -v
		"""

# noinspection SmkAvoidTabWhitespace
rule fetch_msa_ids:
	input:
		post_search_msa = "{run}_{min_ident}/{experiment_id}/clustalo/{experiment_id}.msa.fasta"
	output:
		id_list = "{run}_{min_ident}/{experiment_id}/sequence_ids/id_list.txt"
	conda:
		"envs/dms.yaml"
	message:
		"""
Fetch sequence IDs from {input.post_search_msa}
Export ID list to {output.id_list}	
		"""
	shell:
		"""
		python3 py/id_from_msa.py {input.post_search_msa} {output.id_list}
		"""

# noinspection SmkAvoidTabWhitespace
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
