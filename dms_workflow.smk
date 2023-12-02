# **** Variables ****
# configfile: "config/prywes_dms.yaml"
configfile: "config/prywes_pgym_dms.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile dms_workflow.smk -j 5 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		expand("{run}_{min_ident}/{experiment_id}/processed_inputs/aa_preference.csv",run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/figures/aa_preference.pdf", run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/processed_inputs/nt_alignment_msa.fna", run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/rax_tree/RAxML_bestTree.nt_tree.newick", run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"]),
		expand("{run}_{min_ident}/{experiment_id}/phydmsresults/modelcomparison.md", run=config["run"],experiment_id=config["reference_seq"],min_ident=config["min_ident"])

# noinspection SmkAvoidTabWhitespace
rule convert_enrichment:
	input:
		dms_data = lambda wildcards: glob.glob("{in_dir}/{experiment_id}.csv".format(
			in_dir=config['input_dir'],
			experiment_id=config["reference_seq"]))
	output:
		dms_prefs = "{run}_{min_ident}/{experiment_id}/processed_inputs/aa_preference.csv"
	params:
		position_col = config["position_col"],
		enrichment_col = config["enrichment_col"],
		aminoacid_col = config["aminoacid_col"]
	conda:
		"envs/dms.yaml"
	message:
		"""
		Convert DMS scores in:\n {input.dms_data}
		To AA Preferences in:\n {output.dms_prefs}
		"""
	script:
		"py/enrichm2aa-preference.py"

# noinspection SmkAvoidTabWhitespace
rule inspect_prefs:
	input:
		dms_prefs = "{run}_{min_ident}/{experiment_id}/processed_inputs/aa_preference.csv",
	output:
		radial_prefs = "{run}_{min_ident}/{experiment_id}/figures/aa_preference.png"
	params:
		site_offset = config["site_offset"]
	conda:
		"envs/dms.yaml"
	script:
		"py/inspect_prefs.py"

# noinspection SmkAvoidTabWhitespace
rule aa_preference_logo:
	input:
		dms_prefs = "{run}_{min_ident}/{experiment_id}/processed_inputs/aa_preference.csv"
	output:
		aa_logo = "{run}_{min_ident}/{experiment_id}/figures/aa_preference.pdf"
	conda:
		"envs/dms.yaml"
	shell:
		"""
		phydms_logoplot --prefs {input.dms_prefs} {output.aa_logo}
		"""

# noinspection SmkAvoidTabWhitespace
rule seq_alignment:
	input:
		multi_fasta = lambda wildcards: glob.glob("{in_dir}/multi_fasta.fna".format(in_dir=config['input_dir']))
	output:
		msa = "{run}_{min_ident}/{experiment_id}/processed_inputs/nt_alignment_msa.fna"
	params:
		reference_seq = lambda wildcards: config["reference_seq"][wildcards.experiment_id],
		impute_ref = lambda wildcards: glob.glob("{in_dir}/fasta/{reference_seq}.fasta".format(in_dir=config['input_dir'],
			reference_seq=config["reference_seq"][wildcards.experiment_id])),
		ref_imputed = "{run}_{min_ident}/{experiment_id}/processed_inputs/multi_fasta_ref_imputed.fna"
	conda:
		"envs/dms.yaml"
	message:
		"""
		Align sequences in:\n {input.multi_fasta}
		Output MSA to:\n {output.msa}
		Use Reference Sequence: {params.reference_seq}
		"""
	shell:
		"""
		cat {params.impute_ref} {input.multi_fasta} > {params.ref_imputed}
		phydms_prepalignment {params.ref_imputed} {output.msa} {params.reference_seq} --minidentity {wildcards.min_ident}
		"""

# noinspection SmkAvoidTabWhitespace
rule infer_tree:
	input:
		msa = "{run}_{min_ident}/{experiment_id}/processed_inputs/nt_alignment_msa.fna"
	output:
		phylo_tree = "{run}_{min_ident}/{experiment_id}/rax_tree/RAxML_bestTree.nt_tree.newick"
	params:
		output_path = "{run}_{min_ident}/{experiment_id}/rax_tree"
	conda:
		"envs/rax.yaml"
	threads:
		config["threads"]
	shell:
		"""
		rm -rf {params.output_path}/* || true
		raxmlHPC-PTHREADS-SSE3 -s {input.msa} -w {params.output_path} -n nt_tree.newick -m 'GTRCAT' -p1 -T {threads} 
		"""

# noinspection SmkAvoidTabWhitespace
rule phydms:
	input:
		phylo_tree = "{run}_{min_ident}/{experiment_id}/rax_tree/RAxML_bestTree.nt_tree.newick",
		msa= "{run}_{min_ident}/{experiment_id}/processed_inputs/nt_alignment_msa.fna",
		dms_prefs = "{run}_{min_ident}/{experiment_id}/processed_inputs/aa_preference.csv"
	output:
		dms_models = "{run}_{min_ident}/{experiment_id}/phydmsresults/modelcomparison.md"
	params:
		outdir = "{run}_{min_ident}/{experiment_id}/phydmsresults",
		outprefix = "phydms_run_"
	conda:
		"envs/rax.yaml"
	threads:
		config["threads"]
	shell:
		"phydms_comprehensive --omegabysite --diffprefsbysite --tree {input.phylo_tree} --ncpus {threads} {params.outdir}/{params.outprefix} {input.msa} {input.dms_prefs}"
