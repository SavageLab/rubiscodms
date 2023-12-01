from Bio import AlignIO
from sys import argv


def main():

	# DEBUG
	# seq_path = "/Users/bellieny/projects/team_resources/toolbox/phylo_dms_workflow/protein_gym/MSA_files/A0A1I9GEU1_NEIME_full_11-26-2021_b08.a2m"

	seq_path = argv[1]
	output_id_list = argv[2]

	with open(seq_path, 'r') as f:
		records = AlignIO.read(f, "fasta")

	str_out = ""
	for rec in records:
		seq_id = rec.id.split("_")[1]
		seq_id = seq_id.split("/")[0]
		str_out += f"{seq_id}\n"

	with open(output_id_list, 'w') as out:
		out.write(str_out)


if __name__ == "__main__":
	main()
