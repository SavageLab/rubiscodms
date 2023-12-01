# Native modules
from argparse import ArgumentParser as argp
import copy
import re
# External modules
import pandas as pd
import yaml
from Bio import SeqIO
from bioservices import UniProt
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import urllib.request
import urllib.error


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='fetch_region',
		description='Uses RefSeq/Uniprot protein IDs to retrieve their associated genomic regions. Returns a multi FASTA nt',
		usage='%(prog)s [options]')
	# Define arguments
	parser.add_argument('id_list',
	                    help='Path to a single column list of RefSeq IDs containing a header line (column name)')  # positional argument
	parser.add_argument('login',
	                    help='Entrez login')
	parser.add_argument('-o',
	                    dest='output',
	                    default='sequences.fa',
	                    help='Path to output file [default=sequences.fa]')
	parser.add_argument('-c',
	                    dest='config',
	                    help='Path to config file in the dms_workflow.smk format. '
	                         'The program will save the output on the path specified by input_dir.')
	parser.add_argument('-d',
	                    dest='database',
	                    default='nuccore',
	                    help='NCBI database to search for entries associated with source PTN ID [default=nuccore')

	# positional argument

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def id_parsing(id_string, sep_dict):
	loop_id = id_string
	for level in sep_dict:
		try:
			loop_id = loop_id.split(sep_dict[level][0])[sep_dict[level][1]]
		except IndexError:
			continue
	return loop_id


def ukb2ncbi(uid):
	u = UniProt(verbose=False)
	# gbk_id = u.mapping("UniProtKB_AC-ID", "EMBL-GenBank-DDBJ", query=uid, polling_interval_seconds=3, max_waiting_time=100)["results"][0]["to"]
	try:
		prot_id = u.mapping("UniProtKB_AC-ID", "EMBL-GenBank-DDBJ_CDS", query=uid, polling_interval_seconds=3, max_waiting_time=100)["results"][0]["to"]
	except (TypeError, IndexError):
		prot_id = ''
	return prot_id


def elink_routine(db, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	link_record = ""
	server_attempts = 0
	handle = ""
	try:
		handle = Entrez.elink(dbfrom="protein", db=db, id=f"{hit_uid}")
	except urllib.error.HTTPError as err:
		if err.code == 500:
			print(f'An internal server error occurred while handling the accession {hit_uid}')
			not_found = hit_uid
			return linked, hit_uid, not_found
	try:
		link_record = Entrez.read(handle)
	except (RuntimeError, ValueError):
		not_found = hit_uid
	if link_record:
		try:
			linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
			if linked not in dup_check:
				dup_check.append(linked)
		except (IndexError, KeyError):
			not_found = hit_uid
	handle.close()
	return linked, hit_uid, not_found


def ptn_to_nuc(id_list, db_name):
	progress = 0
	nucleotide_uid_list = []
	source2target = {}
	not_found_list = []

	for seq_id in id_list:
		progress += 1
		# Standardize protein identifiers to NCBI UIDs through ESearch
		handle = Entrez.esearch(db="protein", term=f"{seq_id}", idtype="acc")
		search_record = Entrez.read(handle)
		try:
			uid = search_record['IdList'][0]
		except IndexError:
			print(f"Entrez.esearch could not find results for {seq_id}")
			not_found_list.append(seq_id)
			continue
		handle.close()

		# Grab Nuccore UIDs from the source database
		loop_nuc_gi, loop_nuc_acc, not_found_hit = elink_routine(db_name, uid)
		if not_found_hit:
			loop_nuc_gi, loop_nuc_acc, c_not_found_hit = elink_routine(db_name,
			                                                           ukb2ncbi(not_found_hit))

			if not loop_nuc_gi:
				not_found_list.append(c_not_found_hit)
				continue
		if loop_nuc_gi:
			source2target.setdefault(loop_nuc_gi, seq_id)
			nucleotide_uid_list.append(loop_nuc_gi)
			print(f"Nuccore GI and ptn accession: {loop_nuc_gi} {seq_id}")

	# Ouputs nuccore uids and the ptn->nuc uid links
	print(f"Total n of not found hits: {len(not_found_list)}")
	return nucleotide_uid_list, source2target, list(set(not_found_list))


def nuc_to_gb(uid_list):
	# Get Genbank records for each Nuccore UID
	gb_records = {}
	for uid in uid_list:
		handle = Entrez.efetch(db="nucleotide", id=f"{uid}", rettype = 'gbwithparts', retmode="text")
		record = SeqIO.read(handle, "genbank")
		gb_records.setdefault(uid, record)
		# Returns a list  of Genbank SeqRecords objects
	return gb_records


def gb_plier(uid_to_gb, uid_to_acc, win_size):
	gbk_target = {}
	prot_dict = {}
	uid_to_gb_dict = uid_to_gb.copy()
	for hit_uid in uid_to_gb_dict:
		gbk = uid_to_gb_dict[hit_uid]
		for seq_feature in gbk.features:
			# Avoid blank feature that may occur in GenBank entries
			try:
				qualifiers = seq_feature.qualifiers
			except AttributeError:
				continue
			# Restrict search to protein-containing features
			if "protein_id" in qualifiers:
				prot_id = qualifiers["protein_id"][0]
				# Search for the protein-ids of interest
				if re.search(prot_id, uid_to_acc[hit_uid]):
					# Process feature information for future ref
					f_start = seq_feature.location.start.real
					f_end = seq_feature.location.end.real
					f_strand = seq_feature.strand
					f_seq = qualifiers["translation"][0]
					f_len = len(f_seq)
					highlight_feature = copy.deepcopy(seq_feature)
					highlight_feature.type = "highlight"
					# Set start/end coords using window size
					start = max(int(min([f_start, f_end])) - win_size, 0)
					end = min(int(max([f_start, f_end])) + win_size, len(gbk.seq))

					# Account for strand orientation
					feature_seq_format = gbk.seq[start:end]
					if f_strand < 0:
						feature_seq_format = feature_seq_format.reverse_complement()

					# Create a SeqRecord object with the feature of interest
					gbk_focused = SeqRecord(
						id=gbk.id,
						annotations=gbk.annotations,
						dbxrefs=gbk.dbxrefs,
						seq=feature_seq_format,
						description=gbk.description
					)
					gbk_focused.features.append(highlight_feature)
					# Gather protein data for reference
					prep_prot_dict = {
					                  "nuccore_acc": gbk.id,
					                  # "region_seq": gbk.seq[start:end + 1],
					                  "window_start": start,
					                  "window_end": end,
					                  "feature_start": f_start,
					                  "feature_end": f_end,
					                  "strand": f_strand,
					                  "feature_len": f_len,
					                  "blastp_hit": prot_id,
					                  "ptn_sequence": f_seq
					                  }
					prot_dict.setdefault(uid_to_acc[hit_uid], prep_prot_dict)
					gbk_target.setdefault(f"{gbk.id}_{start}-{end}", gbk_focused)

	return gbk_target, prot_dict


def main():
	# Call argument parsing function
	args = parse_arguments()
	# Load values into variables
	input_col = args.id_list
	entrez_login = args.login
	db = args.database
	config = args.config
	output_path = args.output
	# Set blast arguments
	# DEBUG
	df = pd.read_csv(input_col)
	hit_list = df.iloc[:, 0].tolist()

	if config:
		with open(config, 'r') as config_handle:
			config_dms = yaml.safe_load(config_handle)
		output_path = f"{config_dms['input_dir']}/multi_fasta.fna"

	# Entrez authentication
	print(f"Entrez login: {entrez_login}")
	Entrez.email = entrez_login

	# Query NCBI to get nuccore UIDs associated with the protein hits using ESearch/ELink
	print("Linking protein hit ids to Nuccore entries")
	nuc_uid_list, hit_to_link, hits_not_found = ptn_to_nuc(hit_list, db)

	# Get GenBank entries through EFetch
	print(f"Retrieving GenBank objects. {len(nuc_uid_list)} in total")
	gb_seqrec = nuc_to_gb(nuc_uid_list)

	# Narrow down Genbank files based on hit UIDs
	print("Extracting relevant features from GenBank objects")
	target_gb_dict, target_hit_dict = gb_plier(gb_seqrec, hit_to_link, 0)

	for item in target_gb_dict:
		with open(output_path, "a") as f:
			target_gb_dict[item].id = target_gb_dict[item].features[0].qualifiers["protein_id"][0]
			SeqIO.write(target_gb_dict[item], f, "fasta")

	for entry in hits_not_found:
		with open("failed_hits.txt", 'a') as f:
			f.write(f"{entry}\n")


if __name__ == "__main__":
	main()
