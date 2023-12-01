# PhyDMS Workflow Utilizing DMS data

## Summary
+ This document describes the application of `phydms` to carry out phylogenetic analyses using Experimentally Informed Codon Models (ExpCM) 
### Materials
+ The starting point of this workflow was a list of RefSeq Protein IDs contained in the `/groups/doudna/projects/daniel_projects/prywes_n/input_data/bacterialFormIIs.csv` file
+ These sequences were the basis to build the nucleotide sequence alignment and subsequent phylogenetic reconstruction to support the `phydms` analysis
+ The reference sequence used in this workflow was [WP_011390153.1](https://www.ncbi.nlm.nih.gov/protein/WP_011390153.1?report=genbank&log$=protalign&blast_rank=2&RID=KJP5J24C016)

### Methods
+ The RefSeq IDs were used in a programmatic Entrez search to gather both amino acid and respective encoding nucleotide sequences of each entry.
+ A total of 522 were used to generate the nucleotide alignment using MAFFT (https://doi.org/10.1093/nar/gki198) via `phydms_prepalignment` from the `phydms` package 
+ The `phydms`(https://doi.org/10.7717/peerj.3657) was applied to the set of sequences

### Workflow 
+ Two major branches were conducted in this workflow:
  + One to inspect the Rubisco sequence evolution patterns using different settings
  + Antoher to look at a broader scope what patterns arise in different proteins using the Protein Gym data sets
+ RUBISCO
  + The scope of this analysis is based on the list of protein IDs in the spreadsheet`scratch/Supplemental_File_3_Rubisco_sequence_table_updated.csv`
  + Using the `notebooks/narrow_down_IdList.ipynb` notebook the appropriate list of IDs is filtered and stored in plan text
  + The list generated above goes into `py/ptn2locus.py` that will generate a multi_fasta.fna input for `dms_workflow.smk`
+ PROTEIN GYM