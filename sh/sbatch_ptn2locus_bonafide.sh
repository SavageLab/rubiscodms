#!/bin/bash
#SBATCH --job-name=bonafide_rubs
#SBATCH --time=12:00:00
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate blast2region
# shellcheck disable=SC2164
cd /groups/doudna/projects/daniel_projects/rubiscodms
mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/input_data/all_bona_fide
python py/ptn2locus.py /groups/doudna/projects/daniel_projects/prywes_n/input_data/id_lists/all_bona_fide.txt "thedoudnalab@gmail.com" -o /groups/doudna/projects/daniel_projects/prywes_n/input_data/all_bona_fide
