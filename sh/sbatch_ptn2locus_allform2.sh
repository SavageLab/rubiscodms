#!/bin/bash
#SBATCH --job-name=allform2
#SBATCH --time=6:00:00
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate blast2region
# shellcheck disable=SC2164
cd /groups/doudna/projects/daniel_projects/rubiscodms
mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/input_data/all_form_2
python py/ptn2locus.py /groups/doudna/projects/daniel_projects/prywes_n/input_data/id_lists/all_form_2.txt "thedoudnalab@gmail.com" -o /groups/doudna/projects/daniel_projects/prywes_n/input_data/all_form_2
