#!/bin/bash
#SBATCH --job-name=form2_3
#SBATCH --time=7:00:00
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate blast2region
# shellcheck disable=SC2164
cd /groups/doudna/projects/daniel_projects/rubiscodms
mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/input_data/form2_form3
python py/ptn2locus.py /groups/doudna/projects/daniel_projects/prywes_n/input_data/id_lists/form2_form3.txt "thedoudnalab@gmail.com" -o /groups/doudna/projects/daniel_projects/prywes_n/input_data/form2_form3