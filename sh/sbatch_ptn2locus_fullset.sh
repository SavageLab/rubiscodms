#!/bin/bash
#SBATCH --job-name=phydms_trueRubs
#SBATCH --nodes=2
#SBATCH --thread-spec=48
#SBATCH --mem=0
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate blast2region
# shellcheck disable=SC2164
cd /groups/doudna/projects/daniel_projects/rubiscodms
mkdir -p /groups/doudna/projects/daniel_projects/prywes_n/input_data/fullset
python py/ptn2locus.py /groups/doudna/projects/daniel_projects/prywes_n/input_data/id_lists/fullset.txt "thedoudnalab@gmail.com" -o /groups/doudna/projects/daniel_projects/prywes_n/input_data/fullset
