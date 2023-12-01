#!/bin/bash
#SBATCH --job-name=phydms_trueRubs
#SBATCH --nodes=2
#SBATCH --thread-spec=48
#SBATCH --mem=0
# ACTIVATE ANACONDA
eval "$(conda shell.bash hook)"
conda activate blast2region
# shellcheck disable=SC2164
cd /groups/doudna/team_resources/toolbox/phylo_dms_workflow/
python py/ptn2locus.py /groups/doudna/projects/daniel_projects/prywes_n/input_data/rubisco_bonafide.csv "thedoudnalab@gmail.com" -c /groups/doudna/team_resources/toolbox/phylo_dms_workflow/config/prywes_dms.yaml
