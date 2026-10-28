#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J finder
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=13
#SBATCH --output=job.01CASfinder_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
# ---------------------------------------------------------------------------

source activate cctyper
cd ${LUSTRE_HOME}/06mobilome_predators/61_CRISPRTyper

cctyper ${LUSTRE_HOME}/06mobilome_predators/61_CRISPRTyper/all_CRISPR_level4.fa ${LUSTRE_HOME}/06mobilome_predators/61_CRISPRTyper/res -t 13 --prodigal meta --no_grid
