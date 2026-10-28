#!/bin/bash
#SBATCH --job-name=vs2_020
#SBATCH --partition=cpu
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=14
#SBATCH --output=job-vs2_020.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
# ---------------------------------------------------------------------------

source activate vs2


cd ${LUSTRE_HOME}/06mobilome_predators/71VirSorter2


virsorter run --keep-original-seq -i ${LUSTRE_HOME}/06mobilome_predators/71VirSorter2/split_results/all_4696_predators_genome_renamed.part_020.fa -w vs2-pass1_020 --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 14 all

