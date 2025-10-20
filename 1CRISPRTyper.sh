#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J finder
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=13
#SBATCH --output=job.01CASfinder_lpw.%j.%N.out

source activate cctyper
cd /lustre/home/lipengwei2024phd/06mobilome_predators/61_CRISPRTyper

cctyper /lustre/home/lipengwei2024phd/06mobilome_predators/61_CRISPRTyper/all_CRISPR_level4.fa /lustre/home/lipengwei2024phd/06mobilome_predators/61_CRISPRTyper/res -t 13 --prodigal meta --no_grid
