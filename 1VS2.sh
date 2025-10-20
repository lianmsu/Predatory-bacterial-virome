#!/bin/bash
#SBATCH --job-name=vs2_020
#SBATCH --partition=cpu
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=14
#SBATCH --output=job-vs2_020.%j.%N.out

source activate vs2


cd /lustre/home/lipengwei2024phd/06mobilome_predators/71VirSorter2


virsorter run --keep-original-seq -i /lustre/home/lipengwei2024phd/06mobilome_predators/71VirSorter2/split_results/all_4696_predators_genome_renamed.part_020.fa -w vs2-pass1_020 --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j 14 all

