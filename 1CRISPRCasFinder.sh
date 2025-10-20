#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J finder
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=job.01CASfinder_lpw.%j.%N.out

source activate crisprcasfinder
cd /lustre/home/lipengwei2024phd/06mobilome_predators/6_CRISPRCasFinder/CRISPRCasFinder-master

perl CRISPRCasFinder.pl -in /lustre/home/lipengwei2024phd/06mobilome_predators/6_CRISPRCasFinder/minced_crispr_contig_renamed.fasta -levelMin 4 -q

#perl CRISPRCasFinder.pl -in /lustre/home/lipengwei2024phd/06mobilome_predators/6_CRISPRCasFinder/minced_crispr_contig.fasta -levelMin 4 -cas -meta -q -cpuMacSyFinder 18
