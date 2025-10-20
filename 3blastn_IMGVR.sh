#!/bin/bash
#SBATCH --job-name=blastnVR
#SBATCH --partition=fat
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=48
#SBATCH --output=job-virus.%j.%N.out

source activate blast2.16
cd /lustre/home/lipengwei2024phd/06mobilome_predators/70_spacer2IMGVR

mkdir -p ./2spacer2virus

blastn -db ./blastdb_IMGVR/IMGVR.db -query ./all_CRISPR_level4_spacers.fa -out ./2spacer2virus/spacers_virus.blastn.out -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' -num_threads 48 -mt_mode 1 -evalue 1e-5 -word_size 8 -task blastn-short

cat ./2spacer2virus/spacers_virus.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > spacers_virus.blastn.1mis95cov.out
#cat spacers_virus.blastn.1mis95cov.out | cut -f2 | sort | uniq > targeted.virus.id


