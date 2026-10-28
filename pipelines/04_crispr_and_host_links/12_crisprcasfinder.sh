#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J finder
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=job.01CASfinder_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
# ---------------------------------------------------------------------------

source activate crisprcasfinder
cd ${LUSTRE_HOME}/06mobilome_predators/6_CRISPRCasFinder/CRISPRCasFinder-master

perl CRISPRCasFinder.pl -in ${LUSTRE_HOME}/06mobilome_predators/6_CRISPRCasFinder/minced_crispr_contig_renamed.fasta -levelMin 4 -q

#perl CRISPRCasFinder.pl -in ${LUSTRE_HOME}/06mobilome_predators/6_CRISPRCasFinder/minced_crispr_contig.fasta -levelMin 4 -cas -meta -q -cpuMacSyFinder 18
