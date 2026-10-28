#!/bin/sh
#SBATCH -o jobdia1060.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J diamond
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mail-type=end
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
DB_ROOT="${DB_ROOT:-${LUSTRE_HOME}/database}"       # reference databases
# ---------------------------------------------------------------------------

source activate eggnog

orf="${LUSTRE_HOME}/06mobilome_predators/74prophage_filter/checkv_high_quality_genomadfiltered_provirus_ORF.faa"
out=${LUSTRE_HOME}/06mobilome_predators/74prophage_filter/diamond_res

#各种类型功能基因的数据库
NRdb=${DB_ROOT}/NR/nr.dmnd

query=${orf}

diamond blastp -k 1 -e 0.000001 -p $SLURM_NTASKS -d ${NRdb} -q ${query} -o ${out}/diamond_NR_out --id 60 --subject-cover 90

#--id 80 比较好
#The host system is detected to have 404 GB of RAM. It is recommended to increase the block size for better performance using these parameters : -b8 -c1
#--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#

#--id 80, 
echo "successfully"


