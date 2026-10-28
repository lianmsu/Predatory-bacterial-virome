#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J nr_blast6
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --output=job.blast6_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
DB_ROOT="${DB_ROOT:-${LUSTRE_HOME}/database}"       # reference databases
# ---------------------------------------------------------------------------

cd ${LUSTRE_HOME}/06mobilome_predators/62_NR_filtering
# 设置BLAST数据库路径
source activate blast2.16
export BLASTDB=${DB_ROOT}/NR/nr
# 蛋白序列使用blastp
blastp -query ${LUSTRE_HOME}/06mobilome_predators/62_NR_filtering/split_ORF_dir/all_CRISPR_level4_contig_ORF.part_006.faa \
       -db nr \
	   -task blastp-fast \
	   -taxids "2" \
       -out protein_vs_nr.blast_part_006 \
       -outfmt "6 std staxids stitle" \
       -evalue 1e-5 \
	   -mt_mode 1 \
       -num_threads 14 \
	   -max_target_seqs 1 \
	   -subject_besthit
	   
