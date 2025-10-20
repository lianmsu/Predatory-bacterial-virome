#!/bin/sh
#SBATCH --partition=fat
#SBATCH -J genomad
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH -o job.genomad_lpw.%j.%N.out

# 开始计时
start_time=$(date +%s)


source activate genomad
genomad_db=/lustre/home/lipengwei2024phd/database/genomad_db
#########################################################################
cd /lustre/home/lipengwei2024phd/06mobilome_predators/72genomad
MAG_contig=/lustre/home/lipengwei2024phd/06mobilome_predators/1genome_collections/all_4696_predators_genome_renamed.fa
workdir=/lustre/home/lipengwei2024phd/06mobilome_predators/72genomad

genomad end-to-end --cleanup ${MAG_contig} ${workdir}/res ${genomad_db}

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "脚本运行时间：$runtime 秒"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0



