#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J prodigal
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=job.01prodigal_lpw.%j.%N.out

source activate prokka
# 开始计时
start_time=$(date +%s)

cd /lustre/home/lipengwei2024phd/06mobilome_predators/4_4696genome_prodigal
mags=/lustre/home/lipengwei2024phd/06mobilome_predators/1genome_collections/all_4696_predators_genome_renamed.fa
res_path=/lustre/home/lipengwei2024phd/06mobilome_predators/4_4696genome_prodigal

echo "Running Prodigal"
prodigal -p meta -i "${mags}" \
		-f gff -o "${res_path}/all_4696_predators_genome_renamed_ORF.gff" \
		-m -d "${res_path}/all_4696_predators_genome_renamed_ORF.fna" \
		-a "${res_path}/all_4696_predators_genome_renamed_ORF.faa"




# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0