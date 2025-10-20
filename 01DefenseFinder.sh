#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J df2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48


# 开始计时
start_time=$(date +%s)

source activate defensefinder
#########################################################################

#########################################################################

genome_path=/lustre/home/lipengwei2024phd/06mobilome_predators/1genome_collections/all_4696_predators_genome_renamed.fa
defense_workdir=/lustre/home/lipengwei2024phd/06mobilome_predators/64_DF2
output_dir="${defense_workdir}/res"

mkdir -p ${output_dir}
cd ${defense_workdir}

defense-finder run -a --preserve-raw -o "$output_dir" ${genome_path}

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
