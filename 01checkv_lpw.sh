#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J cv
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end


#进入checkv并导入数据库目录
source activate checkv
export CHECKVDB=/lustre/home/lipengwei2024phd/database/checkv_db/checkv-db-v1.5

virus_contig=/lustre/home/lipengwei2024phd/06mobilome_predators/72genomad/res/all_4696_predators_genome_renamed_summary/all_4696_predators_genome_renamed_virus.fna
workdir=/lustre/home/lipengwei2024phd/06mobilome_predators/73checkV

#end-to-end：污染度、完整度分析和去除原噬菌体污染
checkv end_to_end ${virus_contig} ${workdir} -t $SLURM_NTASKS

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"


