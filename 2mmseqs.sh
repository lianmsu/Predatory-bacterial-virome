#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J mmseq2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --output=job.01CASfinder_lpw.%j.%N.out

source activate mmseqs2
cd /lustre/home/lipengwei2024phd/06mobilome_predators/76_mmseqs2
# 创建输出目录

mkdir -p mmseqs2_clustering_phage

# 使用 easy-cluster 进行聚类
mmseqs easy-cluster ./all_phage.fa nr_phage mmseqs2_clustering_phage \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --cluster-mode 1 \
    --threads 48 \
    -e 1e-5
# 95% 序列一致性
# 90% 覆盖度（AF）
# 覆盖度模式：90% of bothe q and t
# 贪婪聚类
# 线程数
# e-value 阈值

mkdir -p mmseqs2_clustering_plasmid
mmseqs easy-cluster ./all_plasmid.fa nr_plasmid mmseqs2_clustering_plasmid \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --cluster-mode 1 \
    --threads 48 \
    -e 1e-5
# 95% 序列一致性
# 90% 覆盖度（AF）
# 覆盖度模式：90% of bothe q and t
# 贪婪聚类
# 线程数
# e-value 阈值
