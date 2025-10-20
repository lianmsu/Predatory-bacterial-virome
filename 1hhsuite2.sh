#!/bin/bash
#SBATCH --partition=fat
#SBATCH -J hhsearch_job
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=112
#SBATCH -o job.hhsearch.%j.%N.out
#SBATCH -e job.hhsearch.%j.%N.err
source activate hhsuite3
# 设置变量
cd /lustre/home/lipengwei2024phd/06mobilome_predators/80HHsuit
input_fasta="/lustre/home/lipengwei2024phd/06mobilome_predators/80HHsuit/nr_plasmid_rep_seq.faa"
output_dir="/lustre/home/lipengwei2024phd/06mobilome_predators/80HHsuit/nr_plasmid_protein"
database="/lustre/home/lipengwei2024phd/database/HHsuit"
hhsuite_out="/lustre/home/lipengwei2024phd/06mobilome_predators/80HHsuit/nr_plasmid_hhsuite"

# 创建输出目录
mkdir -p "$output_dir" "$hhsuite_out"
cd "$output_dir"

# 记录开始时间
echo "任务开始时间: $(date)"
echo "输入文件: $input_fasta"
echo "输出目录: $output_dir"
echo "HHsuite输出目录: $hhsuite_out"

# 1. 分割多序列FASTA文件
echo "开始分割FASTA文件..."
awk -v output_dir="$output_dir" '
BEGIN { count = 0 }
/^>/ {
    if (count > 0) close(file)
    count++
    file = output_dir "/seq" count ".fasta"
    print > file
    next
}
{
    print > file
}' "$input_fasta"

echo "分割完成，共生成 $(ls -1 "$output_dir"/seq*.fasta | wc -l) 个序列文件"

# 2. 使用 xargs 并行处理
echo "开始使用 xargs 并行处理..."

# 定义处理函数
process_sequence() {
    local seq_file="$1"
    local base=$(basename "$seq_file" .fasta)
    
    echo "处理序列: $base - 开始时间: $(date)"
    
    # 保持原有的hhblits参数不变
    hhblits -i "$seq_file" \
        -o "$hhsuite_out/${base}.hhr" \
        -d "${database}/UniRef30_2023_02/UniRef30_2023_02" \
        -oa3m "$hhsuite_out/${base}.a3m" \
        -n 3 \
        -e 1e-6
    
    # 保持原有的hhsearch参数不变
    hhsearch -i "$hhsuite_out/${base}.a3m" \
        -d "${database}/NCBI_CD_v3.19/NCBI_CD" \
        -d "${database}/pdb70/pdb70" \
        -d "${database}/Pfam35/pfama" \
        -d "${database}/phrogs/phrogs_v4" \
        -d "${database}/scope70/scope70" \
        -d "${database}/uniprot_sprot_vir70/uniprot_sprot_vir70" \
        -Z 250 \
        -loc \
        -z 1 \
        -b 1 \
        -B 250 \
        -ssm 2 \
        -sc 1 \
        -seq 1 \
        -norealign \
        -maxres 32000 \
        -atab "$hhsuite_out/${base}.tsv" \
        -aliw 180
    
    # 处理结果文件
    if [ -f "$hhsuite_out/${base}.tsv" ]; then
        grep "^>" "$hhsuite_out/${base}.tsv" > "${hhsuite_out}/${base}.tmp" && \
        mv "${hhsuite_out}/${base}.tmp" "$hhsuite_out/${base}.tsv"
    fi
    
    echo "完成序列: $base - 结束时间: $(date)"
}

# 导出函数和变量
export -f process_sequence
export hhsuite_out database

# 使用 xargs 进行并行处理
find "$output_dir" -name "seq*.fasta" -print0 | \
xargs -0 -P 50 -I {} bash -c 'process_sequence "$@"' _ {}

# 3. 生成结果摘要
echo "生成结果摘要..."
total_files=$(find "$output_dir" -name "seq*.fasta" | wc -l)
processed_files=$(find "$hhsuite_out" -name "*.tsv" | wc -l)
hhr_files=$(find "$hhsuite_out" -name "*.hhr" | wc -l)
a3m_files=$(find "$hhsuite_out" -name "*.a3m" | wc -l)

echo "=== 处理结果摘要 ==="
echo "总序列文件数: $total_files"
echo "生成的TSV文件数: $processed_files"
echo "生成的HHR文件数: $hhr_files"
echo "生成的A3M文件数: $a3m_files"

# 检查是否有失败的任务
if [ "$total_files" -ne "$processed_files" ]; then
    echo "警告: 有些任务可能失败"
    echo "失败数量: $((total_files - processed_files))"
fi

# 记录结束时间
echo "任务结束时间: $(date)"
echo "所有处理完成!"