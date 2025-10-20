#!/bin/bash
# 脚本功能：从diamond比对结果中通过NCBI Taxonomy数据库生成标准物种注释
cd /lustre/home/lipengwei2024phd/06mobilome_predators/74prophage_filter/diamond_res
# ---------- 参数设置 ----------
diamond_RESULT="diamond_NR_out"          # 您的BLAST比对结果文件
NR_TAXID_MAP="/lustre/home/lipengwei2024phd/database/NCBI_tax/prot.accession2taxid.FULL"  # NCBI蛋白acc2taxid文件
TAXDUMP_DIR="/lustre/home/lipengwei2024phd/database/NCBI_tax/"  # taxdump解压目录
OUTPUT="protein_annotation_besthit.tsv"     # 最终注释结果文件

# ---------- 步骤1：从BLAST结果提取best hit ----------
echo "Step 1: 提取BLAST best hit..."
awk '!seen[$1]++' "$diamond_RESULT" > best_hits.tmp

# 提取匹配的蛋白accession (从sseqid中提取如XP_001234.1部分)
awk -F"\t" '{
    split($2, arr, "|"); 
    if (arr[1]=="ref" || arr[1]=="gb" || arr[1]=="dbj") print arr[2]
    else print $2  # 如果格式不符合|分割，直接使用原字段
}' best_hits.tmp > nr_accessions.tmp

# ---------- 步骤2：获取TaxID ----------
echo "Step 2: 查询蛋白accession对应的TaxID..."
#awk 'NR==FNR{a[$1]=$0; next} $1 in a{print a[$1]}' \
#    "$NR_TAXID_MAP" nr_accessions.tmp > accession_taxid.tmp

# 仅加载 nr_accessions.tmp 中存在的 Accession 到内存
nohup awk 'NR==FNR{ids[$1]; next} $1 in ids' nr_accessions.tmp "$NR_TAXID_MAP" | \
    cut -f1,2 > accession_taxid.tmp 2>> nohup_accession2taxid2.log &

# ---------- 步骤3：获取完整分类信息 ----------
echo "Step 3: 获取TaxID对应的分类信息..."
taxonkit lineage --data-dir "$TAXDUMP_DIR" -i 2 accession_taxid.tmp | \
    taxonkit reformat -i 3 -f "{p}\t{c}\t{o}\t{f}\t{g}\t{s}" | \
    cut -f1,2,4,5,6,7,8,9 > taxid_lineage.tmp

# ---------- 步骤4：合并结果 ----------
echo "Step 4: 生成最终注释文件..."
awk 'BEGIN {FS=OFS="\t"} 
NR==FNR {
    # 读取第二个文件（taxid_lineage.tmp），以第一列作为键
    lineage[$1] = $0
    next
}
{
    # 处理第一个文件（best_hits.tmp）
    if ($2 in lineage) {
        # 如果第二列在第二个文件中找到匹配，输出合并结果
        print $0, lineage[$2]
    } else {
        # 如果没有匹配，设置为unknown
        print $0, "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown"
    }
}' taxid_lineage.tmp best_hits.tmp > merged_results.tsv

# ------统计每个contig的物种----------------
awk -F"\t" '
{
    # 提取contig名称
    contig = $1
    sub(/_[^_]*$/, "", contig)
    
    # 处理各分类级别信息

    phylum[contig][$15]++    # 第14列: Phylum
    class[contig][$16]++     # 第15列: Class
    order[contig][$17]++     # 第16列: Order
    family[contig][$18]++    # 第17列: Family
    genus[contig][$19]++     # 第18列: Genus
    species[contig][$20]++   # 第19列: Species
    
    total_genes[contig]++
}
END {
    # 输出标题行
    print "Contig\tTotal_Genes\tDominant_Phylum\tDominant_Class\tDominant_Order\tDominant_Family\tDominant_Genus\tDominant_Species"
    
    # 遍历每个contig
    for (contig in total_genes) {
        # 为每个级别找到最主要分类
        dom_phylum = find_dominant(phylum[contig])
        dom_class = find_dominant(class[contig])
        dom_order = find_dominant(order[contig])
        dom_family = find_dominant(family[contig])
        dom_genus = find_dominant(genus[contig])
        dom_species = find_dominant(species[contig])
        
        print contig "\t" total_genes[contig] "\t" dom_phylum "\t" dom_class "\t" dom_order "\t" dom_family "\t" dom_genus "\t" dom_species
    }
}

function find_dominant(taxa_array) {
    max_count = 0
    dominant = "Unknown"
    for (taxa in taxa_array) {
        if (taxa != "" && taxa != "unknown" && taxa_array[taxa] > max_count) {
            max_count = taxa_array[taxa]
            dominant = taxa
        }
    }
    return dominant
}' /lustre/home/lipengwei2024phd/06mobilome_predators/74prophage_filter/diamond_res/merged_results.tsv > contig_dominant_taxa.tsv

# --------手动整理contig_dominant_taxa.tsv -----------


# ---------- 清理临时文件 ----------
rm *.tmp

echo "完成！结果已保存到 $OUTPUT"
echo "前5条注释："
head -n 5 "$OUTPUT"
