cd /lustre/home/lipengwei2024phd/06mobilome_predators/70_spacer2IMGVR

cp /lustre/home/lipengwei2024phd/06mobilome_predators/62_NR_filtering/blastp_res_tax/contig_dominant_taxa_manuscript_filter.tsv ./

#过滤文件spacers_virus.blastn.1mis95cov.out，只保留该文件中包含文件contig_dominant_taxa_manuscript_filter.tsv第一列字符串的行，不用严格匹配。
# 提取contig_dominant_taxa_manuscript_filter.tsv的第一列内容
cut -f1 contig_dominant_taxa_manuscript_filter.tsv | tail -n +2 > filter_patterns.txt

# 使用grep进行包含匹配
grep -f filter_patterns.txt spacers_virus.blastn.1mis95cov.out > filtered_spacers_virus.blastn.1mis95cov.out

# 使用grep进行包含匹配
grep -f filter_patterns.txt spacers_plasmid.blastn.1mis95cov.out > filtered_spacers_plasmid.blastn.1mis95cov.out

cut -f2 filtered_spacers_virus.blastn.1mis95cov.out | sort | uniq | seqkit grep -f - /lustre/home/lipengwei2024phd/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna > filtered_spacers_virus.blastn.1mis95cov.fa

cut -f2 filtered_spacers_plasmid.blastn.1mis95cov.out | sort | uniq | seqkit grep -f - /lustre/home/lipengwei2024phd/06mobilome_predators/2IMG_PR/IMGPR_nucl.fna > filtered_spacers_plasmid.blastn.1mis95cov.fa

