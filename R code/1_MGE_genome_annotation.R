
setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/08_special_MGE")

all_mge_gff <- 
  read.delim2("selected_MGE.gff",header = F, sep = '\t')

library(dplyr)
library(stringr)
all_mge_gff <- all_mge_gff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列
all_mge_gff <- all_mge_gff[,c(1,3,4,5,7,10)]

colnames(all_mge_gff) <- c("contig","type","start","stop","strand",
                           "name")

#type和name可以改为需要的内容（type决定颜色，name决定标签）

write.table(all_mge_gff,"all_mge_gff.tsv",
            sep = "\t",row.names = F,
            quote = F)











