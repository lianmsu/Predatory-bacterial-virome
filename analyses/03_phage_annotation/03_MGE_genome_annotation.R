# =============================================================================
# Predatory-bacterial-virome -- phage annotation / Gene map of a representative MGE
# =============================================================================
# Gene-arrow annotation of the representative mobile genetic element
# carrying cargo defense/secretion genes.
#
# Manuscript      : Fig. S2
# Main outputs    : MGE gene annotation figure
# Original script : 08_special_MGE/1_MGE_genome_annotation.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)


setwd(file.path(GTR_DATA, "08_special_MGE"))

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











