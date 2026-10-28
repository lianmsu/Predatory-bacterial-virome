# =============================================================================
# Predatory-bacterial-virome -- vcontact3 analysis.R / vContact3 analysis (inherited from the ESKAPE-virome project)
# =============================================================================
# vContact3 phage/provirus analysis of a phage set; the working directory
# still points at the ESKAPE-virome project (GTR_WORK/2ESKAPE_MGE), so the
# script is kept for reference only.
#
# Manuscript      : -
# Main outputs    : 1_vcontact3_analysis.RData
# Original script : 02_phage_vcontact3/1_vcontact3_analysis.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
GTR_WORK <- normalizePath(Sys.getenv("GTR_WORK", unset = "."), mustWork = FALSE)

setwd(file.path(GTR_WORK, "2ESKAPE_MGE/5phage_vcontact3"))

ESKAPE_all_phage_info <- read.delim2("ESKAPE_all_phage_info_filter.tsv")

ESKAPE_complete_phage_info <- 
  ESKAPE_all_phage_info[ESKAPE_all_phage_info$checkv_quality %in% 
                                                      c("Complete"),]

table(ESKAPE_complete_phage_info$Host_tax)
# Ab    Ef    En En,Kp    Kp    Pa    Sa 
# 138    29   102    47   443    99    35

multi_host_complete_phage <- 
  ESKAPE_complete_phage_info[ESKAPE_complete_phage_info$Host_tax %in% c("En,Kp"),]


vcontact3_result <- read.csv("part1.cyjs default node.csv")
length(unique(vcontact3_result$id))

potential_novel_complete_phage <- 
  ESKAPE_complete_phage_info[! ESKAPE_complete_phage_info$seq_name %in% vcontact3_result$id,]
head(potential_novel_complete_phage)

multi_host_complete_phage_novel <- 
  multi_host_complete_phage[!multi_host_complete_phage$seq_name %in% vcontact3_result$id,]

write.table(potential_novel_complete_phage, 
            "potential_novel_complete_phage_ID.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )

write.table(multi_host_complete_phage, 
            "multi_host_complete_phage_ID.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )
save.image("1_vcontact3_analysis.RData")

load("1_vcontact3_analysis.RData")




