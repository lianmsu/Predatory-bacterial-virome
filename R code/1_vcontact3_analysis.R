setwd("D:/pku/creative_work/2ESKAPE_MGE/5phage_vcontact3")

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




