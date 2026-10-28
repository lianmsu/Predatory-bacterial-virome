# =============================================================================
# Predatory-bacterial-virome -- marker phylogeny / DJR major-capsid protein identifiers
# =============================================================================
# Build the short-to-full identifier table of the DJR major capsid proteins
# used as the CLANS input.
#
# Manuscript      : Fig. S5
# Main outputs    : merged_DJRMCP_seqname table
# Original script : 1mobilome_predatory_bacteria/9special_protein/DJR_MCP/1_seqname_treatment.R
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

setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/9special_protein/DJR_MCP"))

seqname <- read.delim2("merged_DJRMCP_seqname2.txt",header = F)
colnames(seqname) <- c("shortid","fullid")

seqname <- seqname %>% 
  mutate(fullid = trimws(fullid))


network_file <- read.delim2("clans_net")

group_metadata <- readxl::read_xlsx("41564_2024_1839_MOESM13_ESM.xlsx")

library(dplyr)
# 首先将group_metadata转换为长格式，便于匹配
group_long <- group_metadata %>%
  select(Node_A, `Node_A group`) %>%
  rename(fullid = Node_A, group = `Node_A group`) %>%
  bind_rows(
    group_metadata %>%
      select(Node_B, `Node_B group`) %>%
      rename(fullid = Node_B, group = `Node_B group`)
  ) %>%
  distinct(fullid, group)  # 去除重复的匹配

# 如果同一个fullid有多个group，可以选择保留第一个或处理冲突
group_long <- group_long %>%
  group_by(fullid) %>%
  slice(1) %>%  # 如果有重复，保留第一个出现的group
  ungroup()

# 将group信息合并到seqname数据框
seqname_with_group <- seqname %>%
  left_join(group_long, by = "fullid")

# 查看结果
head(seqname_with_group)

write.table(seqname_with_group,"seqname_with_group.tsv",quote = F, sep = "\t")

