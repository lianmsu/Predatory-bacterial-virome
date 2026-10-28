# =============================================================================
# Predatory-bacterial-virome -- phage catalog / Plasmid metadata and host assignment
# =============================================================================
# Assemble the per-plasmid metadata (source, size, predicted host) of the
# non-redundant plasmid set.
#
# Manuscript      : Plasmid catalog (Fig. 1I)
# Main outputs    : self_selected_plasmid.tsv and related tables
# Original script : 1mobilome_predatory_bacteria/3vcontact2_plasmid/1_plasmid_metadata.R
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


setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/3vcontact2_plasmid"))

self_selected_plasmid <- read.delim2("seq_info_table.tsv",sep = "\t")

library(dplyr)

self_selected_plasmid$source <- "NCBI"

self_selected_plasmid$source[364:400] <- "This study"

self_selected_plasmid$topology <- "Uncomplete"


plasmid_genomad_metadata <- read.delim2("nr_plasmid_rep_seq_plasmid_summary.tsv")
# 复制数据框以避免修改原数据（可选）
plasmid_processed <- plasmid_genomad_metadata


# 首先将两个数据框通过seq_name和Node列关联
self_selected_plasmid <- self_selected_plasmid %>%
  left_join(plasmid_processed %>% select(seq_name, topology),
            by = c("seq_id" = "seq_name")) %>%
  mutate(topology = case_when(
    topology.y %in% c("DTR", "ITR") ~ "Complete",
    TRUE ~ topology.x  # 保持原来的值
  )) %>%
  select(-topology.x, -topology.y)  # 删除临时列

self_selected_plasmid$MGE <- "Unknown"
self_selected_plasmid <- self_selected_plasmid %>%
  left_join(plasmid_processed %>% select(seq_name, n_hallmarks),
            by = c("seq_id" = "seq_name")) %>%
  mutate(MGE = ifelse(n_hallmarks > 0, "Plasmid", MGE)) %>%
  select(-n_hallmarks)  # 删除临时列

self_selected_plasmid$MGE[is.na(self_selected_plasmid$MGE)] <- "Unknown"
self_selected_plasmid$MGE[1:363] <- "Plasmid"


#宿主信息####
IMGPR <- read.delim2("filtered_spacers_plasmid.blastn.1mis95cov.out",header = F)
devired_from_genome <- read.delim2("spacers_plasmid.blastn.1mis95cov.out",header = F)

library(stringr)
IMGPR$genome <- str_replace(IMGPR$V1, "-[^-]*$", "")
IMGPR_2 <- select(IMGPR,V2,genome)
IMGPR_2 <- IMGPR_2 %>%
  distinct(V2, .keep_all = TRUE)
colnames(IMGPR_2) <- c("Plasmid","Host_genome")

devired_from_genome$genome <- str_replace(devired_from_genome$V1, "-[^-]*$", "")
devired_from_genome_2 <- select(devired_from_genome,V2,genome)
devired_from_genome_2 <- devired_from_genome_2 %>%
  distinct(V2, .keep_all = TRUE)
colnames(devired_from_genome_2) <- c("Plasmid","Host_genome")

all_plasmid_host <- rbind(IMGPR_2, devired_from_genome_2)

all_4696genome_metadata <- read.delim2("all_4649_genome_metadata.tsv")

all_plasmid_host <- all_plasmid_host %>%
  left_join(all_4696genome_metadata %>% select(Genome, phylum),
            by = c("Host_genome" = "Genome"))

all_plasmid_host <- all_plasmid_host %>%
  distinct(Plasmid, .keep_all = TRUE)


self_selected_plasmid$host_phylum <- "Unknown"
self_selected_plasmid <- self_selected_plasmid %>%
  left_join(all_plasmid_host %>% select(Plasmid, phylum),
            by = c("seq_id" = "Plasmid")) %>%
  mutate(host_phylum = phylum) %>%
  select(-phylum)  # 删除临时列

self_selected_plasmid$host[364:400] <- self_selected_plasmid$host_phylum[364:400]
self_selected_plasmid <- self_selected_plasmid[,-6]

write.table(self_selected_plasmid,"self_selected_plasmid.tsv",sep = "\t",
            quote = F)






