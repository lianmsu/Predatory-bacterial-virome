# =============================================================================
# Predatory-bacterial-virome -- phage catalog / Phage metadata and host assignment
# =============================================================================
# Assemble the per-phage metadata from the vContact2 run, including the
# predicted host genome of each phage.
#
# Manuscript      : Fig. 3C
# Main outputs    : self_selected_phage_node.tsv and related tables
# Original script : 1mobilome_predatory_bacteria/3vcontact2/1_phage_metadata.R
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


setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/3vcontact2"))

self_selected_phage <- read.delim2("self_prepared_phage.csv",sep = ",")

library(dplyr)
self_selected_phage_node <- unique(self_selected_phage$contig_id) %>%
  as.data.frame() %>% setNames("Node")

self_selected_phage_node$source <- "This_study"

self_selected_phage_node$source[1984:1992] <- "Other_study"

self_selected_phage_node$topology <- "Uncomplete"


phage_genomad_metadata <- read.delim2("nr_phage_rep_seq_virus_summary.tsv")
# 复制数据框以避免修改原数据（可选）
phage_processed <- phage_genomad_metadata
# 对topology为"Provirus"的行进行处理
provirus_mask <- phage_processed$topology == "Provirus"
phage_processed$seq_name[provirus_mask] <- sub("\\|[^|]*$", "", phage_processed$seq_name[provirus_mask])
phage_processed <- phage_processed %>%
  distinct(seq_name, .keep_all = TRUE)


# 首先将两个数据框通过seq_name和Node列关联
self_selected_phage_node <- self_selected_phage_node %>%
  left_join(phage_processed %>% select(seq_name, topology),
            by = c("Node" = "seq_name")) %>%
  mutate(topology = case_when(
    topology.y %in% c("DTR", "ITR") ~ "Complete",
    TRUE ~ topology.x  # 保持原来的值
  )) %>%
  select(-topology.x, -topology.y)  # 删除临时列

self_selected_phage_node$MGE <- "Unknown"
self_selected_phage_node <- self_selected_phage_node %>%
  left_join(phage_processed %>% select(seq_name, n_hallmarks),
            by = c("Node" = "seq_name")) %>%
  mutate(MGE = ifelse(n_hallmarks > 0, "Phage", MGE)) %>%
  select(-n_hallmarks)  # 删除临时列

self_selected_phage_node$MGE[is.na(self_selected_phage_node$MGE)] <- "Unknown"


#宿主信息####
IMGVR <- read.delim2("filtered_spacers_virus.blastn.1mis95cov.out",header = F)
devired_from_genome <- read.delim2("spacers_virus.blastn.1mis95cov.out",header = F)
provirus <- read.delim2("filtered_validated_proviruses.tsv")

library(stringr)
IMGVR$genome <- str_replace(IMGVR$V1, "-[^-]*$", "")
IMGVR_2 <- select(IMGVR,V2,genome)
IMGVR_2 <- IMGVR_2 %>%
  distinct(V2, .keep_all = TRUE)
colnames(IMGVR_2) <- c("Phage","Host_genome")

devired_from_genome$genome <- str_replace(devired_from_genome$V1, "-[^-]*$", "")
devired_from_genome_2 <- select(devired_from_genome,V2,genome)
devired_from_genome_2 <- devired_from_genome_2 %>%
  distinct(V2, .keep_all = TRUE)
colnames(devired_from_genome_2) <- c("Phage","Host_genome")

provirus$genome <- str_replace(provirus$source_seq, "-[^-]*$", "")
provirus_2 <- select(provirus, seq_name, genome)
colnames(provirus_2) <- c("Phage","Host_genome")


all_phage_host <- rbind(IMGVR_2, devired_from_genome_2, provirus_2)

all_4696genome_metadata <- read.delim2("all_4649_genome_metadata.tsv")

all_phage_host <- all_phage_host %>%
  left_join(all_4696genome_metadata %>% select(Genome, phylum),
            by = c("Host_genome" = "Genome"))

all_phage_host <- all_phage_host %>%
  distinct(Phage, .keep_all = TRUE)


self_selected_phage_node$host_phylum <- "Unknown"
self_selected_phage_node <- self_selected_phage_node %>%
  left_join(all_phage_host %>% select(Phage, phylum),
            by = c("Node" = "Phage")) %>%
  mutate(host_phylum = phylum) %>%
  select(-phylum)  # 删除临时列


self_selected_phage_node$host_phylum[1984:1992] <- c("Bdellovibrionota",
                                                     "Myxococcota",
                                                     "Myxococcota",
                                                     "Myxococcota",
                                                     "Myxococcota",
                                                     "Myxococcota",
                                                     "Bdellovibrionota",
                                                     "Bdellovibrionota",
                                                     "Bdellovibrionota")
self_selected_phage_node$topology[1984:1992] <- "Complete"
self_selected_phage_node$MGE[1984:1992] <- "Phage"

write.table(self_selected_phage_node,"self_selected_phage_node.tsv",sep = "\t",
            quote = F)






