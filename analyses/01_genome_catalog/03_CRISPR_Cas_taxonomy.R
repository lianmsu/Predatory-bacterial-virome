# =============================================================================
# Predatory-bacterial-virome -- genome catalog / CRISPR-Cas type and subtype assignment
# =============================================================================
# Assign every detected CRISPR-Cas system to its type/subtype from the
# CRISPRCasTyper output and tabulate it per genome.
#
# Manuscript      : Fig. 2A
# Main outputs    : per-genome CRISPR-Cas type table
# Original script : 1mobilome_predatory_bacteria/2CRISPRTyper/1_CRISPRCas_tax.R
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


setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/2CRISPRTyper"))

CRISPR_Cas <- read.delim2("CRISPR_Cas.tab")

library(stringr)
CRISPR_Cas$genome <- str_extract(CRISPR_Cas$Contig, "^[^-]*")

all_4649genome_metadata <- read.delim2("all_4649_genome_metadata.tsv")

CRISPR_Cas_genome <- all_4649genome_metadata[all_4649genome_metadata$Genome %in% CRISPR_Cas$genome,]

# 按照CRISPR_Cas$genome的顺序重新排序
CRISPR_Cas_genome2 <- CRISPR_Cas_genome[
  match(CRISPR_Cas$genome, CRISPR_Cas_genome$Genome),
]

CRISPR_Cas_genome2 <- merge(CRISPR_Cas_genome, CRISPR_Cas, by.x = "Genome",by.y = "genome")

#I-F挑选
table(CRISPR_Cas$Prediction)

I_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-A",]
I_B <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-B",]
I_C <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-C",]
I_D <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-D",]
I_E <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-E",]
I_F <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-F",]
I_G <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="I-G",]

II_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="II-A",]
II_B <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="II-B",]
II_C <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="II-C",]

III_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-A",]
III_B <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-B",]
III_C <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-C",]
III_D <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="III-D",]

V_A <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="V-A",]
V_B1 <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="V-B1",]
V_M <- CRISPR_Cas_genome2[CRISPR_Cas_genome2$Prediction=="V-M",]

Hybrid <- CRISPR_Cas_genome2[grepl("Hybrid", CRISPR_Cas_genome2$Prediction), ]
