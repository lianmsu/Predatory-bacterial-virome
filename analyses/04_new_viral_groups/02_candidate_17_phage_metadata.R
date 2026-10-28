# =============================================================================
# Predatory-bacterial-virome -- new viral groups / Metadata of the 17 candidate phages
# =============================================================================
# Combine the 17 candidate phages with the geNomad virus summary to obtain
# their genome size, topology and completeness.
#
# Manuscript      : Fig. 4
# Main outputs    : 17-phage metadata table
# Original script : 1mobilome_predatory_bacteria/6phage_phylo/1_all_17_phage_metadata.R
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

setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/6phage_phylo"))

my_17_phage <- read.delim2("Candidate_17_new_phages.tsv")


phage_genomad <- read.delim2("nr_phage_rep_seq_virus_summary.tsv")

library(dplyr)
my_17_phage <- my_17_phage %>% 
  left_join(phage_genomad,by=c("seq_name"="seq_name"))





