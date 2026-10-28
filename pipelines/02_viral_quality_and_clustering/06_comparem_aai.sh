# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
# ---------------------------------------------------------------------------

cd ${LUSTRE_HOME}/06mobilome_predators/83phage_phylogeny

seqkit split -i ${LUSTRE_HOME}/06mobilome_predators/83phage_phylogeny/all_48_phage.fasta -O ./split_48_phages

comparem aai_wf --cpus 10 ${LUSTRE_HOME}/06mobilome_predators/83phage_phylogeny/split_48_phages aai_output --file_ext fasta
