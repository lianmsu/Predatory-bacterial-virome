# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
# ---------------------------------------------------------------------------
cd ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP

cat ./*.fa > merged_DJRMCP.fa

mkdir -p clans_blast_res

cd ${LUSTRE_HOME}/software/CLANS
python -m clans -nogui -infile ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP.fa -saveto ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/clans_blast_res -output_format delimited -eval 0.0001

python -m clans -nogui -infile ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP.fa -saveto ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/clans_blast_res2 -eval 0.0001

python -m clans -nogui -load ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/clans_blast_res2 -dorounds 3 -saveto ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/clans_net -output_format delimited

cd ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP
grep '>' ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP.fa > merged_DJRMCP_seqname.txt

sed -i 's/^>//' ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP_seqname.txt

awk '{
    prefix = substr($0, 1, 15)
    print prefix "_" (NR-1) "\t" $0
}' ${LUSTRE_HOME}/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP_seqname.txt > output.txt
