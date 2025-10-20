cd /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP

cat ./*.fa > merged_DJRMCP.fa

mkdir -p clans_blast_res

cd /lustre/home/lipengwei2024phd/software/CLANS
python -m clans -nogui -infile /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP.fa -saveto /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/clans_blast_res -output_format delimited -eval 0.0001

python -m clans -nogui -infile /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP.fa -saveto /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/clans_blast_res2 -eval 0.0001

python -m clans -nogui -load /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/clans_blast_res2 -dorounds 3 -saveto /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/clans_net -output_format delimited

cd /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP
grep '>' /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP.fa > merged_DJRMCP_seqname.txt

sed -i 's/^>//' /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP_seqname.txt

awk '{
    prefix = substr($0, 1, 15)
    print prefix "_" (NR-1) "\t" $0
}' /lustre/home/lipengwei2024phd/06mobilome_predators/86special_protein/DJR_MCP/merged_DJRMCP_seqname.txt > output.txt
