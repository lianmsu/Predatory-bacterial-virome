
cd /lustre/home/lipengwei2024phd/06mobilome_predators/83phage_phylogeny

seqkit split -i /lustre/home/lipengwei2024phd/06mobilome_predators/83phage_phylogeny/all_48_phage.fasta -O ./split_48_phages

comparem aai_wf --cpus 10 /lustre/home/lipengwei2024phd/06mobilome_predators/83phage_phylogeny/split_48_phages aai_output --file_ext fasta