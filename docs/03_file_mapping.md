# File mapping

Every script of this repository, its path in the previous (flat) layout and its
purpose. The scripts were **moved, renamed and made path-portable**; the analysis logic
is unchanged -- only the path statements and the newly added headers differ from the
original files.

Code edits applied to every file:

1. the site-specific absolute paths were replaced by variables -- `GTR_DATA`,
   `GTR_WORK` and `GTR_DEF` in the R scripts, `LUSTRE_HOME` / `DB_ROOT` in the pipeline
   scripts (see [below](#site-specific-paths-and-how-they-are-hidden));
2. a header was added to every R script (purpose, manuscript figure, outputs, original
   path);
3. line endings were normalised to LF (`.gitattributes`) and R/python/job-log artefacts
   were added to `.gitignore`.

## Upstream pipelines

| New path | Previous path | Purpose |
| --- | --- | --- |
| `pipelines/01_genome_and_phage_prediction/01_prodigal_genomes.sh` | `01prodigal_2.sh` | Prodigal gene prediction (meta mode) of the 4,696 predator genomes, producing the ORF GFF/FNA/FAA used by every downstream annotation. |
| `pipelines/01_genome_and_phage_prediction/02_genomad_genomes.sh` | `01genomad_lpw.sh` | geNomad end-to-end on the 4,696 predator genomes: plasmid and virus/provirus detection. |
| `pipelines/01_genome_and_phage_prediction/03_virsorter2.sh` | `1VS2.sh` | VirSorter2 on split genome files as a second viral-sequence caller. |
| `pipelines/02_viral_quality_and_clustering/04_checkv.sh` | `01checkv_lpw.sh` | CheckV quality assessment of the geNomad viral contigs (completeness, contamination, provirus removal). |
| `pipelines/02_viral_quality_and_clustering/05_mmseqs_cluster_phages.sh` | `2mmseqs.sh` | mmseqs2 easy-cluster of all phage genomes (95% nucleotide identity, 90% coverage) to build the non-redundant phage set. |
| `pipelines/02_viral_quality_and_clustering/06_comparem_aai.sh` | `1comparem.sh` | CompareM AAI workflow on the 48 phages (each genome split and compared pairwise) for the AAI heatmap. |
| `pipelines/03_functional_annotation/07_blastp_nr.sh` | `1blast_nr.sh` | BLASTP (blastp-fast) of CRISPR-contig ORFs against NCBI NR with taxonomy (`-taxids 2`, `staxids`/`stitle` in the output). |
| `pipelines/03_functional_annotation/08_diamond_nr_prophage.sh` | `2orf_diamond_li.sh` | DIAMOND blastp of the high-quality prophage ORFs against the NR diamond database (>=60% identity, >=90% subject coverage). |
| `pipelines/03_functional_annotation/09_diamond_taxonomy.sh` | `diamond_res_tax.sh` | Turn the DIAMOND best hits into standardised taxonomy with the NCBI prot.accession2taxid and taxdump files. |
| `pipelines/03_functional_annotation/10_hhsuite_plasmid.sh` | `1hhsuite2.sh` | HH-suite (hhblits/hhsearch) screening of the non-redundant plasmid proteins against the HH-suite databases. |
| `pipelines/03_functional_annotation/11_clans_djr_mcp.sh` | `1CLANS.sh` | CLANS all-against-all clustering of the DJR major capsid proteins, with the post-processing of the sequence names. |
| `pipelines/04_crispr_and_host_links/12_crisprcasfinder.sh` | `1CRISPRCasFinder.sh` | CRISPRCasFinder (levelMin 4) on the minced CRISPR contigs. |
| `pipelines/04_crispr_and_host_links/13_crisprcaster.sh` | `1CRISPRTyper.sh` | CRISPRCasTyper (cctyper) typing of the level-4 CRISPR arrays. |
| `pipelines/04_crispr_and_host_links/14_makeblastdb_imgvr.sh` | `2makeblastdb.sh` | Build the IMG/VR BLAST database (deduplication of repeated sequence IDs first). |
| `pipelines/04_crispr_and_host_links/15_blastn_spacers_imgvr.sh` | `3blastn_IMGVR.sh` | BLASTN of the CRISPR spacers against IMG/VR; keeps hits with <=1 mismatch and >=95% coverage. |
| `pipelines/04_crispr_and_host_links/16_filter_spacer_hits.sh` | `4res_filter.sh` | Filter the spacer-virus and spacer-plasmid hits to the host contigs of interest and extract the targeted sequences. |
| `pipelines/05_defense_systems/17_defensefinder.sh` | `01DefenseFinder.sh` | DefenseFinder (`--preserve-raw`) on the 4,696 predator genomes; the resulting systems table feeds the host-defense analyses. |
| `pipelines/06_utilities/gff_to_genbank.py` | `1mobilome_predatory_bacteria/5special_plasmid/gff_to_genbank.py` | Convert a GFF3 plus its FASTA into GenBank format for deposition of the special plasmids/phages. |

## R analyses

| New path | Previous working-copy path | Manuscript | Title |
| --- | --- | --- | --- |
| `analyses/01_genome_catalog/01_genome_metadata.R` | `01_genomes/1_genome_metadata.R` | Fig. S1 | Myxobacterial genome catalog |
| `analyses/01_genome_catalog/02_genome_tree_annotation.R` | `01_genomes/1_MAGs_itol_anno.R` | Fig. 2A | 565-genome phylogeny annotation |
| `analyses/01_genome_catalog/03_CRISPR_Cas_taxonomy.R` | `1mobilome_predatory_bacteria/2CRISPRTyper/1_CRISPRCas_tax.R` | Fig. 2A | CRISPR-Cas type and subtype assignment |
| `analyses/02_phage_catalog/01_viral_genome_info.R` | `01_phage_genome_info/1_virus_genome_info.R` | Fig. 3A-F | Myxobacterial viral catalog |
| `analyses/02_phage_catalog/02_biosample_info_extract.R` | `01_phage_genome_info/2_biosample_info_extract.R` | Fig. 3D | Biosample metadata extraction |
| `analyses/02_phage_catalog/03_phage_metadata.R` | `1mobilome_predatory_bacteria/3vcontact2/1_phage_metadata.R` | Fig. 3C | Phage metadata and host assignment |
| `analyses/02_phage_catalog/04_phage_geographical_distribution.R` | `1mobilome_predatory_bacteria/3vcontact2/2_phage_geographical_distribution.R` | Fig. 3D | Geographical distribution of the phages |
| `analyses/02_phage_catalog/05_AMG_statistics.R` | `06_allphage_AMGs/1_AMG_res_statistatistics.R` | Fig. 3G | Auxiliary metabolic genes |
| `analyses/02_phage_catalog/06_plasmid_metadata.R` | `1mobilome_predatory_bacteria/3vcontact2_plasmid/1_plasmid_metadata.R` | Plasmid catalog (Fig. 1I) | Plasmid metadata and host assignment |
| `analyses/02_phage_catalog/07_plasmid_geographical_distribution.R` | `1mobilome_predatory_bacteria/3vcontact2_plasmid/2_plasmid_geographical_distribution.R` | Plasmid catalog | Geographical distribution of the plasmids |
| `analyses/03_phage_annotation/01_phage_gene_annotation.R` | `03_phage_gene_anno/1_phage_gene_annotation.R` | Fig. 5 | Gene maps of the proposed myxobacterial viruses |
| `analyses/03_phage_annotation/02_candidate_phage_gene_annotation.R` | `1mobilome_predatory_bacteria/7phage_gene_annotation/1_phage_genome_annotation.R` | Fig. 5 | Gene maps of the 17 candidate phages |
| `analyses/03_phage_annotation/03_MGE_genome_annotation.R` | `08_special_MGE/1_MGE_genome_annotation.R` | Fig. S2 | Gene map of a representative MGE |
| `analyses/03_phage_annotation/04_phage_genome_features.R` | `07_phage_molecular_chara/2_phage_genome_features.R` | Fig. S3 | Genomic and proteomic characteristics of the viruses |
| `analyses/03_phage_annotation/05_gene_presence_dotplot.R` | `1mobilome_predatory_bacteria/7phage_gene_annotation/phage_gene_anno_plot/create_virus_dotplot.R` | Fig. S6 | Gene presence/absence dot-plot |
| `analyses/03_phage_annotation/06_gene_dotplot_simple.R` | `1mobilome_predatory_bacteria/7phage_gene_annotation/phage_gene_anno_plot/virus_dotplot_simple.R` | Fig. S6 (variant) | Gene dot-plot, simple variant |
| `analyses/03_phage_annotation/07_gene_dotplot_full.R` | `1mobilome_predatory_bacteria/7phage_gene_annotation/phage_gene_anno_plot/virus_protein_dotplot.R` | Fig. S6 (variant) | Gene dot-plot, full variant |
| `analyses/04_new_viral_groups/01_aai_heatmap.R` | `04_new_phage/1_aai_heatmap.R` | Fig. 4D | AAI heatmap of the new viral groups |
| `analyses/04_new_viral_groups/02_candidate_17_phage_metadata.R` | `1mobilome_predatory_bacteria/6phage_phylo/1_all_17_phage_metadata.R` | Fig. 4 | Metadata of the 17 candidate phages |
| `analyses/05_marker_phylogeny/01_marker_gene_annotation.R` | `09_special_phage_marker_phylogeny/1_phage_gene_annotation.R` | Fig. S4 | Marker-gene phylogeny annotation |
| `analyses/05_marker_phylogeny/02_Cas12m_tree_annotation.R` | `1mobilome_predatory_bacteria/9special_protein/cas12m/cas12_phylo/1_Cas12m_itol_anno.R` | Fig. 6D | Type V effector (Cas12m) phylogeny |
| `analyses/05_marker_phylogeny/03_DJR_MCP_sequence_names.R` | `1mobilome_predatory_bacteria/9special_protein/DJR_MCP/1_seqname_treatment.R` | Fig. S5 | DJR major-capsid protein identifiers |
| `analyses/06_host_function/01_myxobacteria_functional_annotation.R` | `10_myxo_function_anno/1_myxo_functional_anno.R` | Fig. 1 | Functional repertoire of the myxobacterial genomes |
| `analyses/06_host_function/02_host_defense_systems.R` | `1mobilome_predatory_bacteria/10host_DS/1_Host_Defense.R` | Fig. 1A | Defense systems of the phage hosts |

## Misc

| New path | Previous path | Note |
| --- | --- | --- |
| `misc/01_vcontact3_analysis.R` | `02_phage_vcontact3/1_vcontact3_analysis.R` | vContact3 phage/provirus analysis of a phage set; the working directory still points at the ESKAPE-virome project (GTR_WORK/2ESKAPE_MGE), so the script is kept for reference only. |

## Files that were deliberately left out

| File | Why |
| --- | --- |
| `1mobilome_predatory_bacteria/1genome_metadata/1_genome_metadata.R` | 474-line draft of `01_genomes/1_genome_metadata.R`; the 541-line version is kept |
| `1mobilome_predatory_bacteria/4genome_phylo/1_MAGs_itol_anno.R` | 305-line draft of `01_genomes/1_MAGs_itol_anno.R`; the 315-line version is kept |
| `1mobilome_predatory_bacteria/6phage_phylo/comparem_aai/1_aai_heatmap.R` | 71-line draft of `04_new_phage/1_aai_heatmap.R`; the 196-line version is kept |
| `try`, `R code/try` | empty placeholder files |

## Site-specific paths and how they are hidden

The absolute paths of the working copies -- user home directories, a local Windows
analysis tree, an earlier working copy and a colleague's defensome folder -- have been
removed from every file and replaced by variables:

| Variable | Replaces | Used in |
| --- | --- | --- |
| `GTR_DATA` | the current analysis root (`06_predatory_bacteria_virome`) | every script in `analyses/` |
| `GTR_WORK` | the earlier working root that holds `1mobilome_predatory_bacteria/` and the sibling `2ESKAPE_MGE/` project | `02_phage_catalog/03-07`, `03_phage_annotation/02`, `04_new_viral_groups/02`, `05_marker_phylogeny/02-03`, `06_host_function/02`, `misc/01` |
| `GTR_DEF` | the shared defensome data folder | `03_phage_annotation/07` |
| `LUSTRE_HOME` | home directory on the large-storage filesystem | all shell scripts |
| `DB_ROOT` | reference databases (`${LUSTRE_HOME}/database` by default) | the scripts that call geNomad, CheckV, BLAST, HH-suite, DIAMOND or NCBI taxonomy |

Each R script declares only the variables it uses:

```r
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
GTR_WORK <- normalizePath(Sys.getenv("GTR_WORK", unset = "."), mustWork = FALSE)
```

`normalizePath()` makes the paths absolute at load time, so they stay valid after the
`setwd()` calls inside the scripts.

