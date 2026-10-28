# Data availability and input files

The repository contains **code only**. Every table read by the R scripts is listed
below, so the analyses can be re-run once the data are in place.

Three roots are used:

* `GTR_DATA` — the current analysis root (`06_predatory_bacteria_virome`); default: the
  current working directory, override with
  `Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")`.
* `GTR_WORK` — the earlier working root holding `1mobilome_predatory_bacteria/` (and the
  sibling `2ESKAPE_MGE/` project used by `misc/01_vcontact3_analysis.R`).
* `GTR_DEF` — the shared defensome data folder used by the dot-plot variant.

The sub-directory names are the ones used during the study:

```
<GTR_DATA>/
├── 01_genomes/                  selected MAGs list
├── 01_phage_genome_info/        viral summary, iphop, CheckV, lifecycle tables
├── 02_phage_vcontact3/          vContact3 output (misc only)
├── 03_phage_gene_anno/          GFF of the selected phages
├── 04_new_phage/                CompareM AAI summary
├── 06_allphage_AMGs/            AMG annotation tables
├── 07_phage_molecular_chara/    per-phage molecular features
├── 08_special_MGE/              GFF of the representative MGE
├── 09_special_phage_marker_phylogeny/   ViPTree GFF + tree
├── 10_myxo_function_anno/       DS / TXSS / AMP / BGC annotation of the genomes
└── 1mobilome_predatory_bacteria/
    ├── 1genome_metadata/        genome records of the public sources
    ├── 2CRISPRTyper/            CRISPRCasTyper output
    ├── 3vcontact2/              vContact2 phage metadata
    ├── 3vcontact2_plasmid/      vContact2 plasmid metadata
    ├── 4genome_phylo/           tree + iTOL input
    ├── 6phage_phylo/            candidate phage lists
    ├── 7phage_gene_annotation/  GFF of the candidate phages
    ├── 9special_protein/        Cas12m and DJR MCP data
    └── 10host_DS/              DefenseFinder output of the 4,696 genomes

<GTR_WORK>/
├── 1mobilome_predatory_bacteria/   the same tree as above, earlier state
└── 2ESKAPE_MGE/5phage_vcontact3/   ESKAPE-virome vContact3 output (misc/01 only)
```

## Where the data come from

| Data | Source |
| --- | --- |
| Myxobacterial genomes | NCBI GenBank + the CRBC/GEM/GOMC/GROWdb/SMAG collections, plus the MAGs recovered here |
| Raw metagenomic reads of the new MAGs | NCBI SRA, accession given in the paper |
| Virus and plasmid calls | regenerate with `pipelines/01-02` |
| Functional annotation (NR, taxonomy, HH-suite) | regenerate with `pipelines/03` |
| CRISPR arrays and spacer-host links | regenerate with `pipelines/04` |
| Host defense systems | regenerate with `pipelines/05` |
| Reference viral genomes (RefSeq, IMG/VR, VIRE) | public databases, cited in the paper |

## Per-script inputs


### `analyses/01_genome_catalog/01_genome_metadata.R`

**Myxobacterial genome catalog** — Merge the genome records of all sources (GenBank, CRBC, GEM, GOMC, GROWdb, SMAG) with the newly recovered MAGs, then build the environmental-origin, genome-size/GC and taxonomic-composition panels.

* Manuscript: Fig. S1
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/1genome_metadata`, `01_genomes`
* Reads:
    * `gtdb_r226_predators_metadata_new_2.tsv`
    * `CRBC_merged_results.tsv`
    * `GEM_merged_results.tsv`
    * `GOMC_merged_results.tsv`
    * `GROWdb_merged_results.tsv`
    * `SMAG_merged_results.tsv`
    * `TG2G_merged_results.tsv`
    * `myxococcota_mag_metadata.tsv`
* Writes: `1_MAG_tax_sankey.pdf`, `1_genome_metadata.pdf`, `1_genome_metadata.RData`

### `analyses/01_genome_catalog/02_genome_tree_annotation.R`

**565-genome phylogeny annotation** — iTOL annotation of the 565-genome maximum-likelihood tree: class-level taxonomy, method used for virus-host assignment and the CRISPR-Cas type/subtype rings.

* Manuscript: Fig. 2A
* Working directory: `01_genomes`
* Reads:
    * `Selected_phylo_MAGs_filtered_Myxoco.txt`
    * `myxococcota_mag_metadata.tsv`
    * `CRISPR_Cas.tab`
    * `spacer_target_virus_MAGs.txt`
    * `Contain_provirus_MAGs.txt`
    * `./1_MAGs_itol_anno.RData`
* Writes: `1_MAGs_itol_anno.RData`, `iTOL annotation files`

### `analyses/01_genome_catalog/03_CRISPR_Cas_taxonomy.R`

**CRISPR-Cas type and subtype assignment** — Assign every detected CRISPR-Cas system to its type/subtype from the CRISPRCasTyper output and tabulate it per genome.

* Manuscript: Fig. 2A
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/2CRISPRTyper`
* Reads:
    * `CRISPR_Cas.tab`
    * `all_4649_genome_metadata.tsv`
* Writes: `per-genome CRISPR-Cas type table`

### `analyses/02_phage_catalog/01_viral_genome_info.R`

**Myxobacterial viral catalog** — Core viral catalog: retain the viral sequences assigned to Myxococcota hosts, then build the CheckV quality, vOTU taxonomy, host-class, ecosystem/geography and lifestyle panels.

* Manuscript: Fig. 3A-F
* Working directory: `01_phage_genome_info`
* Reads:
    * `all_4649_genome_metadata.tsv`
    * `spacers_IMGVR.blastn.1mis100cov.out2`
    * `spacers_VIRE.blastn.1mis100cov.out2`
    * `genome_derived_virus.txt`
    * `filtered_validated_proviruses.tsv`
    * `phage_detailed_info_trim.tsv`
    * `iphop_phage2hostgenome.tsv`
    * `myxo_IMGVR_phage.tsv`
    * `IMG.xlsx`
    * `myxo_vire_phage.tsv`
    * `vire_matched_rows.tsv`
    * `./1_virus_genome_info.RData`
* Writes: `1_virus_overview_info.pdf`, `1_p_phage_quality2.pdf`, `1_p_phage_length_density.pdf`, `1_p_phage_sankey_4level.pdf`, `1_phage_geographical_distribution.pdf`, `1_p_major_phage.pdf`, `1_virus_genome_info.RData`

### `analyses/02_phage_catalog/02_biosample_info_extract.R`

**Biosample metadata extraction** — Parse the NCBI Biosample records of the reference phages (isolation source, country, coordinates) into a tidy metadata table used for the geography panel.

* Manuscript: Fig. 3D
* Working directory: `01_phage_genome_info`
* Reads: -
* Writes: `biosample metadata table`

### `analyses/02_phage_catalog/03_phage_metadata.R`

**Phage metadata and host assignment** — Assemble the per-phage metadata from the vContact2 run, including the predicted host genome of each phage.

* Manuscript: Fig. 3C
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/3vcontact2`
* Reads:
    * `self_prepared_phage.csv`
    * `nr_phage_rep_seq_virus_summary.tsv`
    * `filtered_spacers_virus.blastn.1mis95cov.out`
    * `spacers_virus.blastn.1mis95cov.out`
    * `filtered_validated_proviruses.tsv`
    * `all_4649_genome_metadata.tsv`
* Writes: `self_selected_phage_node.tsv and related tables`

### `analyses/02_phage_catalog/04_phage_geographical_distribution.R`

**Geographical distribution of the phages** — Link each phage to the isolation site of its host and summarise the global distribution by country and ecosystem.

* Manuscript: Fig. 3D
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/3vcontact2`
* Reads:
    * `self_selected_phage_node.tsv`
    * `filtered_IMGVR_sequences.tsv`
    * `all_4649_genome_metadata.tsv`
    * `IMG_genomeCart_metadata.xlsx`
    * `gtdb_r226_predators_metadata_new_2.tsv`
    * `GEM_genome_metadata.tsv`
    * `self_selected_phage_3.xlsx`
* Writes: `2_phage_geographical_distribution.pdf`

### `analyses/02_phage_catalog/05_AMG_statistics.R`

**Auxiliary metabolic genes** — Count the auxiliary metabolic genes per phage, group them by KEGG category and cross-tabulate them with the virulent/temperate lifestyle.

* Manuscript: Fig. 3G
* Working directory: `06_allphage_AMGs`
* Reads:
    * `output_annotated_all_AMG.ORF2GENE_with_meta.tsv`
    * `Myxococcota_all_phage_info_filterd.fa.bacphlip`
    * `1_AMG_res_statistatistics.RData`
* Writes: `1_AMG_res_statistatistics.RData`, `AMG summary tables`

### `analyses/02_phage_catalog/06_plasmid_metadata.R`

**Plasmid metadata and host assignment** — Assemble the per-plasmid metadata (source, size, predicted host) of the non-redundant plasmid set.

* Manuscript: Plasmid catalog (Fig. 1I)
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/3vcontact2_plasmid`
* Reads:
    * `seq_info_table.tsv`
    * `nr_plasmid_rep_seq_plasmid_summary.tsv`
    * `filtered_spacers_plasmid.blastn.1mis95cov.out`
    * `spacers_plasmid.blastn.1mis95cov.out`
    * `all_4649_genome_metadata.tsv`
* Writes: `self_selected_plasmid.tsv and related tables`

### `analyses/02_phage_catalog/07_plasmid_geographical_distribution.R`

**Geographical distribution of the plasmids** — Map the plasmids onto the isolation sites of their hosts and summarise the distribution.

* Manuscript: Plasmid catalog
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/3vcontact2_plasmid`
* Reads:
    * `self_selected_plasmid.tsv`
    * `filtered_IMGPR_sequences.tsv`
    * `all_4649_genome_metadata.tsv`
    * `IMG_genomeCart_metadata.xlsx`
    * `gtdb_r226_predators_metadata_new_2.tsv`
    * `GEM_genome_metadata.tsv`
    * `self_selected_plasmid_3.xlsx`
* Writes: `2_plasmid_geographical_distribution3.pdf`

### `analyses/03_phage_annotation/01_phage_gene_annotation.R`

**Gene maps of the proposed myxobacterial viruses** — Gene-arrow maps of the selected phages with their functional categories, plus the protospacer positions.

* Manuscript: Fig. 5
* Working directory: `03_phage_gene_anno`, `GTR_WORK/1mobilome_predatory_bacteria/7phage_gene_annotation`
* Reads:
    * `selected_phage_treated.gff`
    * `phage_annotation_summary.tsv`
    * `all_phage_gff_anno.tsv`
    * `Candidate_17_new_phages.tsv`
    * `nr_phage_rep_seq_virus_summary.tsv`
    * `filtered_candidate17.out`
    * `all_4649_genome_metadata.tsv`
    * `selected_phage_anno_withtype.tsv`
    * `novel_6_complete_phage_info.xlsx`
    * `./1_phage_gene_annotation.RData`
* Writes: `p_gene_dot.pdf`, `1_phage_gene_annotation.RData`

### `analyses/03_phage_annotation/02_candidate_phage_gene_annotation.R`

**Gene maps of the 17 candidate phages** — Gene-arrow maps of the 17 newly identified candidate phages, one panel per phage.

* Manuscript: Fig. 5
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/7phage_gene_annotation`
* Reads:
    * `Candidate_17_new_phages_treated.gff`
    * `phage_annotation_summary.tsv`
    * `all_17_phage_gff_anno.tsv`
    * `Candidate_17_new_phages.tsv`
    * `nr_phage_rep_seq_virus_summary.tsv`
    * `filtered_candidate17.out`
    * `all_4649_genome_metadata.tsv`
* Writes: `17-phage gene map figures`

### `analyses/03_phage_annotation/03_MGE_genome_annotation.R`

**Gene map of a representative MGE** — Gene-arrow annotation of the representative mobile genetic element carrying cargo defense/secretion genes.

* Manuscript: Fig. S2
* Working directory: `08_special_MGE`
* Reads:
    * `selected_MGE.gff`
* Writes: `MGE gene annotation figure`

### `analyses/03_phage_annotation/04_phage_genome_features.R`

**Genomic and proteomic characteristics of the viruses** — Genome size, GC content, amino-acid usage, average protein length and predicted pI of the viruses grouped by host class, and the GC comparison between virus and host.

* Manuscript: Fig. S3
* Working directory: `07_phage_molecular_chara`
* Reads:
    * `myxoco_all_phage_info.tsv`
    * `AAU_matrix.csv`
    * `pI_distribution_matrix.csv`
    * `Molecular_Final_Full_Stats.txt`
* Writes: `2_phage_genome_features.RData`, `virus feature figures`

### `analyses/03_phage_annotation/05_gene_presence_dotplot.R`

**Gene presence/absence dot-plot** — Dot-plot of the annotated genes across a set of viruses (filled circle = orthologue present); the script installs its packages and builds the plot from a prepared table.

* Manuscript: Fig. S6
* Working directory: -
* Reads:
    * `all_17_phage_gff_withtype2_simplifiedname.xlsx`
    * `Candidate_14_new_phages.xlsx`
* Writes: `p_phage_gene_anno2.pdf`

### `analyses/03_phage_annotation/06_gene_dotplot_simple.R`

**Gene dot-plot, simple variant** — Simplified variant of the gene dot-plot, driven by embedded example data; useful as a template for panel styling.

* Manuscript: Fig. S6 (variant)
* Working directory: `GTR_DEF/代码调试`
* Reads: -
* Writes: `virus_protein_dotplot.pdf`

### `analyses/03_phage_annotation/07_gene_dotplot_full.R`

**Gene dot-plot, full variant** — Full variant of the gene dot-plot with the complete styling; reads and writes inside an external prototyping folder (GTR_DEF).

* Manuscript: Fig. S6 (variant)
* Working directory: `GTR_DEF/代码调试`
* Reads: -
* Writes: `virus_protein_dotplot.pdf`

### `analyses/04_new_viral_groups/01_aai_heatmap.R`

**AAI heatmap of the new viral groups** — Build the average-amino-acid-identity matrix from the CompareM output and draw the clustered heatmaps of the phage set.

* Manuscript: Fig. 4D
* Working directory: `04_new_phage`
* Reads:
    * `aai_summary.tsv`
    * `VIRIDIC_sim-dist_table.tsv`
    * `./1_aai_heatmap.RData`
* Writes: `1_ht1.pdf`, `1_ht2.pdf`, `1_aai_heatmap.RData`

### `analyses/04_new_viral_groups/02_candidate_17_phage_metadata.R`

**Metadata of the 17 candidate phages** — Combine the 17 candidate phages with the geNomad virus summary to obtain their genome size, topology and completeness.

* Manuscript: Fig. 4
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/6phage_phylo`
* Reads:
    * `Candidate_17_new_phages.tsv`
    * `nr_phage_rep_seq_virus_summary.tsv`
* Writes: `17-phage metadata table`

### `analyses/05_marker_phylogeny/01_marker_gene_annotation.R`

**Marker-gene phylogeny annotation** — Annotate the 39-phage ViPTree maximum-likelihood tree: gene maps aligned to the tips plus the marker-gene presence strips.

* Manuscript: Fig. S4
* Working directory: `09_special_phage_marker_phylogeny`
* Reads:
    * `viptree_selected_39_phage_clean_treated.gff`
    * `39phage_annotation_summary.tsv`
    * `selected_phage_anno_withtype.tsv`
    * `ESKAPE_complete_phage_treated.gff`
    * `phage_annotation_summary.tsv`
    * `ESKAPE_all_phage_info_filter.tsv`
* Writes: `1_phage_gene_annotation.RData`, `tree annotation figures`

### `analyses/05_marker_phylogeny/02_Cas12m_tree_annotation.R`

**Type V effector (Cas12m) phylogeny** — iTOL annotation of the collected Cas12m proteins and TnpB nucleases (taxonomy and protein features).

* Manuscript: Fig. 6D
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/9special_protein/cas12m/cas12_phylo`
* Reads:
    * `collected_cas12.xlsx`
* Writes: `iTOL annotation files`

### `analyses/05_marker_phylogeny/03_DJR_MCP_sequence_names.R`

**DJR major-capsid protein identifiers** — Build the short-to-full identifier table of the DJR major capsid proteins used as the CLANS input.

* Manuscript: Fig. S5
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/9special_protein/DJR_MCP`
* Reads:
    * `merged_DJRMCP_seqname2.txt`
    * `clans_net`
    * `41564_2024_1839_MOESM13_ESM.xlsx`
* Writes: `merged_DJRMCP_seqname table`

### `analyses/06_host_function/01_myxobacteria_functional_annotation.R`

**Functional repertoire of the myxobacterial genomes** — The largest analysis of the study: defense systems, secretion systems (TXSS), antimicrobial peptides (AMP) and biosynthetic gene clusters (BGC) across classes, their prevalence per class and the functional cargo of the MGEs.

* Manuscript: Fig. 1
* Working directory: `10_myxo_function_anno`
* Reads:
    * `myxococcota_mag_metadata.tsv`
    * `merged_myxo_mag_defense_finder_genes.tsv`
    * `DF_genes.gff`
    * `merged_myxo_mag_defense_finder_systems.tsv`
    * `macrel.out.prediction`
    * `AMP_genes.gff`
    * `merged_best_solutions.tsv`
    * `TXSS_genes.gff`
    * `all_merged_BGC.tsv`
    * `BGC_genes.gff`
    * `merged_myxo_mag_virus_summary.tsv`
    * `merged_myxo_mag_plasmid_summary.tsv`
    * `3229MAGs_filtered_completeIS.tsv`
    * `merged_myxo_mag.integrons`
    * `processed_ice_sequences_ID.txt`
    * `all_mge_gff.tsv`
    * `./1_myxo_functional_anno.RData`
* Writes: `1_myxo_functional_anno.RData`, `DS/TXSS/AMP/BGC and MGE-cargo figures`

### `analyses/06_host_function/02_host_defense_systems.R`

**Defense systems of the phage hosts** — Summarise the DefenseFinder systems of the 4,696 predator genomes and of the phage hosts, and plot the DS distribution.

* Manuscript: Fig. 1A
* Working directory: `GTR_WORK/1mobilome_predatory_bacteria/10host_DS`
* Reads:
    * `all_4696_predators_genome_renamed_defense_finder_systems.tsv`
    * `phage_protospacer_anno.tsv`
* Writes: `1_p_Host_DS.pdf`

### `misc/01_vcontact3_analysis.R`

**vContact3 analysis (inherited from the ESKAPE-virome project)** — vContact3 phage/provirus analysis of a phage set; the working directory still points at the ESKAPE-virome project (GTR_WORK/2ESKAPE_MGE), so the script is kept for reference only.

* Manuscript: -
* Working directory: `GTR_WORK/2ESKAPE_MGE/5phage_vcontact3`
* Reads:
    * `ESKAPE_all_phage_info_filter.tsv`
    * `part1.cyjs default node.csv`
    * `1_vcontact3_analysis.RData`
* Writes: `1_vcontact3_analysis.RData`

