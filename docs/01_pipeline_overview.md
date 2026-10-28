# Workflow overview

The study has two layers. The **upstream layer** (`pipelines/`) turns the collected
myxobacterial genomes into a viral catalogue and the annotation tables with Slurm job
scripts; the **downstream layer** (`analyses/`) turns those tables into the statistics and
figures of the paper with R. The layers communicate only through files.

```
   myxobacterial genomes (GenBank, CRBC, GEM, GOMC, GROWdb, SMAG + new MAGs)
        |
   [1] genome and phage prediction -----------------------------------------+
       Prodigal ORFs, geNomad, VirSorter2                                   |
        |                                                                    |
   [2] viral quality and clustering --+                                     |
       CheckV, mmseqs2, CompareM AAI   |   [3] functional annotation -----+  |
        |                              |       BLASTP NR, DIAMOND NR      |  |
        |                              |       taxonomy, HH-suite, CLANS  |  |
        |                              |                                  |  |
   [4] CRISPR and host links ----------+   [5] defense systems ----------+  |
       CRISPRCasFinder, CRISPRCasTyper      DefenseFinder                  |
       IMG/VR BLAST + hit filtering                                        |
        |                              |                                  |  |
        v                              v                                  v  v
   [6] R analyses: analyses/01_genome_catalog ... analyses/06_host_function
       statistics, figures, supplementary tables
```

## Stage 1 — genome and phage prediction

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 1.1 | `01_prodigal_genomes.sh` | Prodigal (meta) | ORF GFF/FNA/FAA of the 4,696 predator genomes |
| 1.2 | `02_genomad_genomes.sh` | geNomad | plasmid and virus/provirus calls |
| 1.3 | `03_virsorter2.sh` | VirSorter2 | second viral-sequence caller (split genome files) |

## Stage 2 — viral quality and clustering

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 2.1 | `04_checkv.sh` | CheckV (db 1.5) | completeness/contamination, quality tiers |
| 2.2 | `05_mmseqs_cluster_phages.sh` | mmseqs2 | non-redundant phage set (95% ANI, 90% coverage) |
| 2.3 | `06_comparem_aai.sh` | CompareM | pairwise AAI of the 48 phages |

## Stage 3 — functional annotation

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 3.1 | `07_blastp_nr.sh` | BLAST+ 2.16 | ORF hits against NCBI NR with taxonomy |
| 3.2 | `08_diamond_nr_prophage.sh` | DIAMOND | prophage ORF hits against NR |
| 3.3 | `09_diamond_taxonomy.sh` | NCBI taxdump | standardised taxonomic annotation of the DIAMOND hits |
| 3.4 | `10_hhsuite_plasmid.sh` | HH-suite | structural/remote-homology screening of the plasmid proteins |
| 3.5 | `11_clans_djr_mcp.sh` | CLANS | all-against-all clustering of the DJR capsid proteins |

## Stage 4 — CRISPR arrays and historical host links

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 4.1 | `12_crisprcasfinder.sh` | CRISPRCasFinder | CRISPR arrays and Cas loci |
| 4.2 | `13_crisprcaster.sh` | CRISPRCasTyper | type/subtype assignment of the arrays |
| 4.3 | `14_makeblastdb_imgvr.sh` | BLAST+ | IMG/VR nucleotide database |
| 4.4 | `15_blastn_spacers_imgvr.sh` | BLAST+ | spacer -> IMG/VR hits (<=1 mismatch, >=95% coverage) |
| 4.5 | `16_filter_spacer_hits.sh` | seqkit/grep | hits restricted to the host contigs of interest |

## Stage 5 — defense systems

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 5.1 | `17_defensefinder.sh` | DefenseFinder (`--preserve-raw`) | defense systems of the 4,696 genomes |

## Stage 6 — utility

| Script | Purpose |
| --- | --- |
| `gff_to_genbank.py` | convert a GFF3 + FASTA into GenBank for deposition |

## R analyses

| Module | Manuscript | Content |
| --- | --- | --- |
| `01_genome_catalog` | Fig. S1, Fig. 2A | genome metadata, 565-genome tree annotation, CRISPR-Cas type assignment |
| `02_phage_catalog` | Fig. 3 | viral catalogue (quality, taxonomy, hosts, geography), AMGs, plasmid catalog |
| `03_phage_annotation` | Fig. 5; Figs. S2, S3, S6 | gene maps, MGE map, genomic/proteomic features, gene presence dot-plots |
| `04_new_viral_groups` | Fig. 4 | AAI heatmaps of the new viral groups, candidate phage metadata |
| `05_marker_phylogeny` | Fig. S4, S5; Fig. 6D | marker-gene tree annotation, Cas12m tree, CLANS input names |
| `06_host_function` | Fig. 1 | DS/TXSS/AMP/BGC repertoire of the genomes and the MGE cargo |
| `misc` | - | one vContact3 script inherited from the sibling ESKAPE project |

Per-script inputs are listed in
[04_data_availability.md](04_data_availability.md); the script -> figure -> original path
mapping is in [03_file_mapping.md](03_file_mapping.md).
