# Directory layout

```
Predatory-bacterial-virome/
├── README.md                  English overview
├── README_zh.md               中文说明
├── docs/                      workflow, layout, mapping, data, software
├── pipelines/                 layer 1: bioinformatics, Slurm + conda
├── analyses/                  layer 2: statistics and figures, R
└── misc/                      scripts kept for reference, outside the manuscript pipeline
```

Both layers use the same convention: a two-digit directory prefix gives the execution /
narrative order, and a two-digit file prefix gives the order inside the directory.

## `pipelines/` — upstream, ordered by stage

| Directory | Scope |
| --- | --- |
| `01_genome_and_phage_prediction/` | Prodigal on the 4,696 genomes, geNomad, VirSorter2 |
| `02_viral_quality_and_clustering/` | CheckV, mmseqs2 clustering, CompareM AAI |
| `03_functional_annotation/` | BLASTP NR, DIAMOND, NCBI taxonomy, HH-suite, CLANS |
| `04_crispr_and_host_links/` | CRISPRCasFinder, CRISPRCasTyper, IMG/VR BLAST and hit filtering |
| `05_defense_systems/` | DefenseFinder |
| `06_utilities/` | GFF -> GenBank conversion |

The scripts are the Slurm job scripts as they were submitted, with the `#SBATCH` headers
preserved. Site-specific paths are variables declared in the `# Site-specific settings`
block below each `#SBATCH` header:

```bash
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
DB_ROOT="${DB_ROOT:-${LUSTRE_HOME}/database}"       # reference databases
```

File numbers are re-sequenced (`01` … `17`) so that they follow the pipeline order; the
original names are listed in [03_file_mapping.md](03_file_mapping.md).

## `analyses/` — downstream, ordered by the manuscript

| Directory | Manuscript | Content |
| --- | --- | --- |
| `01_genome_catalog/` | Fig. S1, Fig. 2A | genome collection, 565-genome tree, CRISPR-Cas types |
| `02_phage_catalog/` | Fig. 3 | viral catalogue, host assignment, geography, AMGs, plasmids |
| `03_phage_annotation/` | Fig. 5; Figs. S2, S3, S6 | gene maps, MGE map, molecular features, dot-plots |
| `04_new_viral_groups/` | Fig. 4 | AAI heatmaps, candidate phage metadata |
| `05_marker_phylogeny/` | Fig. S4, S5; Fig. 6D | marker-gene tree, Cas12m tree, CLANS names |
| `06_host_function/` | Fig. 1 | DS/TXSS/AMP/BGC repertoire, host defense systems |

## Script header

Every R script starts with the same block:

```r
# Predatory-bacterial-virome -- phage catalog / Myxobacterial viral catalog
# Core viral catalog: retain the viral sequences assigned to Myxococcota hosts, ...
#
# Manuscript      : Fig. 3A-F
# Main outputs    : 1_virus_overview_info.pdf, ...
# Original script : 01_phage_genome_info/1_virus_genome_info.R
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
```

Three roots are available and declared only when used:

| Variable | Meaning |
| --- | --- |
| `GTR_DATA` | the current analysis root (`06_predatory_bacteria_virome`) |
| `GTR_WORK` | the earlier working root (`1mobilome_predatory_bacteria/`, sibling `2ESKAPE_MGE/`) |
| `GTR_DEF` | the shared defensome data folder used by the dot-plot variant |

## `misc/`

`misc/01_vcontact3_analysis.R` is a vContact3 analysis whose working directory still
points at the sibling ESKAPE-virome project. It is kept because it documents how the
vContact3 output was handled, but it is not part of the manuscript pipeline.

## `docs/`

| File | Content |
| --- | --- |
| `01_pipeline_overview.md` | end-to-end workflow, stage by stage |
| `02_directory_layout.md` | this file |
| `03_file_mapping.md` | previous path -> new path for every script |
| `04_data_availability.md` | working directory and input tables of every R script |
| `05_software_and_environment.md` | conda environments, R packages, cluster notes |
