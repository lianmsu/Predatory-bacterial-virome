# Predatory-bacterial-virome

Code accompanying **"The virome of globally distributed myxobacteria"**
(Li, P. & Ni, J.).

The repository contains every script of the study: the upstream bioinformatic pipeline
that goes from the collected myxobacterial genomes to the viral catalogue and its
annotation, and the R analysis scripts that produce the statistics and the figures.

> 中文版说明见 **[README_zh.md](README_zh.md)**。

---

## What is in here

```
Predatory-bacterial-virome/
├── pipelines/          # upstream bioinformatics (Slurm + conda, shell/python)
│   ├── 01_genome_and_phage_prediction/    Prodigal, geNomad, VirSorter2
│   ├── 02_viral_quality_and_clustering/   CheckV, mmseqs2, CompareM AAI
│   ├── 03_functional_annotation/          BLASTP NR, DIAMOND, NCBI taxonomy, HH-suite, CLANS
│   ├── 04_crispr_and_host_links/          CRISPRCasFinder, CRISPRCasTyper, IMG/VR BLAST
│   ├── 05_defense_systems/                DefenseFinder
│   └── 06_utilities/                      GFF -> GenBank
├── analyses/           # downstream statistics and figures (R)
│   ├── 01_genome_catalog/                 Fig. S1, Fig. 2A
│   ├── 02_phage_catalog/                  Fig. 3
│   ├── 03_phage_annotation/               Fig. 5; Figs. S2, S3, S6
│   ├── 04_new_viral_groups/               Fig. 4
│   ├── 05_marker_phylogeny/               Fig. S4, S5; Fig. 6D
│   └── 06_host_function/                  Fig. 1
├── misc/               # one script kept for reference, outside the manuscript pipeline
└── docs/               # workflow, layout, file mapping, data and software notes
```

Every R script is named after the analysis it performs and carries a header with its
purpose, the manuscript figure, its main output files and its original path. A
one-to-one table of script -> figure -> original file name is in
[docs/03_file_mapping.md](docs/03_file_mapping.md).

## Documentation

| Document | Content |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | End-to-end workflow, stage by stage |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | What lives where, and the naming conventions |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | Every script: previous path -> new path, purpose, figure |
| [docs/04_data_availability.md](docs/04_data_availability.md) | Working directory and input tables of every R script |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | Cluster setup, conda environments, R packages |

## Reproducing the analyses

**1. Upstream pipelines.** The scripts under `pipelines/` are the Slurm job scripts used
for the production runs (`#SBATCH` headers included). The site-specific absolute paths
have been replaced by variables: export `LUSTRE_HOME` and `DB_ROOT`, adapt the queue
names, then submit in the order given in [docs/01](docs/01_pipeline_overview.md).

**2. R analyses.** The scripts require R >= 4.2 and the packages listed in
[docs/05](docs/05_software_and_environment.md). Input tables are read relative to
`GTR_DATA`, which defaults to the current working directory:

```bash
cd /path/to/06_predatory_bacteria_virome
Rscript /path/to/Predatory-bacterial-virome/analyses/02_phage_catalog/01_viral_genome_info.R
# or
Rscript -e 'Sys.setenv(GTR_DATA = "/path/to/data"); source("analyses/02_phage_catalog/01_viral_genome_info.R")'
```

Some scripts read the earlier working copy of the project instead (`GTR_WORK`), and one
dot-plot variant reads a shared defensome folder (`GTR_DEF`); the header of each script
lists the variables it needs.

## Data availability

The repository is code-only. The myxobacterial genome collection, the viral catalogue and
the annotation tables are provided as Supplementary Data of the paper; reads of the newly
recovered MAGs are available from the NCBI SRA under the accession given in the paper.
The directory layout the scripts expect is documented in
[docs/04_data_availability.md](docs/04_data_availability.md).

## Note on the previous layout

The previous version of this repository was a flat list of shell scripts plus an `R code/`
folder. It has been reorganized into the layered structure above, and **updated from the
current working copy**: 15 R scripts that were missing from the old version were added
(candidate-phage gene maps, CRISPR-Cas typing, phage/plasmid metadata and geography, the
Cas12m tree, host defense systems and the dot-plot variants), plus the GFF -> GenBank
utility that had not been published. Three superseded drafts and two empty `try`
placeholders were dropped.

No absolute path, user home directory or local Windows folder remains in the repository:
everything is handled by the `GTR_DATA`, `GTR_WORK`, `GTR_DEF`, `LUSTRE_HOME` and
`DB_ROOT` variables.
