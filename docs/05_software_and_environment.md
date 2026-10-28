# Software and environment

## Cluster / job scheduler

All `pipelines/` scripts were submitted with **Slurm** (`#SBATCH` headers preserved).
The production runs used two queues, `cpu` and `fat`; job names, memory requests and
thread counts are kept in the scripts as a record of the actual runs.

The site-specific absolute paths are variables, declared in the
`# Site-specific settings` block just below each `#SBATCH` header:

| Variable | Meaning | Default in the script |
| --- | --- | --- |
| `LUSTRE_HOME` | home directory on the large-storage filesystem | `/path/to/lustre/home` |
| `DB_ROOT` | reference databases (geNomad, CheckV, BLAST, HH-suite, NCBI taxonomy) | `${LUSTRE_HOME}/database` |

To re-run a script elsewhere:

```bash
export LUSTRE_HOME=/your/lustre/home
export DB_ROOT=/your/databases
sbatch pipelines/05_defense_systems/17_defensefinder.sh
```

Also adapt `#SBATCH --partition=...` and replace `source activate <env>` with
`conda activate <env>` on modern conda.

## Conda environments (as used)

| Environment | Used by | Software |
| --- | --- | --- |
| `prokka` | `01_prodigal_genomes.sh` | Prodigal (meta mode) |
| `genomad` | `02_genomad_genomes.sh` | geNomad |
| `vs2` | `03_virsorter2.sh` | VirSorter2 |
| `checkv` | `04_checkv.sh` | CheckV (db 1.5, `CHECKVDB`) |
| `mmseqs2` | `05_mmseqs_cluster_phages.sh` | mmseqs2 |
| *(comparem)* | `06_comparem_aai.sh` | CompareM |
| `blast2.16` | `07_blastp_nr.sh`, `15_blastn_spacers_imgvr.sh` | BLAST+ 2.16 |
| `eggnog` | `08_diamond_nr_prophage.sh` | DIAMOND |
| *(taxdump)* | `09_diamond_taxonomy.sh` | NCBI prot.accession2taxid + taxdump |
| `hhsuite3` | `10_hhsuite_plasmid.sh` | HH-suite 3 |
| *(CLANS)* | `11_clans_djr_mcp.sh` | CLANS (`python -m clans`) |
| `crisprcasfinder` | `12_crisprcasfinder.sh` | CRISPRCasFinder |
| `cctyper` | `13_crisprcaster.sh` | CRISPRCasTyper |
| `defensefinder` | `17_defensefinder.sh` | DefenseFinder (`--preserve-raw`) |

Databases (geNomad, CheckV, NR, HH-suite, IMG/VR, taxdump) must be downloaded separately;
their locations are the `DB_ROOT`-based paths in the scripts.

## R

The analyses are plain R scripts (no `renv`, no project file). They were developed under
R 4.x; all 25 scripts parse cleanly under R 4.6.0. Packages used, by frequency:

* data handling: `dplyr`, `tidyr`, `tidyverse`, `stringr`, `data.table`, `tibble`, `readxl`
* plotting: `ggplot2`, `cowplot`, `Cairo`, `ggpubr`, `ggsankeyfier`, `ggsignif`, `ggsci`,
  `gghalves`, `ggbeeswarm`, `MetBrewer`, `RColorBrewer`, `viridis`, `scales`, `pheatmap`,
  `ComplexHeatmap`, `circlize`, `grid`
* trees / trees annotation: `ape`, `itol.toolkit`
* statistics: `rstatix`

A minimal install:

```r
install.packages(c("dplyr", "tidyr", "tidyverse", "stringr", "ggplot2", "cowplot",
                   "Cairo", "ggpubr", "rstatix", "pheatmap", "RColorBrewer",
                   "viridis", "scales", "readxl", "data.table"))
install.packages(c("itol.toolkit", "ggsankeyfier", "ggsignif", "ggsci", "gghalves",
                   "ggbeeswarm", "MetBrewer", "ComplexHeatmap", "circlize"))
```

`Cairo` needs the system `libcairo`; `ComplexHeatmap` comes from Bioconductor.

## Reproducibility notes

* Scripts are self-contained: no shared helper file, no `source()` between analyses.
* Figures are written as PDF with explicit `width`/`height` in inches; several scripts
  write the underlying tables (`.RData`) next to the input data.
* `GTR_DATA`, `GTR_WORK` and `GTR_DEF` are resolved with `normalizePath()` at load time,
  so the `setwd()` calls inside the scripts cannot invalidate them.
