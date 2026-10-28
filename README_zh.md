# Predatory-bacterial-virome

> **The virome of globally distributed myxobacteria** 论文的代码仓库
> （Li, P. & Ni, J.）

本仓库收录了研究的全部代码：上游生物信息学流程（从收集的黏细菌基因组到病毒目录及其注释），
以及下游 R 统计分析与绘图脚本。

## 目录结构

```
Predatory-bacterial-virome/
├── pipelines/          # 上游流程：Slurm 作业脚本 + conda 环境（shell / python）
│   ├── 01_genome_and_phage_prediction/    Prodigal、geNomad、VirSorter2
│   ├── 02_viral_quality_and_clustering/   CheckV、mmseqs2、CompareM AAI
│   ├── 03_functional_annotation/          BLASTP NR、DIAMOND、NCBI 物种注释、HH-suite、CLANS
│   ├── 04_crispr_and_host_links/          CRISPRCasFinder、CRISPRCasTyper、IMG/VR BLAST
│   ├── 05_defense_systems/                DefenseFinder
│   └── 06_utilities/                      GFF 转 GenBank
├── analyses/           # 下游分析：R 统计与绘图
│   ├── 01_genome_catalog/                 图 S1、图 2A
│   ├── 02_phage_catalog/                  图 3
│   ├── 03_phage_annotation/               图 5；图 S2、S3、S6
│   ├── 04_new_viral_groups/               图 4
│   ├── 05_marker_phylogeny/               图 S4、S5；图 6D
│   └── 06_host_function/                  图 1
├── misc/               # 仅作参考保留、不属于论文流程的脚本
└── docs/               # 流程、目录、文件对照表、数据与软件说明
```

每个 R 脚本按分析内容命名，头部写明用途、对应图件、主要输出文件与原始路径；完整的
"脚本 → 图件 → 原文件名"对照见 [docs/03_file_mapping.md](docs/03_file_mapping.md)。

## 如何复现

**1. 上游流程**：`pipelines/` 下是当时实际投递的 SLURM 作业脚本（保留 `#SBATCH` 头）。
与运行环境相关的绝对路径已改为变量：导出 `LUSTRE_HOME`、`DB_ROOT`，改好队列名，按
[docs/01](docs/01_pipeline_overview.md) 的顺序提交。

**2. R 分析**：需要 R ≥ 4.2 及 [docs/05](docs/05_software_and_environment.md) 中列出的 R 包。
输入数据通过 `GTR_DATA` 定位，默认取当前工作目录：

```bash
cd /path/to/06_predatory_bacteria_virome
Rscript /path/to/Predatory-bacterial-virome/analyses/02_phage_catalog/01_viral_genome_info.R
```

部分脚本读取项目早期工作副本 `GTR_WORK`；一个点图变体读取共享的 defensome 目录
`GTR_DEF`——每个脚本头部都注明了它需要的变量。

## 数据

仓库只放代码。黏细菌基因组集合、病毒目录与注释结果随论文以 Supplementary Data 提供；新
MAG 的测序数据在 NCBI SRA 的登录号见论文。脚本所需的数据目录结构见
[docs/04_data_availability.md](docs/04_data_availability.md)。

## 关于旧版仓库

旧版本是平铺的 shell 脚本加一个 `R code/` 目录。本次按上面的分层结构整理，并**根据当前
工作目录的代码做了更新**：补充了旧版缺失的 15 个 R 脚本（候选噬菌体基因图、CRISPR-Cas
分型、噬菌体/质粒元数据与地理分布、Cas12m 系统发育、宿主防御系统、点图变体等），以及此前
未发布的 GFF 转 GenBank 工具；删除了 3 个被取代的旧稿和 2 个空的 `try` 占位文件。

仓库中已不含任何绝对路径、个人主目录或本地 Windows 目录：全部由 `GTR_DATA`、`GTR_WORK`、
`GTR_DEF`、`LUSTRE_HOME`、`DB_ROOT` 变量控制。
