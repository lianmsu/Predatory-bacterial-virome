# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home
# ---------------------------------------------------------------------------
cd ${LUSTRE_HOME}/06mobilome_predators/70_spacer2IMGVR
mkdir -p ./blastdb_IMGVR

#运行以下命令查找重复的 ID：
grep "^>" IMGVR_all_nucleotides.fna | sort | uniq -d
#使用 seqkit 工具自动去重（保留第一条出现的序列）：
seqkit rmdup -s IMGVR_all_nucleotides.fna -o IMGVR_all_nucleotides_dedup.fna

#构建成功
nohup makeblastdb -in ${LUSTRE_HOME}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna -dbtype nucl -out ./blastdb_IMGVR/IMGVR.db -max_file_sz '4GB' > nohup_makeblastdb_IMGVR.log 2>&1 &

#报错序列id太长
nohup makeblastdb -in ${LUSTRE_HOME}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna -dbtype nucl -parse_seqids -out ./blastdb_IMGVR/IMGVR.db -max_file_sz '4GB' > nohup_makeblastdb_IMGVR.log 2>&1 &

#报错存在重复序列，
nohup makeblastdb -in ${LUSTRE_HOME}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna -dbtype nucl -out ./blastdb_IMGVR/IMGVR.db -hash_index -max_file_sz '4GB' > nohup_makeblastdb_IMGVR.log 2>&1 &

#sz服务器无法makeblastdb建库，在bb建好copy过来的。

mkdir -p ./blastdb_IMGPR

nohup makeblastdb -in ${LUSTRE_HOME}/06mobilome_predators/2IMG_PR/IMGPR_nucl.fna -dbtype nucl -out ./blastdb_IMGPR/IMGPR.db -hash_index -max_file_sz '4GB' > nohup_makeblastdb_IMGPR.log 2>&1 &

