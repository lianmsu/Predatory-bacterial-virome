# =============================================================================
# Predatory-bacterial-virome -- phage annotation / Gene maps of the 17 candidate phages
# =============================================================================
# Gene-arrow maps of the 17 newly identified candidate phages, one panel per
# phage.
#
# Manuscript      : Fig. 5
# Main outputs    : 17-phage gene map figures
# Original script : 1mobilome_predatory_bacteria/7phage_gene_annotation/1_phage_genome_annotation.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
GTR_WORK <- normalizePath(Sys.getenv("GTR_WORK", unset = "."), mustWork = FALSE)


setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/7phage_gene_annotation"))

all_17_phage_gff <- read.delim2("Candidate_17_new_phages_treated.gff",
                                header = F, sep = '\t')

library(dplyr)
library(stringr)
all_17_phage_gff <- all_17_phage_gff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列
all_17_phage_gff <- all_17_phage_gff[,c(1,3,4,5,7,10)]

colnames(all_17_phage_gff) <- c("contig","type","start","stop","strand",
                           "name")

#type和name可以改为需要的内容（type决定颜色，name决定标签）

write.table(all_17_phage_gff,"all_17_phage_gff.tsv",sep = "\t",quote = F)


primary_anno <- read.delim2("phage_annotation_summary.tsv")

all_17_phage_gff$seq <- sprintf("seq%d", 1:nrow(all_17_phage_gff))


all_17_phage_gff <- all_17_phage_gff %>%
  left_join(primary_anno, by = "seq")

#write.table(all_17_phage_gff,"all_17_phage_gff_anno.tsv",sep = "\t",quote = F)

all_17_phage_gff_manually <- read.delim2("all_17_phage_gff_anno.tsv")

library(dplyr)
library(stringr)

# 假设你的数据框名为 all_17_phage_gff_manually
classified_gff <- all_17_phage_gff_manually %>%
  mutate(
    # 将原type列重命名为其他列（例如feature），因为我们即将覆盖type列
    # 如果不需要保留原CDS信息，可以省略这一步，直接创建新的type列
    feature = type, 
    # 开始根据name列的内容创建新的type列
    type = case_when(
      # 1.  Virion structure and morphogenesis
      # 查找如 capsid, tail, fiber, portal, terminase, virion, morphogenesis, head, collar, sheath, baseplate 等关键词
      str_detect(name, regex("MCP|major capsid|capsid|tail|fiber|portal|terminase|virion|morphogenesis|head|collar|sheath|baseplate|ATPase|tape|neck|prohead|Holin|Haemolysin", ignore_case = TRUE)) ~ "Virion structure and morphogenesis",
      
      # 2. DNA replication, repair, transcription or translation
      # 查找如 replicase, polymerase, helicase, nuclease, ligase, repair, transcriptase, transcription, translation, ribosom, recombinase, primase, ssb 等关键词
      # 注意：'Phosphoadenosine phosphosulfate reductase' 是代谢相关，通常不属于这一类，但根据你的要求，它被下面的“Other function”捕获
      #str_detect(name, regex("replicase|polymerase|helicase|nuclease|ligase|repair|transcriptase|transcription|translation|ribosom|recombinase|primase|ssb|acetyltransferase|methyltransferase|synthetase", ignore_case = TRUE)) ~ "DNA replication, repair, transcription or translation",
      # 特别注意：上面的规则也捕获了“acetyltransferase”和“methyltransferase”，根据你的要求，它们被分到了这一类。
      # 但根据你的要求，DNA甲基化应单独分类。所以我们需要更精确的规则。
      
      # 2b. 更精确地捕获“DNA methylation and glycosylation” (优先级更高，所以放在前面)
      # 但case_when按顺序执行，所以我们需要把更精确的规则（DNA甲基化）放在更通用的规则（转录翻译）前面？
      # 不，因为目前“methyltransferase”同时被两个规则匹配。我们需要先定义DNA甲基化。
      
      # 修正：我们先定义DNA甲基化和糖基化，这样它们就不会被下面的通用规则捕获
      str_detect(name, regex("DNA.*methyl|methyl.*DNA|glycosylation|Glycosyltransferase|Glucosyltransferase|Acetylgalactosaminyltransferase|Mannosyltransferase", ignore_case = TRUE)) ~ "DNA methylation and glycosylation",
      
      # 3. 重新运行第2条规则（现在更通用的“转录翻译”规则不会捕获DNA甲基化相关基因了）
      str_detect(name, regex("replicase|polymerase|helicase|nuclease|ligase|repair|transcriptase|transcription|translation|ribosom|recombinase|primase|ssb|synthetase|acetyltransferase|Ribonuclease|replication|DNA-binding|RNA-binding|antigen|HNH endonuclease|DNA binding|RNA binding|nucleotidyltransferase|endonuclease|Rnase|Dnase|Fur|RNA polymerase|exonuclease|initiation|Primase–polymerase", ignore_case = TRUE)) ~ "DNA replication, repair, transcription or translation",
      
      # 4. Other function
      # 查找已知是其他功能的，比如代谢酶（如 reductase）
      # 或者明确不属于上述任何一类的已知功能
      str_detect(name, regex("reductase|metabol|auxiliary|host takeover|moron", ignore_case = TRUE)) ~ "Other function",
      
      # 5. Unannotated protein
      # 如果name是空的、NA、或包含 unannotated, unknown, hypothetical, NA 等
      is.na(name) | name == "" | str_detect(name, regex("unannotated|unknown|hypothetical|^NA$", ignore_case = TRUE)) ~ "Unannotated protein",
      
      # 6. 默认情况：如果以上都不匹配，保留原来的type（比如CDS）或标记为“Other function”
      # 根据你的需求，这里我们设为“Other function”
      TRUE ~ "Other function"
    )
  )

# 查看修改后的数据的前几行
head(classified_gff)



write.table(classified_gff,"all_17_phage_gff_withtype.tsv",sep = "\t",quote = F)









#######

setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/7phage_gene_annotation"))

my_17_phage <- read.delim2("Candidate_17_new_phages.tsv")


phage_genomad <- read.delim2("nr_phage_rep_seq_virus_summary.tsv")

library(dplyr)
my_17_phage <- my_17_phage %>% 
  left_join(phage_genomad,by=c("seq_name"="seq_name"))


phage_protospacer <- read.delim2("filtered_candidate17.out",header = F)

phage_protospacer <- phage_protospacer[,c(1,2,15,16)]
colnames(phage_protospacer) <- c("host","contig","start","stop")
phage_protospacer$type <- "Protospacer"

library(stringr)
phage_protospacer <- phage_protospacer %>% 
  mutate(
    host      = str_remove(host, "-[^-]*$"),      # 删掉最后一个 “-” 及之后
    contig = str_remove(contig, "\\|.*")   # 删掉第一个 “|” 及之后
  )

phage_protospacer <- unique(phage_protospacer)
phage_protospacer$strand <- "+"

#host info

all_4649_genome_metadata <- read.delim2("all_4649_genome_metadata.tsv")

phage_protospacer <- phage_protospacer %>%
  left_join(all_4649_genome_metadata, by = c("host"="Genome"))

phage_protospacer <- phage_protospacer[,c(2:21,1)]

write.table(phage_protospacer,"phage_protospacer_anno.tsv",quote = F,
            sep = '\t')






