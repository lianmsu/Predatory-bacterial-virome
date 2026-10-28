# =============================================================================
# Predatory-bacterial-virome -- host function / Defense systems of the phage hosts
# =============================================================================
# Summarise the DefenseFinder systems of the 4,696 predator genomes and of
# the phage hosts, and plot the DS distribution.
#
# Manuscript      : Fig. 1A
# Main outputs    : 1_p_Host_DS.pdf
# Original script : 1mobilome_predatory_bacteria/10host_DS/1_Host_Defense.R
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

setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/10host_DS"))

All_DF_res <- read.delim2("all_4696_predators_genome_renamed_defense_finder_systems.tsv")

phage_host_info <- read.delim2("phage_protospacer_anno.tsv")
phage_host_info <- phage_host_info[!phage_host_info$contig %in% c("IMGVR_UViG_3300048838_000477",
                                                                  "IMGVR_UViG_3300042453_004253",
                                                                  "IMGVR_UViG_3300025936_000002"),]

library(dplyr)
Host_MAG <- table(phage_host_info$host) %>% as.data.frame()

library(stringr)
All_DF_res <- All_DF_res %>% 
  mutate(
    host  = str_remove(sys_id, "-[^-]*$")      # 删掉最后一个 “-” 及之后
  )

All_DF_res <- All_DF_res %>% 
  mutate(
    host  = str_remove(host, "-[^-]*$")      # 删掉最后一个 “-” 及之后
  )
#重复多次

Host_DF_res <- All_DF_res[All_DF_res$host %in% Host_MAG$Var1,]
Host_DF_res <- Host_DF_res[Host_DF_res$activity == "Defense",]
Host_DF_res_DS <- table(Host_DF_res$type) %>% as.data.frame()

#DS分布####
library(RColorBrewer)
# 设置要显示的具体数量
top_n <- 20
DStype_counts <- Host_DF_res %>%
  count(type, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(type = ifelse(row_number() > top_n, "Other", type)) %>%
  group_by(type) %>%
  summarise(count = sum(count)) %>%
  arrange(desc(count))

library(ggplot2)
# 2. 绘制柱状图
p_DS <-
  ggplot(DStype_counts, aes(x = reorder(type, -count), y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(top_n), "gray")) +
  
  labs(title = "Distribution of defense system types",
       x = "",
       y = "Count") +
  geom_text(aes(label = count), vjust = -0.3, size = 3.5) +  # 在柱子上方显示数量
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 10, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 10, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_DS

setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/10host_DS"))
ggsave("1_p_Host_DS.pdf", p_DS, width = 5.5, height = 5.2, units = 'in')








