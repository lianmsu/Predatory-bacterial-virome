# =============================================================================
# Predatory-bacterial-virome -- phage catalog / Geographical distribution of the phages
# =============================================================================
# Link each phage to the isolation site of its host and summarise the global
# distribution by country and ecosystem.
#
# Manuscript      : Fig. 3D
# Main outputs    : 2_phage_geographical_distribution.pdf
# Original script : 1mobilome_predatory_bacteria/3vcontact2/2_phage_geographical_distribution.R
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


setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/3vcontact2"))

self_selected_phage <- read.delim2("self_selected_phage_node.tsv")
self_selected_phage <- self_selected_phage[1:1983,]
library(dplyr)
# 处理第一列：删除|及之后的部分
self_selected_phage <- self_selected_phage %>%
  mutate(Node_clean = sub("\\|.*", "", Node))

IMGVR_phage_metadata <- read.delim2("filtered_IMGVR_sequences.tsv")

self_selected_phage <- self_selected_phage %>%
  left_join(IMGVR_phage_metadata %>% select(UVIG, Ecosystem.classification),
            by = c("Node_clean" = "UVIG"))

all_4649MAG_metadata <- read.delim2("all_4649_genome_metadata.tsv")

library(stringr)
self_selected_phage$Node_clean <- str_replace(self_selected_phage$Node_clean, "-[^-]*$", "")

self_selected_phage <- self_selected_phage %>%
  left_join(all_4649MAG_metadata %>% select(Genome, source_category),
            by = c("Node_clean" = "Genome")) %>%
  mutate(Ecosystem.classification = ifelse(is.na(Ecosystem.classification),source_category,Ecosystem.classification) ) %>%
  select(-source_category)


self_selected_phage_2 <- self_selected_phage %>%
  mutate(
    Ecosystem.classification = case_when(
      # 海洋环境
      grepl("marine|seawater|ocean|hydrothermal vent|coral|sponge|sea water", 
            Ecosystem.classification, ignore.case = TRUE) ~ "Marine",
      
      # 淡水环境
      grepl("lake|river|pond|freshwater|stream|aquifer|groundwater", 
            Ecosystem.classification, ignore.case = TRUE) ~ "Freshwater",
      
      # 陆地环境
      grepl("soil|sediment|rhizosphere|root|plant|forest|grassland|sand|terrestrial", 
            Ecosystem.classification, ignore.case = TRUE) ~ "Terrestrial",
      
      # 宿主相关
      grepl("human|mouse|rat|pig|cow|host|gut|intestine|fecal|feces|Host-associated|Host_associated", 
            Ecosystem.classification, ignore.case = TRUE) ~ "Host_associated",
      
      # 工程环境
      grepl("bioreactor|wastewater|sludge|landfill|leachate|activated sludge|waste water|Engineered", 
            Ecosystem.classification, ignore.case = TRUE) ~ "Engineered",
      
      # 极端环境
      grepl("hot spring|geothermal|alkaline|acidic|halophilic|saline|extreme", 
            Ecosystem.classification, ignore.case = TRUE) ~ "Extreme",
      
      # 默认类别
      TRUE ~ "other"
    )
  )

category_summary <- self_selected_phage_2 %>%
  count(Ecosystem.classification) %>%
  arrange(desc(n))
print(category_summary)

library(ggplot2)
library(dplyr)

nature_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#79706E",
                   "#B07AA1", "#FFBE7D", "#8CD17D", "#86BCB6", "#499894", "#D37295")

# 水平条形图版本
p_virus_unclassified_MGE <- 
ggplot(category_summary, aes(x = n, y = reorder(Ecosystem.classification, n), fill = Ecosystem.classification)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",linewidth = 0.1) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    y = "",
    x = "Number of viruses and unclassified MGEs"
  ) +
  scale_fill_manual(values = c("Freshwater"="#4E79A7","Terrestrial"="#FFBE7D",
                               "Engineered"="#59A14F","Marine"="#00b4d8","other"="#79706E",
                               "Host_associated"="#E15759","Extreme"="#D37295"
                               )) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 12, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_virus_unclassified_MGE
##地理信息##########

IMG_genomeCart_metadata <- readxl::read_excel("IMG_genomeCart_metadata.xlsx")

IMGVR_phage_metadata <- IMGVR_phage_metadata %>%
  left_join(IMG_genomeCart_metadata %>% select(`IMG Genome ID`, Latitude, Longitude),
            by = c("Taxon_oid" = "IMG Genome ID"))

self_selected_phage_3 <- self_selected_phage_2 %>%
  left_join(IMGVR_phage_metadata %>% select(UVIG,  Latitude, Longitude),
            by = c("Node_clean" = "UVIG")) 

GTDB_genome_lat_lon <- read.delim2("gtdb_r226_predators_metadata_new_2.tsv")
GTDB_genome_lat_lon <- GTDB_genome_lat_lon[,c("accession","ncbi_lat_lon")]

library(dplyr)
library(stringr)
GTDB_genome_lat_lon <- GTDB_genome_lat_lon %>%
  mutate(
    Latitude = case_when(
      ncbi_lat_lon == "none" ~ "none",
      str_detect(ncbi_lat_lon, "S") ~ paste0("-", str_extract(ncbi_lat_lon, "^[0-9.]+")),
      str_detect(ncbi_lat_lon, "N") ~ str_extract(ncbi_lat_lon, "^[0-9.]+"),
      TRUE ~ NA_character_
    ),
    Longitude = case_when(
      ncbi_lat_lon == "none" ~ "none",
      str_detect(ncbi_lat_lon, "W") ~ paste0("-", str_extract(ncbi_lat_lon, "[0-9.]+(?= [EW])")),
      str_detect(ncbi_lat_lon, "E") ~ str_extract(ncbi_lat_lon, "[0-9.]+(?= [EW])"),
      TRUE ~ NA_character_
    )
  )

GTDB_genome_lat_lon <- GTDB_genome_lat_lon %>%
  mutate(accession = str_replace_all(accession, "\\.", "_"))

self_selected_phage_3 <- self_selected_phage_3 %>%
  left_join(
    GTDB_genome_lat_lon %>% select(accession, Latitude_new = Latitude, Longitude_new = Longitude),
    by = c("Node_clean" = "accession")
  ) %>%
  mutate(
    Latitude  = ifelse(is.na(Latitude_new), Latitude, Latitude_new),
    Longitude = ifelse(is.na(Longitude_new), Longitude, Longitude_new)
  ) %>%
  select(-Latitude_new, -Longitude_new)


GEM_genome_metadata <- read.delim2("GEM_genome_metadata.tsv")

self_selected_phage_3 <- self_selected_phage_3 %>%
  left_join(
    GEM_genome_metadata %>% select(genome_id, Latitude_new = latitude, Longitude_new = longitude),
    by = c("Node_clean" = "genome_id")
  ) %>%
  mutate(
    Latitude  = ifelse(is.na(Latitude_new), Latitude, Latitude_new),
    Longitude = ifelse(is.na(Longitude_new), Longitude, Longitude_new)
  ) %>%
  select(-Latitude_new, -Longitude_new)


write.csv(self_selected_phage_3, "self_selected_phage_3.csv")


#手动整理部分病毒序列的经纬度
self_selected_phage_4 <- readxl::read_excel("self_selected_phage_3.xlsx")
self_selected_phage_4 <- self_selected_phage_4[!self_selected_phage_4$Latitude %in% c("none","NA","n/a"),]
self_selected_phage_4 <- self_selected_phage_4[!self_selected_phage_4$Longitude %in% c("none","NA","n/a"),]
self_selected_phage_4 <- self_selected_phage_4 %>%
  mutate(
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

redundancy_count <- self_selected_phage_4 %>%
  select(Ecosystem.classification, Latitude, Longitude) %>%   # 只保留最后三列
  group_by(Ecosystem.classification, Latitude, Longitude) %>% # 按三列分组
  summarise(virus_number = n(), .groups = "drop")             # 计数

#参考https://cloud.tencent.com/developer/article/2206435
library(ggplot2)
world<-map_data("world")
world.map <- 
ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=Longitude,y=Latitude, 
                 fill = Ecosystem.classification, size = virus_number),
             shape = 21,      # 带边框的填充圆
             colour = "black",# 边框颜色
             stroke = 0.05     # 边框宽度（mm）
             )+
  scale_size(range = c(2,6))+
  scale_fill_manual(values = c("Freshwater"="#4E79A7","Terrestrial"="#FFBE7D",
                               "Engineered"="#59A14F","Marine"="#00b4d8","other"="#79706E",
                               "Host_associated"="#E15759","Extreme"="#D37295"
  )) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map


world.map2 <- 
  ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=Longitude,y=Latitude, 
                 color = Ecosystem.classification, size = virus_number),
  )+
  scale_color_manual(values = c("Freshwater"="#4E79A7","Terrestrial"="#FFBE7D",
                               "Engineered"="#59A14F","Marine"="#00b4d8","other"="#79706E",
                               "Host_associated"="#E15759","Extreme"="#D37295"
  )) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map2

library(Cairo)
library(cowplot)
#cairo_pdf(filename = "1_genome_metadata.pdf" , width = 8,height = 4,unit = 'in')
g1 <- plot_grid(world.map,p_virus_unclassified_MGE,
                labels = letters[2], ncol = 2,
                rel_widths = c(2, 1)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

ggsave("2_phage_geographical_distribution.pdf",g1,width = 8.46,height = 3.66,units = 'in')








