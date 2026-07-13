setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/04_new_phage")

comparem_aai <- read.delim2("aai_summary.tsv",check.names = F)

# 移除 "all_48_phage.part_" 前缀
comparem_aai$`Genome A` <- gsub("39_phage\\.part_", "", comparem_aai$`Genome A`)

# 同样处理 Genome B 列（如果需要）
comparem_aai$`Genome B` <- gsub("39_phage\\.part_", "", comparem_aai$`Genome B`)

# 删除 "__" 及之后的所有字符
comparem_aai$`Genome A` <- sub("__.*", "", comparem_aai$`Genome A`)
comparem_aai$`Genome B` <- sub("__.*", "", comparem_aai$`Genome B`)

comparem_aai$`Genome A` <- gsub("^viptree_selected_", "", comparem_aai$`Genome A`)
comparem_aai$`Genome B` <- gsub("^viptree_selected_", "", comparem_aai$`Genome B`)

comparem_aai$`Orthologous fraction (OF)` <- 
  as.numeric(comparem_aai$`Orthologous fraction (OF)`)

# 获取所有基因组
all_genomes <- unique(c(comparem_aai$`Genome A`, comparem_aai$`Genome B`))

# 创建空矩阵
of_matrix <- matrix(0, nrow = length(all_genomes), ncol = length(all_genomes),
                    dimnames = list(all_genomes, all_genomes))

# 填充已知的OF值（转换为百分比）
for (i in 1:nrow(comparem_aai)) {
  a <- comparem_aai$`Genome A`[i]
  b <- comparem_aai$`Genome B`[i]
  of_value <- comparem_aai$`Orthologous fraction (OF)`[i] 
  
  of_matrix[a, b] <- of_value
  of_matrix[b, a] <- of_value  # 保持对称
}

# 设置对角线为100
diag(of_matrix) <- 100

# 查看矩阵维度
cat("矩阵维度:", dim(of_matrix), "\n")

# 查看部分结果
print(of_matrix[1:5, 1:5])


library(ggplot2)
library(pheatmap)
p0 <-
  pheatmap(of_matrix,
           fontsize = 8,#字体大小
           # fontface="italic",#斜体,其中'font'和'fontface'两个量只能设定一个
           border=T,#边框
           border_color = "gray100",          # 边框颜色
           #cluster_rows = F,
           #cluster_cols = F,
           treeheight_row = 0,
           treeheight_col = 0,
           # color = c("#DDDDDD", colorRampPalette(c("#EFA5CA", "#AAE4E4", "#7272FF"))(50)),  # 颜色向量
           # breaks = c(-Inf, 0, seq(75, 100, length.out = 51))  # 断点设置
           color = colorRampPalette(c( "#e7ecef","#AAE4E4", "#530050"))(100)#字体颜色
  )#导出为图片
p0
# 
# # 获取聚类后的行顺序
# clustered_row_names <- rownames(ani)[p0$tree_row$order]
# print(clustered_row_names)
# 
# ggsave(plot=p0,'1_48phage_AAI.pdf',width = 7.5,height = 7,units = 'in')

############################################

library(ggplot2)

library(ComplexHeatmap)
library(circlize)
library(grid)

#col_fun <- colorRamp2(
#  c(min(of_matrix, na.rm = TRUE), 50, max(of_matrix, na.rm = TRUE)),
#  c("#e7ecef", "#AAE4E4", "#7272FF")
#)

col_fun <- colorRamp2(
  c(min(of_matrix, na.rm = TRUE), 50, max(of_matrix, na.rm = TRUE)),
  c("#f7f4f9", "#df65b0", "#67001f")
)

ht <- Heatmap(
  of_matrix,
  col = col_fun,
  #cluster_rows = FALSE,
  #cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,   # ← 关键参数
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  #rect_gp = gpar(col = NA),   # 关闭默认边框
  # 关键：默认所有格子白色边框
  rect_gp = gpar(
    col = "white",
    lwd = 0.5
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(of_matrix[i, j]) && of_matrix[i, j] > 10) {
      grid.rect(
        x = x, y = y,
        width = w, height = h,
        gp = gpar(
          fill = NA,
          col = "black",
          lwd = 0.8
        )
      )
    }
  }
)

draw(ht)

# 打开PDF设备
pdf("1_ht1.pdf", width = 7.97, height = 5.17)
# 绘制热图
draw(ht)  # 重要：必须使用draw()函数！
# 关闭设备
dev.off()




###############################################################################
viptree_inter_simi <- read.delim2("VIRIDIC_sim-dist_table.tsv",row.names = "genome")

viptree_inter_simi <- as.matrix(viptree_inter_simi) # 先变成字符矩阵
viptree_inter_simi <- apply(viptree_inter_simi, 2, function(x) {
  as.numeric(gsub(",", ".", x))
})

rownames(viptree_inter_simi) <- colnames(viptree_inter_simi)

col_fun2 <- colorRamp2(
  c(min(viptree_inter_simi, na.rm = TRUE), 50, max(viptree_inter_simi, na.rm = TRUE)),
  c("#e7ecef", "#AAE4E4", "#7272FF")
)

ht2 <- Heatmap(
  viptree_inter_simi,
  col = col_fun2,
  #cluster_rows = FALSE,
  #cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,   # ← 关键参数
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  #rect_gp = gpar(col = NA),   # 关闭默认边框
  # 关键：默认所有格子白色边框
  rect_gp = gpar(
    col = "white",
    lwd = 0.5
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(viptree_inter_simi[i, j]) && viptree_inter_simi[i, j] > 70) {
      grid.rect(
        x = x, y = y,
        width = w, height = h,
        gp = gpar(
          fill = NA,
          col = "black",
          lwd = 0.8
        )
      )
    }
  }
)

ht2

# 打开PDF设备
pdf("1_ht2.pdf", width = 7.97, height = 5.17)
# 绘制热图
draw(ht2)  # 重要：必须使用draw()函数！
# 关闭设备
dev.off()







save.image("./1_aai_heatmap.RData")
load("./1_aai_heatmap.RData")
