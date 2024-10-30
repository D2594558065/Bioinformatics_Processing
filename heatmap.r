###普通的热图
p1 <- pheatmap(data, 
         scale  = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Gene Expression Heatmap",
         cluster_rows = TRUE,
         cluster_cols = FALSE)






p1 <- Heatmap(cor_matrix, 
        name = "mat", 
        col = col_fun,  # 应用颜色函数
        show_row_names = FALSE,  # 不显示行名
        show_column_names = FALSE,  # 不显示列名
        top_annotation = ha)  # 添加箱线图注释


require(ggplotify)
p1 <- as.ggplot(p1)
ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/LINE1_OR_Data/bulk-RNA/Blood/GSE190125/Bloodheat_box.pdf", plot = p1, device = "pdf", width = 12, height = 16)



############
p <- pheatmap(cor_matrix,
         display_numbers = FALSE, # 如果你想在热图上显示相关性数值
         main = "Gene Correlation Heatmap", # 主标题
         color = colorRampPalette(c("blue", "white", "red"))(length(color_breaks) - 1), # 从蓝到白再到红的颜色渐变
         cluster_rows = TRUE, # 聚类行
         cluster_cols = TRUE, # 聚类列
         show_rownames = FALSE, # 不显示行名
         show_colnames = FALSE, # 不显示列名
         border_color = NA, # 去除单元格边框
         breaks = color_breaks # 设置颜色断点确保从-1到1
         )


##############带注释的热图
##参考链接https://blog.csdn.net/qq_42090739/article/details/121788868
############################################################
####其实只要准备一个基因名，描述的dataframe 传入到pheatmap里面就行
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(tibble)  # 确保加载 tibble 包
df_expanded <- go_enrichment1[c(1,10),] %>%
  separate_rows(geneID, sep = "/")  %>% 
  dplyr::select(Description,geneID)

df_combined <- df_expanded %>%
  group_by(geneID) %>%
  summarise(Description = paste(unique(Description), collapse = ", ")) %>%
  ungroup()
df_combined <- as.data.frame(df_combined)
rownames(df_combined) <- df_combined$geneID
df_combined <- df_combined %>%
  dplyr::select(-geneID)

drawdata <- normalized_counts_new[which(normalized_counts_new$SYMBOL %in% df_expanded$geneID),]

rownames(drawdata) <- drawdata$SYMBOL

drawdata1 <- drawdata[,1:6]
colnames(drawdata1) <- c("WT_293_DOX","WT_297_DOX","WT_349_DOXA","H11_280_DOX","H11_335_DOX","H11_13_DOX")


library(ggplot2)
library(pheatmap)
# 分组信息

p1 <- pheatmap(drawdata1,cluster_rows = T,cluster_cols = F,
          color = colorRampPalette(c("navy","white","firebrick3"))(100),
          show_colnames = T, border_color = NA,scale = "row",annotation_row = df_combined)


ggsave(filename = "/home/mnt1/wanghao/project/help/XY_DOX/rowdata/MJ20240822392-MJ-R-20240905082-徐晓钰-真核转录组纯测序-15个样本/plot/cell_different_heatmap.pdf", 
        plot = p1, device = "pdf", width = 10, height = 8)


ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/LINE1_OR_Data/bulk-RNA/Lung/GSE57148/lungheat.pdf", plot = p, device = "pdf", width = 5, height = 10,limitsize = FALSE)

