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


ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/LINE1_OR_Data/bulk-RNA/Lung/GSE57148/lungheat.pdf", plot = p, device = "pdf", width = 5, height = 10,limitsize = FALSE)

