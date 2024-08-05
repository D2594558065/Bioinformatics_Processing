gene_list <- bitr(or6c74$ENS, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)   


go_enrichment1 <- enrichGO(gene = gene_list$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",  # 可以是 "BP", "CC", "MF"
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = TRUE)

up <- as.data.frame(go_enrichment1)
down <- as.data.frame(go_enrichment)
up$Type <- "Upregulated"
down$Type <- "Downregulated"

drawdata <- rbind(up[1:10,],down[1:10,])


# 为 ggplot 设置因子以控制条目的显示顺序
drawdata$Description <- factor(drawdata$Description, levels = unique(drawdata$Description))
#######绘图
# 将Description列转换为因子，以保持顺序
drawdata <- drawdata %>% distinct(Description, .keep_all = TRUE)
drawdata$Description <- factor(drawdata$Description, levels = rev(drawdata$Description))

# 绘制图表
library(ggplot2)

p <- ggplot(drawdata, aes(x = -log10(p.adjust), y = Description, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "#CD5C5C", "Downregulated" = "#1f77b4")) +
  theme_minimal(base_size = 14) +  # 调整基本字体大小
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 1, angle = 45, hjust = 1),  # 调整 X 轴标签角度以改善可读性
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),
    plot.caption = element_text(size = 8),
    panel.background = element_rect(fill = "white", colour = "white"),  # 设置面板背景为白色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.border = element_blank()  # 移除面板边框
  ) +
  labs(
    title = "GO Pathway Enrichment",
    x = "-log10(p.adjust)",
    fill = "Type"
  )

# 显示图表
print(p)

# 绘制基础图形

ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/RNA_seq/plot/OR6C74_GO.pdf", plot = p, device = "pdf", width = 6, height = 6)


