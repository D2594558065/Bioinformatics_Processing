head(L1_result )
                baseMean log2FoldChange     lfcSE       stat       pvalue         padj
HAL1:L1:LINE   338.97557     0.90471529 0.1313271  6.8890215 5.617743e-12 6.300339e-11
HAL1M8:L1:LINE  31.69716    -0.13382053 0.3440152 -0.3889960 6.972791e-01 8.040920e-01
HAL1ME:L1:LINE 118.47863     0.66081129 0.1675267  3.9445135 7.996214e-05 3.097685e-04
HAL1b:L1:LINE   17.60660     0.05515429 0.3910374  0.1410461 8.878335e-01 9.359191e-01
L1HS:L1:LINE    40.24408     1.42630297 0.3173889  4.4938660 6.994165e-06 3.243617e-05
L1M1:L1:LINE    55.70585     0.50173091 0.3954796  1.2686645 2.045608e-01 3.322021e-01

输入数据为DESeq2的输出结果result
#####################################
L1_rows <- grepl("L1:LINE", rownames(result))

# 提取包含 'L1' 的行
L1_result <- result[L1_rows, ]

L1_result_filtered <- L1_result[rowMeans(L1_result, na.rm = TRUE) != 0, ]

library(ggplot2)
df <- L1_result_filtered

df$color <- ifelse(df$padj <= 0.05 & df$log2FoldChange > 0, "red",
                   ifelse(df$padj <= 0.05 & df$log2FoldChange < 0, "blue", "grey"))

# 绘制散点图
p <- ggplot(df, aes(x = row.names(df), y = log2FoldChange, color = color, size = baseMean)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # 添加 y = 0 的横线
  scale_color_identity() +
  scale_size_continuous(range = c(2, 20)) + # 设置点的大小范围
  labs(title = "Scatter Plot Highlighting LINE1 Transposons",
       x = "Transposon",
       y = "Log2 Fold Change") +
  theme_minimal(base_size = 15) + # 设置基础字体大小
  theme(
    panel.background = element_blank(),  # 设置白色背景
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    axis.text.x = element_blank(),  # 去除 X 轴标签
    axis.title.x = element_blank(), # 去除 X 轴标题
    axis.ticks.x = element_blank(), # 去除 X 轴刻度
    panel.border = element_rect(colour = "black", fill=NA, size=1) # 添加黑色方框
  )

ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/RNA_seq/plot/plot.pdf", plot = p, device = "pdf", width = 8, height = 7)

