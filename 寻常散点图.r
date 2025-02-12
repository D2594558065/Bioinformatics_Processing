#####常规的散点图

library(ggplot2)

p <- ggplot(data = drawdata, aes(x = OvsY_log2FoldChange, y = HvsL_log2FoldChange)) +
  geom_point(color = 'grey') +
  geom_smooth(method = "lm", color = 'lightblue', fill = 'lightblue') +  # 回归线和填充设置为浅蓝色
  ggpubr::stat_cor(method = "spearman", color = 'black') +  # 相关性文本为黑色
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )+
  labs(
    x = "OvsY log2FC",  # 设置x轴标签
    y = "HvsL log2FC"  # 设置y轴标签
  )
