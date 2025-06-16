df_top10 <- df_sorted[1:15, ]

library(ggplot2)

p <- ggplot(df_top10, aes(x = reorder(SMTS, -AverageExpression), y = AverageExpression, fill = SMTS)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() +
  theme(
    text = element_text( size = 18, face = "bold"),  # 字体设置
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 黑色边框
    panel.background = element_blank()  # 移除默认背景
  ) +
  labs(
    title = "OR51E2 Expression in Different Tissues",
    x = "SMTS", 
    y = "Average Expression", 
    fill = "SMTS"
  )
ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/plot/OR51E2_expression.pdf", plot = p, device = "pdf", width = 8, height = 8)
