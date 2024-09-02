###########################################生成韦恩图
library(VennDiagram)

list_of_sets <- list(
  pygmy = rownames(pygmystage1),
  sperm = rownames(spermstage1)
)
# 二维韦恩图
# 你之前定义的list_of_sets已经包含了所有需要的集合

# 生成韦恩图
venn.plot <- venn.diagram(
  x = list_of_sets,
  filename = NULL,
  output = TRUE,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.50,
  label.col = "black",
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("cornflowerblue", "darkorchid1"),
  cat.cex = 1,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = 0.2
)


# 载入grid包
library(grid)

# 开启一个新的图形页面
grid.newpage()

# 使用grid.draw绘制韦恩图grobs对象
grid.draw(venn.plot)

pdf(file="/home/mnt1/wanghao/project/Whale_IPS/IPS/IPS_PCA/PCA_plot/EB10_vswipsc_diferr_0.05_1.5_venn.pdf", width = 4, height = 4)
grid.draw(venn.plot)
dev.off()
