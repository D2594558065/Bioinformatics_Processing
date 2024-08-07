#################载荷PCA图片

merged_df <- as.data.frame(merged_df)
merged_df <- merged_df[!duplicated(merged_df$gene_TE), ]
rownames(merged_df) <- merged_df$gene_TE
dfM <- merged_df[,-1]
condition <- factor(c(rep("OR51E2", 3), rep("OR51H1P", 3),rep("OR52B2", 3), rep("OR6C74", 3), rep("PC", 3)), 
                    levels = c("OR51E2", "OR51H1P","OR52B2","OR6C74","PC"))

# 创建列数据框架
colData <- data.frame(condition, row.names = colnames(dfM))
head(merged_df)

library(DESeq2)
merged_df <- as.data.frame(merged_df)
merged_df <- merged_df[!duplicated(merged_df$gene_TE), ]
rownames(merged_df) <- merged_df$gene_TE
dfM <- merged_df[,-1]
# 转换 dfM 中的数据为整数型，处理 NA 值
dfM[] <- lapply(dfM, function(x) as.integer(x))
dfM <- dfM[!rowSums(is.na(dfM)), ]  # 删除包含 NA 的行

# 创建 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(countData = dfM, colData = colData, design = ~ condition)

# 运行 DESeq2 差异表达分析
dds <- DESeq(dds)

vst_data <- vst(dds, blind=FALSE)
vst_matrix <- assay(vst_data)

vst_matrix <- vst_matrix[69224:nrow(vst_matrix), ]

sample_categories <- c(rep("OR51E2",3),rep("OR51H1P",3),rep("OR52B2",3),rep("OR6C74",3),rep("PC",3))


library(ggplot2)

# 进行PCA
pca_result <- prcomp(t(vst_matrix), center=TRUE, scale.=FALSE) #注意我们对转置的矩阵进行PCA，因为我们想对样本进行PCA

# 可视化PCA结果
explained_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

pca_data <- as.data.frame(pca_result$x)

# 将样本类别添加到PCA数据框中
pca_data$category <- sample_categories


# 提取 PCA 载荷
loadings <- pca_result$rotation[, 1:2]  # 提取前两个主成分的载荷

# 创建载荷向量的数据框
loadings_df <- data.frame(Label = rownames(loadings), PC1 = loadings[, "PC1"], PC2 = loadings[, "PC2"])  %>% 
  arrange(-abs(PC1))


library(ggplot2)
library(grid)  # 可能需要用到，用于调整箭头样式

# 绘制 PCA 点和载荷向量
p <- ggplot() +
  geom_point(data = pca_data, aes(x = PC1, y = PC2, color = category), size = 5) +
  geom_segment(data = loadings_df[1:20, ], aes(x = 0, y = 0, xend = PC1*20, yend = PC2*20), 
               arrow = arrow(type = "closed", length = unit(0.02, "inches")), color = "gray") +
  geom_text(data = loadings_df[1:20, ], aes(x = PC1*20, y = PC2*20, label = Label), 
            hjust = -0.1, vjust = -0.1, check_overlap = TRUE) +
  labs(title = "PCA Plot with Loadings",
       x = paste("PC1: ", explained_var[1], "% variance"),
       y = paste("PC2: ", explained_var[2], "% variance")) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = NA)) +
  scale_color_brewer(palette = "Set1")  # 使用预设颜色板来区分不同的类别

ggsave(filename = "/home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/RNA_seq/plot/TE_PCA.pdf", plot = p, device = "pdf", width = 8, height = 8)
