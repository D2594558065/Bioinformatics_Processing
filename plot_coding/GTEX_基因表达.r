###############################脚本文件存放在
/home/mnt1/wanghao/project/help/GTEX/plot.r




library(data.table)
library(dplyr)
setwd("/home/mnt1/wanghao/project/help/基因表达临床融合数据/GTEx")
data <- fread('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct')
phe <- read.csv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep = '\t')  %>% 
    dplyr::select(SAMPID,SMTS)
dath <- read.csv("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",sep = '\t') 
# 假设您的数据框名为df，第一列名为V1

##########
df1 <-  data[grepl("^OR\\d+[A-Za-z]", data$Description), ]

df1 <- data.frame(df1)
df1 <- df1[!duplicated(df1$Description), ]
rownames(df1) <- df1$Description
df1 <- df1[,-c(1,2)]
df2 <- data.frame(t(df1)) #####转换为数据框类型十分的重要


row_names <- gsub("\\.", "-", rownames(df2))
df2$SAMPID <- row_names

df3 <- df2  %>% 
  inner_join(phe,'SAMPID')

sampid <- sapply(strsplit(df3$SAMPID, "-"), function(x) paste(x[1:2], collapse = "-"))
df3$SUBJID <- sampid
age <- sapply(strsplit(dath$AGE, "-"), function(x) mean(as.numeric(x)))
dath$age <- age
df4 <- df3  %>% 
  inner_join(dath,'SUBJID') 

# First, identify the columns that contain gene expression data
gene_expression_columns <- grep("^OR", names(df4), value = TRUE)

library(tidyr)
# 将基因表达数据转换为长格式
df_long <- df4 %>%
  pivot_longer(cols = gene_expression_columns, names_to = "Gene", values_to = "Expression")

# 按组织类型（假设为"SMTS"列）和基因分组，计算平均表达量
average_expression <- df_long %>%
  group_by(SMTS, Gene) %>%
  summarize(AverageExpression = mean(Expression, na.rm = TRUE)) %>%
  ungroup()

# 按组织类型和平均表达量排序
sorted_expression <- average_expression %>%
  arrange(desc(AverageExpression))

# 查看排序后的结果
print(sorted_expression)

 
sorted_expression[which(sorted_expression$Gene == "OR52B2"),]

df <- sorted_expression  %>% 
    filter(Gene == "OR51E2" )

library(ggplot2)

df_sorted <- df[order(-df$AverageExpression), ]

# 只选择前10个条目用于绘图
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
