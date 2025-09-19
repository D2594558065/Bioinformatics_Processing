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
df1 <-  data[data$Description %in% c("POU5F1","TP53"),]
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

df <- df4[df4$SEX == 1,]
library(ggplot2)

# 假设你的数据框叫 df
# AGE 是字符型（比如 "60-69"），先确保它是因子并排序
df$AGE <- factor(df$AGE, levels = sort(unique(df$AGE)))


library(ggplot2)
library(ggpubr)

# 确保 AGE 是因子
df$AGE <- as.factor(df$AGE)

p <- ggplot(df, aes(x = .data$AGE, y = .data$POU5F1, fill = .data$AGE)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  facet_wrap(~SMTS, scales = "free_y") +
  theme_bw() +
  labs(
    title = "POU5F1 expression ",
    x = "Age group",
    y = "POU5F1 expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  # 添加方差分析 (ANOVA) 结果
  stat_compare_means(
    method = "anova",
    label = "p.format",
  )

print(p)


# ggsave(filename = "/home/mnt1/wanghao/project/help/ZYX/plot/OCT4_male.pdf", plot = p, device = "pdf", width = 8, height = 14)
ggsave(filename = "/home/mnt1/wanghao/project/help/ZYX/plot/OCT4_female.pdf", plot = p, device = "pdf", width =12, height = 14)



