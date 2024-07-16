library(GO.db)
library(clusterProfiler)

goterms <- Term(GOTERM)
term2name <- data.frame("GOID"=names(goterms),"term"=goterms )
head(term2name)

#####gaf文件可以在NCBI物种基因组注释FTP下载
gaf_file <- "/home/mnt1/wanghao/project/Whale_IPS/IPS/IPS_PCA/pluripotency/GAF/GCF_002837175.3-RS_2023_04_gene_ontology.gaf"
gaf_data <- read.delim(gaf_file, header = FALSE, comment.char = "!", stringsAsFactors = FALSE)

# #########过滤出生物学过程 P (biological process), F (molecular function) or C (cellular component) 的哪种GO注释
# gaf_data <- gaf_data   %>% 
#     filter(V9 == "P")

gene_go_annotations <- gaf_data[, c(5, 3)]
colnames(gene_go_annotations) <- c("GOID","GeneSYMBOL")

ego <- enricher(gene = result1$Symbol[1:200],
                TERM2GENE = gene_go_annotations,
                TERM2NAME = term2name,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
sperm_pygmy_both <- as.data.frame(ego)
spermOnly <- as.data.frame(ego)
# 以点图展示富集分析结果
dotplot(ego, showCategory = 20) + ggtitle("Dotplot of GO Enrichment Analysis")
