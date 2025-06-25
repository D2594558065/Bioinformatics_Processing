######################如果没有gaf文件
#####通过gaf文件进行非模式物种的GSEA
library(GO.db)
library(clusterProfiler)

goterms <- Term(GOTERM)
term2name <- data.frame("GOID"=names(goterms),"term"=goterms )
head(term2name)

gaf_file <- "/home/mnt1/wanghao/project/Whale_IPS/IPS/IPS_PCA/pluripotency/GAF/GCF_002837175.3-RS_2023_04_gene_ontology.gaf"
gaf_data <- read.delim(gaf_file, header = FALSE, comment.char = "!", stringsAsFactors = FALSE)

#########过滤出生物学过程 P (biological process), F (molecular function) or C (cellular component) 的哪种GO注释
gaf_data <- gaf_data   %>% 
    filter(V9 == "P")

gene_go_annotations <- gaf_data[, c(5, 3)]
colnames(gene_go_annotations) <- c("GOID","GeneSYMBOL")

library(GSEABase)

sperm <- read.table("/home/mnt1/wanghao/project/Whale_IPS/IPS/IPS_PCA/file_out/spermDevelop_DESeq2.out")  %>% 
          dplyr::select(log2FoldChange_EB10_vs_wipsc,padj_EB10_vs_wipsc,log2FoldChange_EB20_vs_EB10,padj_EB20_vs_EB10)
spermprotein <- read.table("/home/mnt1/wanghao/project/Whale_IPS/IPS/IPS_PCA/pluripotency/gtf/sperm_protein_coding_gene_names.txt")

sperm <- sperm[which(rownames(sperm) %in% spermprotein$V1),]


sperm <- sperm  %>% 
  arrange(-log2FoldChange_EB20_vs_EB10)

# 假设你的数据框是sperm
sperm <- subset(sperm, !is.na(log2FoldChange_EB20_vs_EB10))


genelist <-setNames(sperm$log2FoldChange_EB20_vs_EB10,rownames(sperm))

egmt <- GSEA(genelist,TERM2GENE = gene_go_annotations,TERM2NAME = term2name ,minGSSize = 1)
