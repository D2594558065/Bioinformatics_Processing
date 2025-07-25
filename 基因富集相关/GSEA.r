geingclock <- read.table("/home/mnt1/wanghao/project/help/CYcgas/fruit_fly/DEseq2_Output/ageclockGene.txt",header = TRUE)

####制造gmt文件
agegene <- ageingclock  %>% 
    dplyr::inner_join(geneType,by = c("fly_gene" = "V2"))



gset <- c("AgingClock gene","NA",agegene$V1)
gset <- gset%>% 
     as.data.frame() %>% 
    t()
write.table(gset,file = "/home/mnt1/wanghao/project/help/CYcgas/fruit_fly/DEseq2_Output/fruitfly2.gmt",sep = "\t",row.names = F,col.names = F,quote = F)




library(clusterProfiler)
library(gggsea)
library(ggplot2)


gmt <- read.gmt("/home/mnt1/wanghao/project/help/CYcgas/fruit_fly/DEseq2_Output/fruitfly2.gmt")

rnk <-  read.table("/home/mnt1/wanghao/project/help/CYcgas/fruit_fly/DEseq2_Output/N4F_VS_NF.rank",header = TRUE,sep =',')
genes_list_by_term <- split(gmt$gene, gmt$term)
# geneList = rnk$log2FoldChange #把foldchange按照从大到小提取出来
names(geneList) <- rnk$V1 #给上面提取的foldchange加上对应上ENTREZID

geneList <- na.omit(geneList)
gsea <- fgsea::fgsea(genes_list_by_term, geneList, nperm=1000)
# df2 <- gseaCurve(geneList, genes_list_by_term, gsea)


ggplot2::ggplot() + 
  geom_gsea(df2,labelsize = 18) + 
  theme_gsea(7)+ 
  theme_bw()+
  theme(strip.text = element_text(size = 20,face = 'italic'),#增大分面标题字体
        strip.background = element_rect(fill = 'white'))+
  xlab(bquote(italic('Rank')))+ylab(bquote(italic('Enrichment Score')))+
  theme(axis.text.x = element_text(size=10,angle = -30,face = 'plain',hjust = 0.5),
        axis.text.y = element_text(size=10,angle = 0,face = 'plain',vjust = 0.5))+ theme(legend.key = element_rect(size = 50))
