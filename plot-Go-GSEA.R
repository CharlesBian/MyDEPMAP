#3. correlation plot
library(ggstatsplot)
ggscatterstats(data =exprSet, 
               y = GMEB2, 
               x = YTHDF1,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7", 
               yfill = "#009E73", 
               marginal.type = "histogram",
               title = "Relationship between GMEB2 and YTHDF1")

# Go and KEGG analysis for singificant gene

gene <-str_trim(cor_data_sig$SYMBOL,'both')
gene = bitr(gene, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")

go <- enrichGO(gene = gene$ENTREZID, 
               OrgDb = "org.Hs.eg.db", 
               ont = "BP")

barplot(go, split="ONTOLOGY")+ 
  facet_grid(ONTOLOGY~., scale="free")

dotplot(go, split="ONTOLOGY")+ 
  facet_grid(ONTOLOGY~., scale="free")

###KEGG分析  
de2<-gene$ENTREZID
EGG <- enrichKEGG(de2,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
head(EGG)
dotplot(EGG)

###########
#5 聚类分析

#ClusterProfiler analysisL GO.KEGG, Network, GSEA...
# https://blog.csdn.net/weixin_43358851/article/details/83378100
# https://www.jianshu.com/p/feaefcbdf986
# https://cloud.tencent.com/developer/article/1078364

library("DOSE")
library("clusterProfiler") # first load "DOSE" library
library(stringr)
library(org.Hs.eg.db)
library(enrichplot) # gseaplot2

#获得基因列表

library(stringr)

# generation of gelist file for GSEA

d2<- str_trim(cor_data_df$SYMBOL,'both')

d2 = bitr(d2, fromType="SYMBOL", 
          toType="ENTREZID", 
          OrgDb="org.Hs.eg.db")


d3 <- cor_data_df %>% 
  inner_join(d2,by="SYMBOL")%>%
  arrange(desc(correlation)) %>%# sorting by R of correlation
  drop_na()




## Reactome pathway enrichment analysis  

##GSEA analysis

##create genelist file for GSEA
geneList2 = d3[,2]  # R, p.value, or Fold change as genelist
names(geneList2) = as.character(d3[,4]) # gene ID as names.
geneList2 = sort(geneList2, decreasing = TRUE) # sorting according R, p.value, or Fold change as genelist



#GSEA-GO/KEGG mapping
kk1 <- gseGO(geneList     = geneList2,
             ont="ALL",
             OrgDb        = org.Hs.eg.db,
             nPerm        = 1000,
             minGSSize    = 20,
             maxGSSize=500,
             pvalueCutoff = 0.05,
             verbose      = TRUE)
head(kk1)
enrichplot::dotplot(kk1)
write.csv(kk1,"CRC_NCOA1_kk1.csv")


##### 

kk2 <- gseKEGG(geneList     = geneList2,
               organism     = 'hsa',
               nPerm        = 5000,
               minGSSize    = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
kk2[1:40,]

enrichplot::dotplot(kk2)
typeof(kk2)

barplot(kk2,showCategory=20)
#GSEA-KEGG Plot example
gseaplot(kk2,geneSetID="hsa04392",title="hsa04392-Hippo")
gseaplot2(kk2,geneSetID="hsa04392",title="hsa04392-Hippo") # exprt 500X500
gseaplot2(kk2,geneSetID="hsa04340", subplots = 1:5, base_size = 9,rel_heights = c(2, 0.5, 1),title="hsa04340-Hedgehog") # exprt 500X500


# Example:  https://www.kegg.jp/kegg-bin/show_pathway?hsa04340+339745+2932

sortkk2<-kk2[order(kk2$enrichmentScore,decreasing=T)]
head(sortkk2)
sortkk2[1:10,]
barplot(sortkk2[,1:4])
gseaplot2(sortkk2)

gseaplot2(kk2, row.names(sortkk2)[1:4], pvalue_table = TRUE)
gseaplot2(kk2,geneSetID="hsa03050",title="hsa04392-Hippo")

gseaplot2(kk2, row.names(sortkk2)[2], pvalue_table = TRUE)


#gseaplot2(x, geneSetID, by = "all", title = "",color = "black", color.line = "green", color.vline = "#FA5860", ...)

write.csv(kk2, "CRC_KDM4D_kk2.CSV")

### 6 other Plot for GSEA data
emapplot(kk2)
cnetplot(kk2, showCategory=10)# big file
browseKEGG(kk2,"hsa4340")# Hedghog

####STOP
