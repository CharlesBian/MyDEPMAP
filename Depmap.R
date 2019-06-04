library(tidyr)
library(readr)
library(tidyverse)
library(ggsci)
library(ggstatsplot)
library(cowplot)
library(openxlsx)
library(tibble)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(psych)
library(tidyr)
library(data.table)

#create the file_t(optional)
file <- read.csv("gene_effect_corrected.csv",header = TRUE,check.names = FALSE)
rownames(file)<-file[,1]
file <- file[-1]
file_t <- t(file)
write.csv(file_t,"gene_effect_corrected_2.csv")

#1.read the file
GECKO <-read.csv("gene_effect_corrected_2.csv",header = TRUE,row.names = "X")
RNAi <- read.csv("D2_combined_gene_dep_scores.csv",header = TRUE,row.names = "X")

#2.read the gene
a <- as.data.frame(GECKO)
b <- as.data.frame(RNAi)


#input a gene name or a GENEID
gene = "CTH"
another_gene <- "DNMT3L"

#3.substract the gene expression
fread("gene_effect_corrected_2.csv")

expr1 <- as.numeric(a[gene,])
expr2 <- as.numeric(b[gene,])

expr1 <- as.numeric(a[3609,])
expr2 <- as.numeric(b[3609,])

#4.plot the density photo
p<-ggplot()
p+geom_density(aes(expr1),fill="lightblue",show.legend = TRUE,linetype = "blank",alpha = 0.3)+
  geom_density(aes(expr2),fill="mediumpurple1",show.legend = TRUE,linetype = "blank",alpha =0.3) +
  xlab("Dependency Score")+
  ylab("")+
  labs(title="Dependent Cell Lines",
       subtitle = paste("Gene Name:",gene))+
  geom_vline(xintercept = c(-1,0),lty=c(2,7),col=c("red","black"),lwd=c(1.2,0.8))+
  geom_hline(yintercept = 0,lty=7,col="black",lwd=1)+
 # legend()+
  theme(axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

#5.correlation anaysis

gene_rownames <- rownames(a)

cor_data_df <- data.frame(gene_rownames)

#caculate correlation
Dep_cor <- function(x,y)
{
  z <-cor.test(x[1:558],y,type = "pearson")
  return(z$estimate)
}
cor_data_df[,2]<-apply(a, 1, Dep_cor,y=expr1)

#caculate pvalue
Dep_p <- function(x,y)
{
  z <-cor.test(x[1:558],y,type = "pearson")
  return(z$p.value)
}
cor_data_df[,3]<-apply(a, 1, Dep_p, y=expr1)

names(cor_data_df) <- c("SYMBOL","correlation","pvalue")

#6.filter the data

cor_data_posi <- cor_data_df %>% 
  
  filter(pvalue < 0.05) %>% 
  filter(correlation > 0.15) %>%
  arrange(desc(correlation))%>% 
  dplyr::slice(1:1000)

View(cor_data_posi)

cor_data_neg <- cor_data_df %>% 
  filter(pvalue < 0.05) %>% 
  filter(correlation < -0.15) %>%
  arrange(correlation)%>% 
  dplyr::slice(1:1000)

View(cor_data_neg)
#7.plot the correlation photo

another_expr <- as.numeric(a[another_gene,])

corr <- cor.test(expr1,another_expr)

ggplot(data = NULL,
       aes(x=expr1,y=another_expr))+
  labs(title = paste("Correlation Plot:",gene,"VS",another_gene),
       subtitle=paste("Correlation number =",corr$estimate,"\n","Pvalue =",corr$p.value))+
  xlab(gene)+
  ylab(another_gene)+
  geom_point(size=1.5,alpha=3/10)+
  geom_smooth(method="lm",se=TRUE)+
  theme_bw()+ 
  scale_color_aaas()


write.csv(cor_data_posi,paste(gene,"positive correlation.csv"))

write.csv(cor_data_neg,paste(gene,"negative correlation.csv"))






