
library(readr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(cowplot)


#1 import and process your data as:  row.name=Sample,col.name=Gene Symbol  ID
load("exprSet.Rdata")

#2 批量相关性分析, 

y <- as.numeric(exprSet[,"METTL3"])

colnames <- colnames(exprSet)

cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  
  test <- cor.test(as.numeric(exprSet[,i]),y,type="pearson")
  
  cor_data_df[i,2] <- test$estimate
  
  cor_data_df[i,3] <- test$p.value
  
}

names(cor_data_df) <- c("SYMBOL","correlation","pvalue")

#2,2 查看这个数据结构

head(cor_data_df)

#3.筛选最相关的基因
#筛选p值小于0.05，按照相关性系数绝对值选前1000个的基因， 数量可以自己定


library(dplyr)

library(tidyr)

cor_data_sig <- cor_data_df %>% 
  
  filter(pvalue < 0.05) %>% 
  
  arrange(desc(abs(correlation)))%>% 
  
  dplyr::slice(1:1000)

head(cor_data_sig) #查看这个数据结构

# 4.随机选取正的和负的分别作图验证
##4.1 ggplot2

library(ggplot2)
library(ggsci)
p1<-ggplot(data =exprSet, 
           aes(x=METTL3,
               y=METLL14
           ))+
  geom_point(size=1.5,alpha=3/10)+
  geom_smooth(method="lm",se=TRUE)+
  theme_bw()+ 
  scale_color_aaas()
p1
