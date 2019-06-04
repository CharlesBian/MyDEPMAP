genelist<-fread("result/Corrl_Hippo_genelist .csv")
#genelist <-read.csv("result/Corrl_Hippo_genelist .csv",header = TRUE,row.names = "X")

library(reshape2)

convert_data<-melt(genelist)

write.table(convert_data,"result/Corrl_Hippo_genelist_2.txt")

