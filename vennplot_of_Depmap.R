library("VennDiagram")

venn1 <- read.csv(file.choose(),header = TRUE)
venn2 <- read.csv(file.choose(),header = TRUE)
venn3 <- read.csv(file.choose(),header = TRUE)
venn4 <- read.csv(file.choose(),header = TRUE)
venn5 <- read.csv(file.choose(),header = TRUE)

#Two
venn_ints <- as.character(intersect(venn1$SYMBOL,venn3$SYMBOL))

venn_uni <- as.character(union(venn1$SYMBOL,venn2$SYMBOL))


#Three
venn_result <- list(v1 = as.character(intersect(venn1$SYMBOL,venn4$SYMBOL)),
                    v2 = as.character(intersect(venn1$SYMBOL,venn3$SYMBOL)),
                    v3 = as.character(intersect(venn4$SYMBOL,venn3$SYMBOL)))

names(venn_result) <- c(paste(venn1$SYMBOL[1],"vs",venn2$SYMBOL[1]),
                        paste(venn1$SYMBOL[1],"vs",venn3$SYMBOL[1]),
                        paste(venn2$SYMBOL[1],"vs",venn3$SYMBOL[1]))

venn_result

#Five
venn_result <- list(v1 = as.character(intersect(venn1$SYMBOL,venn2$SYMBOL)),
                    v2 = as.character(intersect(venn1$SYMBOL,venn3$SYMBOL)),
                    v3 = as.character(intersect(venn1$SYMBOL,venn4$SYMBOL)),
                    v4 = as.character(intersect(venn1$SYMBOL,venn5$SYMBOL)),
                    v5 = as.character(intersect(venn2$SYMBOL,venn3$SYMBOL)),
                    v6 = as.character(intersect(venn2$SYMBOL,venn4$SYMBOL)),
                    v7 = as.character(intersect(venn2$SYMBOL,venn5$SYMBOL)),
                    v8 = as.character(intersect(venn3$SYMBOL,venn4$SYMBOL)),
                    v9 = as.character(intersect(venn3$SYMBOL,venn5$SYMBOL)),
                    v10 = as.character(intersect(venn4$SYMBOL,venn5$SYMBOL)))
names(venn_result) <- c(paste(venn1$SYMBOL[1],"vs",venn2$SYMBOL[1]),
                        paste(venn1$SYMBOL[1],"vs",venn3$SYMBOL[1]),
                        paste(venn1$SYMBOL[1],"vs",venn4$SYMBOL[1]),
                        paste(venn1$SYMBOL[1],"vs",venn5$SYMBOL[1]),
                        paste(venn2$SYMBOL[1],"vs",venn3$SYMBOL[1]),
                        paste(venn2$SYMBOL[1],"vs",venn4$SYMBOL[1]),
                        paste(venn2$SYMBOL[1],"vs",venn5$SYMBOL[1]),
                        paste(venn3$SYMBOL[1],"vs",venn4$SYMBOL[1]),
                        paste(venn3$SYMBOL[1],"vs",venn5$SYMBOL[1]),
                        paste(venn4$SYMBOL[1],"vs",venn5$SYMBOL[1]))



#Two venn
x<-list(gene1 = venn1$SYMBOL, gene2 = venn2$SYMBOL)
names(x) <- c(paste(venn1$SYMBOL[1]),paste(venn2$SYMBOL[1]))

venn.diagram(x,
             filename = "gene_venn.tif",
             col = "transparent",
             fill = c("cornflowerblue","darkorchid1"),
             label.col = "black",
             cat.col = c("blue","red"),
             cat.dist = c(0.03,0.03),
             cat.pos = c(-8,8),
             cat.cex = 1.5,
             cex = 2,
             main = "Depmap genes cross",
             main.cex = 2)

#Tree venn
x<-list(gene1 = venn1$SYMBOL, gene2 = venn4$SYMBOL, gene3 = venn3$SYMBOL)
names(x) <- c(paste(venn1$SYMBOL[1]),paste(venn2$SYMBOL[1]),paste(venn3$SYMBOL[1]))

venn.diagram(x,
             filename = "three_gene.tif",
             col = "transparent",
             fill =  c("cornflowerblue","darkorchid1","green"),
             label.col = "black",
             cat.col = c("blue","maroon3","darkgreen"),
             cat.pos = c(-27,27,135),
             cat.cex = 1.5,
             cex = 2,
             main = "Genes venn plot",
             main.cex = 2)

#four venn
x<-list(gene1 = venn1$SYMBOL, gene2 = venn2$SYMBOL, gene3 = venn3$SYMBOL, gene4 = venn4$SYMBOL)
names(x) <- c(paste(venn1$SYMBOL[1]),paste(venn2$SYMBOL[1]),paste(venn3$SYMBOL[1]),paste(venn4$SYMBOL[1]))

venn <-venn.diagram(x,
                    filename = "four_gene.tif",
                    col = "transparent",
                    fill =  c("cornflowerblue","darkorchid1","green","yellow"),
                    label.col = "black",
                    cat.col = c("blue","red","darkgreen","orange"),
                    cat.pos = c(0,0,0,0),
                    cat.cex = 1.5,
                    cex = 2,
                    main = "Genes venn plot",
                    main.cex = 2)

#five venn
x<-list(gene1 = venn1$SYMBOL, gene2 = venn2$SYMBOL, gene3 = venn3$SYMBOL, gene4= venn4$SYMBOL, gene5 = venn5$SYMBOL)
names(x) <- c(paste(venn1$SYMBOL[1]),paste(venn2$SYMBOL[1]),paste(venn3$SYMBOL[1]),paste(venn4$SYMBOL[1]),paste(venn5$SYMBOL[1]))

venn.diagram(x,
             filename = "five_gene.tif",
             height = 4500,
             width = 4500,
             col = "transparent",
             fill =  c("cornflowerblue","darkorchid1","green","yellow","wheat"),
             label.col = "black",
             cat.col = c("blue","red","darkgreen","orange","wheat4"),
             cat.pos = c(0,0,-135,135,0),
             cat.cex = 1.5,
             cex = 2,
             main = "Genes venn plot",
             main.cex = 2)



























