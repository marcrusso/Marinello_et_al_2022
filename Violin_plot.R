library(ggpubr)
require(tidyverse)
library(GSVA)
####Violin plot function
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
give.n <- function(x){
  return(c(y = median(x) + 0.1, label = length(x)))
}

plot_Metadati_gene<-function(gene){
  x <- t(tot_matrix_FPKMUQ_tumor)
  x <- x[, gene]
  x <- rownames_to_column(as.data.frame(x),"samples")
  x <- merge(Metadati_lung,x, by="samples")
  p <- ggviolin(x, x = "Histopathology", y = "x",palette = c("darkcyan","forestgreen", "gold3","darkorange3"), fill = "Histopathology", add = "boxplot", add.params = list(fill = "white")) +
    stat_summary(fun.data = give.n, fun = median,position = position_dodge(width = 0.75),col="navy",size=7,geom = "text", inherit.aes = T) +
    stat_compare_means(inherit.aes = T,symnum.args = symnum.args,ref.group = "SCLC",size=4,method = "wilcox.test")+ylab(gene)
  print(p)
}



# NFKB IRF target GSEA
list_targets <- list(lista_IRF3_new, lista_NFKB1_new)
names(lista_targets)[1] <- "IRF3_targets"
names(lista_targets)[2] <- "NFKB_targets"

# GSVA
geni_targets_gsva <- gsva(as.matrix(tot_matrix_FPKMUQ_tumor), list_targets, method="gsva", verbose=TRUE )

plot_Metadati_path2<-function(pathway){
  x <- t(geni_targets_gsva)
  x <- x[, pathway]
  x <- rownames_to_column(as.data.frame(x),"samples")
  x <- merge(Metadati_lung,x, by="samples")
  p <- ggviolin(x, x = "Histopathology", y = "x", color= "black",fill = "Histopathology", palette = c("darkcyan","forestgreen", "gold3","darkorange3"), add = "boxplot", add.params = list(fill = "white", width=0.07, color="black"), ylab=pathway )+
    stat_summary(fun.data = give.n, col="black", size=4, geom = "text", inherit.aes = T) + ylim(-1,1) + stat_compare_means(inherit.aes = T,symnum.args = symnum.args,ref.group = "SCLC",size=4,method = "wilcox.test")
  print(p)
}

plot_Metadati_path2("IRF3_targets")
plot_Metadati_path2("NFKB_targets")
