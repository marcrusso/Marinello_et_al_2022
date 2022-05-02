require(tidyverse)
library(readxl)
library(GSVA)
library(msigdbr)
library(rms)
library(survminer)
library(grid)
library(gridExtra)
# Import matrix SCLC
library("DESeq2")
EGAD00001001244 <- list.files("EGAD00001001244")
EGAD00001001431 <- list.files("EGAD00001001431")
EGAD00001003801 <- list.files("EGAD00001003801")
samples<-"EGAD00001001244"
files <- file.path(samples, EGAD00001001244)
condition<-rep(1,length(EGAD00001001244))
sampleTable <- data.frame(sampleName = EGAD00001001244,
                          fileName = files, condition= condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,design= ~ 1)
EGAD00001001244_matrix<-counts(ddsHTSeq)
EGAD00001001244_matrix<-as.matrix(EGAD00001001244_matrix)

samples<-"EGAD00001001431"
files <- file.path(samples, EGAD00001001431)
condition<-rep(1,length(EGAD00001001431))
sampleTable <- data.frame(sampleName = EGAD00001001431,
                          fileName = files, condition= condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,design= ~ 1)
EGAD00001001431_matrix<-counts(ddsHTSeq)
EGAD00001001431_matrix<-as.matrix(EGAD00001001431_matrix)


samples<-"EGAD00001003801"
files <- file.path(samples, EGAD00001003801)
condition<-rep(1,length(EGAD00001003801))
sampleTable <- data.frame(sampleName = EGAD00001003801,
                          fileName = files, condition= condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,design= ~ 1)
EGAD00001003801_matrix<-counts(ddsHTSeq)
EGAD00001003801_matrix<-as.matrix(EGAD00001003801_matrix)

# Import  HTSeq-counts LUSC LUAD
Tcga_LUAD_htseq.counts <- read.delim("TCGA-LUAD.htseq_counts.tsv", stringsAsFactors = FALSE)
Tcga_LUSC_htseq.counts <- read.delim("TCGA-LUSC.htseq_counts.tsv", stringsAsFactors = FALSE)

# ConvertEnsembl_ID in rownames
Tcga_LUAD_htseq.counts<- column_to_rownames(Tcga_LUAD_htseq.counts, "Ensembl_ID")
Tcga_LUSC_htseq.counts<- column_to_rownames(Tcga_LUSC_htseq.counts, "Ensembl_ID")

# Convert counts in log2(x + 1) in counts non trasformate
Tcga_LUAD_counts<-apply(Tcga_LUAD_htseq.counts,2,function(x){(2^x)-1})
Tcga_LUSC_counts<-apply(Tcga_LUSC_htseq.counts,2,function(x){(2^x)-1})

# Merge matrix gencode V22
EGAD00001001244_matrix_gencode<- as.data.frame(EGAD00001001244_matrix)
EGAD00001001244_matrix_gencode<-rownames_to_column(EGAD00001001244_matrix_gencode, "gene")

EGAD00001001431_matrix_gencode<- as.data.frame(EGAD00001001431_matrix)
EGAD00001001431_matrix_gencode<- rownames_to_column(EGAD00001001431_matrix_gencode, "gene")

EGAD00001003801_matrix_gencode<- as.data.frame(EGAD00001003801_matrix)
EGAD00001003801_matrix_gencode <-rownames_to_column(EGAD00001003801_matrix_gencode, "gene")

EGAD00001001244_matrix_gencode <- merge(gencode.v22.annotation.gene[,c(1,2)], EGAD00001001244_matrix_gencode, by=1)
EGAD00001001431_matrix_gencode <- merge(gencode.v22.annotation.gene[,c(1,2)], EGAD00001001431_matrix_gencode, by=1)
EGAD00001003801_matrix_gencode <- merge(gencode.v22.annotation.gene[,c(1,2)], EGAD00001003801_matrix_gencode, by=1)

EGAD00001001244_matrix_gencode <- EGAD00001001244_matrix_gencode[2:length(colnames(EGAD00001001244_matrix_gencode))]
EGAD00001001431_matrix_gencode <- EGAD00001001431_matrix_gencode[2:length(colnames(EGAD00001001431_matrix_gencode))]
EGAD00001003801_matrix_gencode <- EGAD00001003801_matrix_gencode[2:length(colnames(EGAD00001003801_matrix_gencode))]

Tcga_LUAD_counts_gencode<- as.data.frame(Tcga_LUAD_counts)
Tcga_LUSC_counts_gencode<-as.data.frame(Tcga_LUSC_counts)

Tcga_LUAD_counts_gencode <- rownames_to_column(Tcga_LUAD_counts_gencode, "gene")
Tcga_LUSC_counts_gencode <- rownames_to_column(Tcga_LUSC_counts_gencode, "gene")

Tcga_LUAD_counts_gencode <- merge(gencode.v22.annotation.gene[,c(1,2)], Tcga_LUAD_counts_gencode, by=1)
Tcga_LUSC_counts_gencode <-merge(gencode.v22.annotation.gene[,c(1,2)], Tcga_LUSC_counts_gencode, by=1)

Tcga_LUAD_counts_gencode <- Tcga_LUAD_counts_gencode[2:length(colnames(Tcga_LUAD_counts_gencode))]
Tcga_LUSC_counts_gencode <- Tcga_LUSC_counts_gencode[2:length(colnames(Tcga_LUSC_counts_gencode))]


# Deduplicate
EGAD00001001244_matrix_gencode<- aggregate(.~gene, EGAD00001001244_matrix_gencode, sum)
EGAD00001001431_matrix_gencode<- aggregate(.~gene, EGAD00001001431_matrix_gencode, sum)
EGAD00001003801_matrix_gencode<- aggregate(.~gene, EGAD00001003801_matrix_gencode, sum)

Tcga_LUAD_counts_gencode <- aggregate(.~gene, Tcga_LUAD_counts_gencode, sum)
Tcga_LUSC_counts_gencode <- aggregate(.~gene, Tcga_LUSC_counts_gencode, sum)

# colnames
colnames(Tcga_LUAD_counts_gencode)<-gsub("\\.", "-", colnames(Tcga_LUAD_counts_gencode))
colnames(Tcga_LUSC_counts_gencode)<-gsub("\\.", "-", colnames(Tcga_LUSC_counts_gencode))

names(EGAD00001001244_matrix_gencode)<-gsub(pattern = ".count.txt*", replacement =  '', x=names(EGAD00001001244_matrix_gencode))
names(EGAD00001001431_matrix_gencode)<-gsub(pattern = ".count.txt*", replacement =  '', x=names(EGAD00001001431_matrix_gencode))
names(EGAD00001003801_matrix_gencode)<-gsub(pattern = ".count.txt*", replacement =  '', x=names(EGAD00001003801_matrix_gencode))
names(EGAD00001003801_matrix_gencode)<-gsub(pattern = "^LCNEC_", replacement =  '', x=names(EGAD00001003801_matrix_gencode))


# Transform genes in rownames
EGAD00001001244_matrix_gencode<- column_to_rownames(EGAD00001001244_matrix_gencode,"gene")
EGAD00001001431_matrix_gencode<- column_to_rownames(EGAD00001001431_matrix_gencode,"gene")
EGAD00001003801_matrix_gencode<- column_to_rownames(EGAD00001003801_matrix_gencode,"gene")

Tcga_LUAD_counts_gencode<- column_to_rownames(Tcga_LUAD_counts_gencode,"gene")
Tcga_LUSC_counts_gencode<-column_to_rownames(Tcga_LUSC_counts_gencode, "gene")

# Mrs_Matrix
tot_matrix <- do.call(cbind, list(EGAD00001001244_matrix_gencode,EGAD00001001431_matrix_gencode,EGAD00001003801_matrix_gencode,Tcga_LUAD_counts_gencode,Tcga_LUSC_counts_gencode))
tot_condition<-as.data.frame(pca_test[,"cluster"])
colnames(tot_condition)<-"condition"
colnames(tot_matrix)[60]<- "S00022_2"
rownames(tot_condition)<-colnames(tot_matrix)
tot_matrix<-as.matrix(tot_matrix)
tot_matrix<-apply(tot_matrix,2,as.integer)
rownames(tot_matrix)<-rownames(EGAD00001001244_matrix_gencode)
dds_lung_expression<- DESeqDataSetFromMatrix(countData = tot_matrix,
                              colData = tot_condition,
                              design = ~ condition)

vsd_lung_expression <- vst(dds_lung_expression, blind=FALSE)
head(assay(vsd_lung_expression[,1:3]))
plotPCA(vsd_lung_expression, intgroup=c("condition"),ntop=1000)
batchcorr_matrix<-limma::removeBatchEffect(x =assay(vsd_lung_expression), batch = tot_condition$condition)
vsd_lung_expression_batch<-vsd_lung_expression
assay(vsd_lung_expression_batch)<-batchcorr_matrix
plotPCA(vsd_lung_expression_batch, intgroup=c("condition"))
rld<-rlog(dds_lung_expression, blind=FALSE)


tot_condition_2 <- tot_condition
tot_condition_2[grep("EGAD", tot_condition_2$condition), "macro_condition"] <- "SCLC"
tot_condition_2[grep("TCGA", tot_condition_2$condition), "macro_condition"] <- "NSCLC"
dds_lung_expression_2<- DESeqDataSetFromMatrix(countData = tot_matrix,
                                             colData = tot_condition_2,
                                             design = ~ condition)
vsd_lung_expression_2 <- vst(dds_lung_expression_2, blind=FALSE)
head(assay(vsd_lung_expression_2[,1:3]))
batchcorr_matrix_2<-limma::removeBatchEffect(x =assay(vsd_lung_expression_2), batch = tot_condition_2$macro_condition)
vsd_lung_expression_batch_2<-vsd_lung_expression_2
assay(vsd_lung_expression_batch_2)<-batchcorr_matrix_2
plotPCA(vsd_lung_expression_batch_2, intgroup=c("condition"))


d <- plotCounts(dds_lung_expression_2, gene="IFNB1", intgroup="condition",
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) +
  geom_boxplot() +
  scale_y_log10(breaks=c(25,100,400))

FPKM_UQ_calc<- function(data){
  genelenghtlist<-gencode.gene.info.v22$full_length[match(rownames(data),gencode.gene.info.v22$gene_name)]
  apply(data,2,function(x){
    FPKM<-(x*1e+09)/(as.numeric(summary(x)[5])*genelenghtlist)
    return(FPKM)
  })
}

tot_matrix_FPKMUQ<-log2(FPKM_UQ_calc(tot_matrix)+1)
rm(genelenghtlist)
summary(tot_matrix_FPKMUQ[,1])
tot_matrix_FPKMUQ_t <- t(tot_matrix_FPKMUQ)
