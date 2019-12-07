if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("edgeR")

install.packages("ggfortify")
install.packages("RColorBrewer")

library(limma)
library(Glimma)
library(edgeR)
library(ggfortify)
library(RColorBrewer)
library(dplyr)

setwd("~/Downloads/510finalproject")
data <- read.table("count_table.txt", sep = "\t", stringsAsFactors = FALSE)
colname <- substr(data[1,],1,9)
names(data) <- colname
data <- data[-c(1),]


data <- data[-c(35365,35369),]

data <- cbind(genes, data)
data <- data[,-c(1,3)]



rownames(data) <- data$gene_symbol
data <- data[,-c(1)]

# Investigate on individual instead of sample
data_copy <- t(data)
data_copy <- as.data.frame(data_copy)
data_copy$sampleID <- rownames(data_copy)
data_copy <- data_copy[,c(50282,1:50281)]
data_individual <- merge(x=data_copy,y=samples, by.x="sampleID", by.y="rnaseq_profile_id")
data_individual <- data_individual[,c(1,50283,50284,2:50282)]
data_individual <- data_individual[,-c(1,3)]
data_individual <- data_individual[!duplicated(data_individual$donor_id),]
rownames(data_individual) <- data_individual$donor_id
data_individual <- data_individual[,-c(1)]


dm <- data.matrix(data)
dm <- as.data.frame(dm)
dge <- DGEList(counts=dm)

genes <- read.csv("rows-genes.csv", header = TRUE)
genes <- genes[,c(3,4)]

samples <- read.csv("columns-samples.csv", header = TRUE)
samples <- samples[,c(1,2,9)]

donor_info <- read.csv("new_DonorInformation.csv", header = TRUE)
donor_info <- donor_info[,c(1,4,5,20)]

total <- merge(samples, donor_info, by="donor_id")
sample_info <- total[,c(2,1,3,4,5,6)]
target <- samples$rnaseq_profile_id
sample_info <- sample_info[match(target, sample_info$rnaseq_profile_id),]

structure <- as.factor(sample_info$structure_acronym)
sex <- as.factor((sample_info$sex))
apo <- as.character(sample_info$apo_e4_allele)
apo <- replace(apo, apo=="N/A", "O")
apo <- as.factor(apo)
AD <- as.factor(sample_info$AD_condtion)

dge$samples$structure <- structure
dge$samples$sex <- sex
dge$samples$APO <- apo
dge$samples$AD <- AD

dge$genes <- genes

x <- dge

keep.exprs <- filterByExpr(x)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
x <- calcNormFactors(x, method = "TMM")
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)



par(mfrow=c(1,2))
col.apo <- apo
levels(col.apo) <-  brewer.pal(nlevels(col.apo), "Set1")
col.apo <- as.character(col.apo)
col.structure <- structure
levels(col.structure) <-  brewer.pal(nlevels(col.structure), "Set2")
col.structure <- as.character(col.structure)
plotMDS(lcpm, labels=apo, col=col.apo)
title(main="A. apo")
plotMDS(lcpm, labels=structure, col=col.structure, dim=c(3,4))
title(main="B. structure")

# design <- model.matrix(~0+apo+structure)
# colnames(design) <- gsub("apo", "", colnames(design))

# contr.matrix <- makeContrasts(N-Y, O-Y, N-O, levels = colnames(design))

# par(mfrow=c(1,2))
# v <- voom(x, design, plot=TRUE)
# vfit <- lmFit(v, design)
# vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
# efit <- eBayes(vfit)
# plotSA(efit, main="Final model: Mean-variance trend")


# pca_table <- t(x$counts)
# pca_table <- as.data.frame(pca_table)
# pca_table$AD <- AD
# pca_table_data <- pca_table[c(1:28485)]
# autoplot(prcomp(pca_table_data), data=pca_table, colour = "AD")

df_cpm <- as.data.frame(lcpm)
pca_table <- t(df_cpm)
pca_table <- as.data.frame(pca_table)
pca_table_data <- pca_table[c(1:28485)]

donor_id <- rownames(data_individual)

data_individual <- data.matrix(data_individual)
data_individual <- as.data.frame(data_individual)
full_table <- cbind(donor_id, data_individual)
full_table <- merge(x=full_table,y=donor_info[,c(1,4)], by.x="donor_id", by.y="donor_id")
rownames(full_table) <- donor_id
full_table <- full_table[,-c(1)]

autoplot(prcomp(data_individual), data=full_table, colour = "AD_condtion")
autoplot(prcomp(counts), data = full_table, color = "AD_condtion")
