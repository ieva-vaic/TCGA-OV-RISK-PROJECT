#THIS IS TCGA-OV-RISK-GENES project script No. 7
#normalize mRNA counts
#plan: get joined tcga/gtex count data -> outlier chek with hplots -> lognormalize

# Load packages ##########################################
library(tidyverse)
library(dendextend)
library(WGCNA)
library(GDCRNATools)
library(edgeR)
library(limma)
library(DESeq2)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")
# Load TCGA mRNA data ###################################
counts_gtcga <- readRDS("gtcga_final_counts.RDS") #large numeric df with rows as genes
counts_gtcga <- data.matrix(counts_gtcga) 
dim(counts_gtcga)#19197         596
# Remove Outliers #################################################
# detect outlier genes with gsg
gsg <- goodSamplesGenes(counts_gtcga) #no transpose, wgcna package
summary(gsg) #weather there are outliers
gsg$allOK #false tells you there are outliers
table(gsg$goodGenes) #false is the number of outliers
table(gsg$goodSamples) #no outliers in samples 
# detect outlier samples - hierarchical clustering 
htree <- hclust(dist(t(counts_gtcga)), method = "average") #can take some time 
#save as a big plot
png(file="~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/htree_proteins_gtex_tcga_20250408.png",
    height=3000, width=6000) 
plot(htree) #dstant samples should be excluded
dev.off()
##remove the most outside case!
which(colnames(counts_gtcga)=="TCGA-13-1499-01A-01R-1565-13") #297
counts_gtcga <- counts_gtcga[, -297] #19197   595 left 
dim(counts_gtcga)#19197         595 left
#chek the htree again
dend <- as.dendrogram(htree)
col_aa_red <- ifelse(grepl("GTEX", labels(dend)), "red", "blue")
dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col") #
png(file="~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/htree_protein_gtex_tcga_red20250408.png", height=500, width=3000) #išsaugijimui didesniu formatu
plot(dend2) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()

# Normalize, I will use GDC RNA tools ######################################
#GDC RNA tools, with filter, model added (with GDC RNA tools)
counts <- counts_gtcga
expr = DGEList(counts = counts)
expr = calcNormFactors(expr) #may take some time
snames = colnames(counts_gtcga) #get case names
group = substr(snames, 1, 4); #Sets up level information for samples.

#unfiltered version if needed:
# v_unfiltered <- voom(expr, design= model.matrix(~0 + group), plot = TRUE)$E

#filtered normalisation:
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(counts)
nGenes <- as.numeric(summary(keepALL)[2]) + 
  as.numeric(summary(keepALL)[3])
nKeep <- summary(keepALL)[3]

cat (paste('Total Number of genes: ', nGenes, '\n', sep='')) #19212
cat (paste('Number of genes for downstream analysis: ', nKeep,'\n', sep='')) #13674

exprALL <- expr[keepALL,,keep.lib.sizes = TRUE]
v_filtered <- voom(exprALL, design= model.matrix(~0 + group), plot = TRUE)$E

saveRDS(v_filtered, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/mrna_voom_protein2025.RDS")

# Chek normalization with another clustering plot ################
mRNA_voomt <- t(v_filtered)
htree_norm <- hclust(dist(mRNA_voomt), method = "average") #can take some time
dend3 <- as.dendrogram(htree_norm)
col_aa_red <- ifelse(grepl("GTEX", labels(dend3)), "red", "blue")
dend4 <- assign_values_to_leaves_edgePar(dend=dend3, value = col_aa_red, edgePar = "col") 
png(file="~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/htree_protein_gtex_tcga_norm_red20250408.png", height=500, width=3000) #išsaugijimui didesniu formatu
plot(dend4) #dstant samples should be excluded #save portrait 20x66 inch pdf
dev.off()