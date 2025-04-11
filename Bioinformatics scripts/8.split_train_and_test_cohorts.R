#THIS IS TCGA-OV-RISK-GENES project script No. 8
#split train and test cohorts

# Load packages ##########################################
library(tidyverse)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")
# Load TCGA and GTEX mRNA data, normlised ###################################
gtcga_counts <- readRDS("mrna_voom_protein2025.RDS")
gtcga_counts <- t(gtcga_counts)
gtcga_counts <- as.data.frame(gtcga_counts)


# Set (random-number-generator) seed so that results are consistent between runs
set.seed(18) #choose favorite number
train_ids <- rbinom(nrow(gtcga_counts), size = 1, prob = 0.8) ==1 #choose persentage
#we split clinical data first (probably should be indicated by y instead of x)
gtex_counts_train = gtcga_counts[train_ids, ] 
gtex_counts_test  = gtcga_counts[!train_ids, ] 

dim(gtex_counts_train) #489 13674
dim(gtex_counts_test) #106 13674

snames_train = rownames(gtex_counts_train);
group_train = as.factor(substr(snames_train, 1, 4))
summary(group_train) #153gtex  336tcga

snames_test = rownames(gtex_counts_test);
group_train = as.factor(substr(snames_test, 1, 4))
summary(group_train) #27gtex   79 tcga

#save
saveRDS(gtex_counts_train, "train_gtcga_normcounts_prot_2025.RDS")
saveRDS(gtex_counts_test, "test_gtcga_normcounts_prot_2025.RDS")
