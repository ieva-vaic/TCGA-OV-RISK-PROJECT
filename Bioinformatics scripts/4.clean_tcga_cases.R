#THIS IS TCGA-OV-RISK-GENES project script No. 4
#Tidy cases 

# Load packages ##########################################
library(tidyverse)
library(SummarizedExperiment)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")

# Splice clinical data with mRNA counts ################################
#get full TCGA and XENA database clinical data:
pheno_final <- readRDS("pheno.RDS") 
#get full TCGA data
tcga_data <- readRDS("tcga_data.RDS") #quite large!
#get mRNA counts:
tcga_data <- assay(tcga_data) 
# transform mRNA counts:
tcga_counts_t <- t(tcga_data) 
tcga_counts_t <- as.data.frame(tcga_counts_t) #transform to df
tcga_counts_t$barcode <- rownames(tcga_counts_t) #add barcode column for joining
#join tcga mRNA counts with clinical 
tgca_pheno <- right_join(pheno_final, tcga_counts_t, by = "barcode")
dim(tgca_pheno) #429 zmones ir 60723 (60660 genes + 45 clinical)
rownames(tgca_pheno) <- tgca_pheno$barcode #get back rownames
#Remove unwanted cases #################################################
#split unwanted cases from the rest of the data by definition 
table(tgca_pheno$definition, useNA = "a") #7 cases non-primary
split_by_definition <- split(tgca_pheno, f = tgca_pheno$definition, drop = T)
removed_cases <- split_by_definition$`Recurrent Solid Tumor`
mRNA_full <- split_by_definition$`Primary solid Tumor` #new tcga and phenodata
dim(removed_cases) # 7 cases removed
dim(mRNA_full) #422 cases left
#split unwanted cases from the rest of the data by treatment
table(mRNA_full$prior_treatment, useNA = "a") #1 case prior treatment
split_prior_treatment <- split(mRNA_full, f = mRNA_full$prior_treatment, drop = T)
mRNA_full <- split_prior_treatment$No #now the mRNA_full are good cases
removed_cases <- rbind(removed_cases, split_prior_treatment$Yes) 
dim(removed_cases) # 8 cases removed
#split unwanted cases from the rest of the data by icd10 codes
table(mRNA_full$primary_diagnosis, useNA="a")
table(mRNA_full$icd_10_code, useNA="a")
split_non_serous <- split(mRNA_full, f = mRNA_full$icd_10_code, drop = T)
mRNA_full <- split_non_serous$C56.9 #now the mRNA_full are final df of good cases
removed_cases <- rbind(removed_cases, split_non_serous$C48.1)
removed_cases <- rbind(removed_cases, split_non_serous$C48.2) 
dim(removed_cases) # 13  cases removed
dim(mRNA_full) # 416 cases left out of 429, so cheks out
#split back tcga mrna data ##########################
tcga_data <- mRNA_full[, 64:60723] #only mrna counts
pheno_data <-  mRNA_full[, 1:63] #only mrna counts
#saveRDSs ########################################################
saveRDS(mRNA_full, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/tcga_selected_full.RDS") 
saveRDS(tcga_data, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/tcga_selected_counts.RDS") 
#saveRDS(pheno_data, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/tcga_selected_pheno.RDS") 
