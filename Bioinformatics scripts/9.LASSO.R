#THIS IS TCGA-OV-RISK-GENES project script No. 9
#get best genes to separate gtex (normal samples) from tcga (ovarian cancer)

# Load packages ##########################################
library(glmnet)
library(tidyverse)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
gtex_counts_train <- data.matrix(gtex_counts_train)
#perform lasso ########################################
snames = rownames(gtex_counts_train)
group = substr(snames, 1, 4); #Sets up level information for samples.
group = as.factor(group)
#Model selection, lasso (no weak values left)
#using norm.counts
#clinical feature: gtex or TCGA data
set.seed(18)
res_gtex = cv.glmnet(
  x = gtex_counts_train,
  y = group,
  alpha = 0.5,
  family = "binomial") #may take some time
res_gtex #choooses #214
# Getting genes that contribute for the prediction
res_coef_gtex = coef(res_gtex, s="lambda.min") # the "coef" function returns a sparse matrix
head(res_coef_gtex) # in a sparse matrix the "." represents the value of zero
# get coefficients with non-zero values
res_coef_gtex = res_coef_gtex[res_coef_gtex[,1] != 0,] 
res_coef_gtex = res_coef_gtex[-1]
res_coef_gtex_names = names(res_coef_gtex) # get names of the (non-zero) variables.
res_coef_gtex_names  #see the names of the chosen genes
#SAVE ###################################################################
saveRDS(res_gtex, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/elastic_net_model_gtex_2025.RDS")
saveRDS(res_coef_gtex_names, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/gtcga_elastic_2025.RDS")
