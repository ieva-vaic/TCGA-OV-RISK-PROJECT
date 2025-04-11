#THIS IS TCGA-OV-RISK-GENES project script No.10
#get best genes to predict overall survival in tcga (ovarian cancer)

# Load packages ##########################################
library(glmnet)
library(tidyverse)
library(survival)
library(gplots)
library(survminer)
library(survivalROC)
library(gridExtra)
library(timeROC)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")
# Load train data ###################################
#train data
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
#genes selected by lasso
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
#filter for genes selected by lasso
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 
dim(gtex_filtered_counts_train) #489 samples 214 #genes
#filter for only tcga genes
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_train2) #336 samples #214 genes
# Load clinical data ############################################
pheno <- readRDS("joinedTCGA_XENA_clinical2025.RDS") #full clinical from step3
dim(pheno) #416 samples, 63 variables
#filter clinical for the tcga train cohort cases
train_ids <- colnames(gtex_filtered_counts_train2)
pheno_train <- pheno[pheno$barcode %in% train_ids, ]  #336 samples
rownames(pheno_train) <- pheno_train$barcode
# Create a survival df of survival data #########################
# survival df
clin_df = pheno_train[,
                      c("barcode",
                        "vital_status",
                        "days_to_death",
                        "days_to_last_follow_up",
                        #"neoplasmhistologicgrade",
                        "figo_stage")]

# create a new boolean variable that has TRUE for dead patients and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"
# create an "overall survival" variable -
#that is equal to days_to_death for dead patients -
#and to days_to_last_follow_up for patients who are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
#fix other clincial data:
# For stage remove any of the letters "a", "b" or "c", but only if they are at the end
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[ABC]$", "", clin_df$figo_stage)
table(clin_df$tumor_stage, useNA = "a")
#only one stage 1 case - turn na
clin_df[which(clin_df$tumor_stage == "Stage I"), "tumor_stage"] = NA
table(clin_df$tumor_stage, useNA = "a")
#save survival df
#saveRDS(clin_df, "pheno_survival_only_2025.RDS")
#join clinical ant tcga train counts #####################################
#transform counts:
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))
#add all genes to clin_df
colnames(clin_df)
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
#fix gene names 
colnames(clin_df_joined)

##############################################################################
#COXNET!
#remove nas in survival and deceased status
clin_df_joined2 <- clin_df_joined %>% drop_na(overall_survival) %>% drop_na(deceased) 
#code deceased as 0 and 1 
clin_df_joined2$decesed2 <- as.integer(as.logical(clin_df_joined2$deceased))
#create time and status variables
time <- clin_df_joined2$overall_survival
status <- clin_df_joined2$decesed
#create matrix of overall_survival and deceased status
y2 <- clin_df_joined2[ ,c(7,6)]
#rename matrix
names(y2) <- c("time", "status")
y2 <- as.matrix(y2)
head(y2)
#select only the genes
surv_counts <- clin_df_joined2[, 9:222]
#fit cox
cox_fitx <- glmnet(surv_counts, y2, family="cox")
cox_fitx
plot(cox_fitx) #make a plot of the fits
#get counts matrix
surv_counts <- as.matrix(surv_counts)
set.seed(5)
cvfit <- cv.glmnet(surv_counts, y2, family = "cox", type.measure = "C") #takes some time
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
#get coeficients - gene names
coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10

#SAVE ############################################################
write.csv(res_coef_cox_names, "res_coef_coxnet_names2025.csv")
saveRDS(cox_fitx, "coxnet_fit_2025.RDS")
saveRDS(cvfit, "coxnet_cvfit_2025.RDS")



