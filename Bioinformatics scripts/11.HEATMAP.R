#THIS IS TCGA-OV-RISK-GENES project script No. 11
#HEATMAP 
# Load packages ##########################################
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer) 
library(circlize)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
gtex_genes <- readRDS("gtcga_elastic_2025.RDS") #only the genes left after lasso
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 
dim(gtex_filtered_counts_train) #489 samples
#get TCGA df
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_train2)
#336 samples #214 genes - leaves only the genes left at the lasso step
#get GTEX df
gtex_filtered_counts_train3 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train3 <- gtex_filtered_counts_train3 %>%  dplyr::select(starts_with("GTEX")) 
dim(gtex_filtered_counts_train3)
#336 samples #214 genes - leaves only the genes left at the lasso step
#READ FULL CLINICAL DATA#######################################
pheno_full <- readRDS("joinedTCGA_XENA_clinical.RDS")
#reduce to train data
tcga_cases <- colnames(gtex_filtered_counts_train2)
sum(pheno_full$barcode %in% tcga_cases) #336
pheno_full <- pheno_full[pheno_full$barcode %in% tcga_cases, ]
dim(pheno_full)

#LEAVE ONLY FEATURES OF INTEREST IN CLINICAL DATA##############################
pheno_best <- pheno_full[, (names(pheno_full)
                            %in% c("barcode", "STAGE", "race", "vital_status", 
                                   "ageatinitialpathologicdiagnosis", "neoplasmhistologicgrade", 
                                   "lymphaticinvasion", "treatment_type"  ))] # left
#fix up grade by removing GB, GX, and singular Grade 4 case
pheno_best$neoplasmhistologicgrade <- recode(pheno_best$neoplasmhistologicgrade,
                                             "GB" = NA_character_,
                                             "GX" = NA_character_,
                                             "G4" = NA_character_)
table(pheno_best$neoplasmhistologicgrade)
#rename clinical features to more conventional names
pheno_best <- pheno_best %>%
  rename(grade = neoplasmhistologicgrade,
         stage = STAGE,
         `vital status` = vital_status,
         `treatment type` = treatment_type,
         age = ageatinitialpathologicdiagnosis,
         `lymphatic invasion` = lymphaticinvasion)
colnames(pheno_best)
#fix up lymphatic invasion
pheno_best$`lymphatic invasion` <- factor(pheno_best$`lymphatic invasion`, levels = c(0 , 1 , NA ),
                                          labels = c("No invasion",
                                                     "Lymphatic invasion"
                                          ))
pheno_best$`lymphatic invasion` <- as.character(pheno_best$`lymphatic invasion`)
#add gtex to clinical data as NAs
pheno_best <- bind_rows(pheno_best, tibble(barcode = colnames(gtex_filtered_counts_train3)))
rownames(pheno_best) <- pheno_best$barcode
pheno_best <- pheno_best %>%
  arrange(match(barcode, rownames(gtex_filtered_counts_train)))
pheno_best <- pheno_best %>%
  mutate(across(everything(), ~replace(., is.na(.), "NA")))
#fix the numeric variable
pheno_best$age <- as.numeric(pheno_best$age)  

#CREATE HEATMAP CLINICAL DATA FEATURES ######################
col_age <- colorRamp2(c(30, 100), c( "#9cd4c4", "#3c402f"))
row_ha = rowAnnotation(Race  = pheno_best$race, 
                       Age = pheno_best$age, 
                       `Vital status` = pheno_best$`vital status`,
                       `Treatment type` = pheno_best$`treatment type`,
                       Grade = pheno_best$grade,
                       Stage =pheno_best$stage,
                       `Lymphatic invasion` = pheno_best$`lymphatic invasion`,
                       #choose colors
                       col = list(Race = c("american indian or alaska native" = "pink",
                                           "asian" = "#9cd4c4",
                                           "black or african american" = "darkgreen",
                                           "native hawaiian or other pacific islander" = "turquoise",
                                           "white" = "lightgrey",
                                           "not reported" = "darkgrey", 
                                           "NA" = "grey"), 
                                  Age = col_age,
                                  `Vital status` = c("Alive" = "#9cd4c4",  
                                                     "Dead" = "#a89cd4", 
                                                     "NA" = "grey"),
                                  `Treatment type` = c("Pharmaceutical Therapy, NOS" = "#9cd4c4",  
                                                       "Radiation Therapy, NOS" = "#a89cd4", 
                                                       "NA" = "grey"),
                                  Grade = c("G2" = "#9cd4c4",  
                                            "G3" = "#a89cd4", 
                                            "NA" = "grey"),
                                  Stage = c("Stage I" = "#9cd4c4",  
                                            "Stage II" = "#a89cd4",
                                            "Stage III" = "pink",  
                                            "Stage IV" = "turquoise", 
                                            "NA" = "grey"),
                                  `Lymphatic invasion` = c("No invasion" = "#9cd4c4",  
                                                           "Lymphatic invasion" = "#a89cd4", 
                                                           "NA" = "grey")
                       ))

#HIGHLITE THE CHOSEN GENES####################################
#gene names:
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5",
                "ZFPL1","VPS33B", "GRB7","TCEAL4")
expression %in% colnames(gtex_counts_train) #check if they are in the data
#chosen genes highlated in red
col_colors <- sapply(colnames(gtex_counts_train), function(x) {
  if (x %in% expression) {
    "red"  # Color selected columns red
  } else {
    "black"  # Default color
  }
})
#CREATE GROUPINGS ACCORDING TO DATA##############################
snames = rownames(gtex_counts_train)
group = substr(snames, 1, 4)
group = as.factor(group)
levels(group) <- c("GTEx", "TCGA-OV")


#GENERATE HEATMAP##############################################
# Generate heatmap with colored column names
gtex_counts_train <- data.matrix(gtex_filtered_counts_train)
heatmap_tcga <- Heatmap(as.matrix(gtex_filtered_counts_train), 
                        row_split = group,   
                        show_row_names = F,
                        column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                        row_names_gp = gpar(fontsize = 2), # 
                        heatmap_legend_param = list(title = "Gene Expression"),
                        right_annotation = row_ha,
                        cluster_rows = F)

#SAVE HEATMAP###################################################################
png("~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics//gtextcga_heatmap2025-04-010.png",
    width = 6500, height = 4000, res = 400) # width and height in pixels, resolution in dpi
heatmap_tcga #
dev.off() # Close the PNG device

#CREATE SMALLER HEATMAP OF CHOSEN GENES#######################################
small_train_df <- gtex_filtered_counts_train[colnames(gtex_filtered_counts_train) %in% expression]
small_train_matrix <- data.matrix(small_train_df)
heatmap_tcga2 <- Heatmap(as.matrix(small_train_matrix), 
                         row_split = group,   
                         show_row_names = F,
                         column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Gene Expression"),
                         right_annotation = row_ha,
                         cluster_rows = F)
#SAVE SMALL HEATMAP###################################################################
png("~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics//gtextcga_heatmap_small2025-04-010.png",
    width = 2000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
heatmap_tcga2 #
dev.off() # Close the PNG device
