# Load packages ##########################################
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer) 
library(circlize)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load test data ###################################
gtex_counts_test <- readRDS("test_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_test <- gtex_counts_test[colnames(gtex_counts_test) %in% gtex_genes] 
#get TCGA df
gtex_filtered_counts_test2 <- as.data.frame(t(gtex_filtered_counts_test))
gtex_filtered_counts_test2 <- gtex_filtered_counts_test2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_test2)
#336 samples #214 genes - leaves only the genes left at the lasso step
#get GTEX df
gtex_filtered_counts_test3 <- as.data.frame(t(gtex_filtered_counts_test))
gtex_filtered_counts_test3 <- gtex_filtered_counts_test3 %>%  dplyr::select(starts_with("GTEX")) 
dim(gtex_filtered_counts_test3)

#READ FULL CLINICAL DATA#######################################
pheno_full <- readRDS("joinedTCGA_XENA_clinical.RDS")
#reduce to train data
tcga_cases <- colnames(gtex_filtered_counts_test2)
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
pheno_best <- bind_rows(pheno_best, tibble(barcode = colnames(gtex_filtered_counts_test3)))
rownames(pheno_best) <- pheno_best$barcode
pheno_best <- pheno_best %>%
  arrange(match(barcode, rownames(gtex_filtered_counts_test)))
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
gtex_counts_test <- data.matrix(gtex_filtered_counts_test)
#gene names:
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5",
                "ZFPL1","VPS33B", "GRB7","TCEAL4")
expression %in% colnames(gtex_counts_test) #check if they are in the data
#chosen genes highlated in red
col_colors <- sapply(colnames(gtex_counts_test), function(x) {
  if (x %in% expression) {
    "red"  # Color selected columns red
  } else {
    "black"  # Default color
  }
})
#CREATE GROUPINGS ACCORDING TO DATA##############################
snames = rownames(gtex_counts_test)
group = substr(snames, 1, 4)
group = as.factor(group)
levels(group) <- c("GTEx", "TCGA-OV")

#GENERATE HEATMAP##############################################
# Generate heatmap with colored column names
heatmap_tcga <- Heatmap(as.matrix(gtex_filtered_counts_test), 
                        row_split = group,   
                        show_row_names = F,
                        column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                        row_names_gp = gpar(fontsize = 2), # 
                        heatmap_legend_param = list(title = "Gene Expression"),
                        right_annotation = row_ha,
                        cluster_rows = F)

#SAVE HEATMAP###################################################################
png("~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics//gtextcga_heatmap_test2025-04-18.png",
    width = 6500, height = 4000, res = 400) # width and height in pixels, resolution in dpi
heatmap_tcga #
dev.off() # Close the PNG device
#CREATE SMALLER HEATMAP OF CHOSEN GENES#######################################
small_test_df <- gtex_filtered_counts_test[colnames(gtex_filtered_counts_test) %in% expression]
small_test_matrix <- data.matrix(small_test_df)
heatmap_tcga2 <- Heatmap(as.matrix(small_test_matrix), 
                         row_split = group,   
                         show_row_names = F,
                         column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Gene Expression"),
                         right_annotation = row_ha,
                         cluster_rows = F)
#SAVE SMALL HEATMAP###################################################################
png("~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics//gtextcga_heatmap_test_small2025-04-18.png",
    width = 2000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
heatmap_tcga2 #
dev.off() # Close the PNG device

#add GTEx AGE###########################################################
# Load gtex sample data ###################################
#dowloaded https://gtexportal.org/home/downloads/adult-gtex/metadata on 2025-04-18
metadata <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Check the first few rows
head(metadata)

# Optional: view column names: id, sex (2= Female),
#DTHHRDY = hardy sclase of death (1 = trauma (sudden), 2 = suden death (<1h) e.g. miocardal infactrion, 
#3 = patients that vere ill but death was unexpected (>24h),
#4 = long illness such as cancer, 0 = ventilator cases)
colnames(metadata)

#filter for test data#################
#ids
ids <- colnames(gtex_filtered_counts_test3)
metadata$SUBJID
#trim ids in gtex_df_train
# Keep only the first two segments
donor_ids <- sub("^([^-]+-[^-]+).*", "\\1", ids)
# Result
donor_ids %in% metadata$SUBJID #all represented in metadata
#filter metadata
filtered_metadata <- metadata[metadata$SUBJID %in% donor_ids, ]

#look at age
filtered_metadata$AGE
pheno_best$AGE <- cut(pheno_best$age,
                      breaks = seq(30, 80, by = 10),   # From 30 to 80 in 10-year steps
                      right = FALSE,                  # e.g., 60 is in "60-69"
                      labels = c("30-39", "40-49", "50-59", "60-69", "70-79"))

pheno_best$AGE

#fix pheno_best
pheno_best$SUBJID <- sub("^([^-]+-[^-]+).*", "\\1", pheno_best$barcode)
#merge
pheno_best$AGE2 <- filtered_metadata$AGE[match(pheno_best$SUBJID, filtered_metadata$SUBJID)]
pheno_best$AGE2 <- factor(pheno_best$AGE2)
#merge
pheno_best$AGE2 <- ifelse(!is.na(pheno_best$AGE),
                          as.character(pheno_best$AGE),
                          as.character(pheno_best$AGE2))

# Convert back to factor with consistent levels (optional)
pheno_best$AGE2 <- factor(pheno_best$AGE2, levels = levels(pheno_best$AGE))
pheno_best$AGE2
#CREATE HEATMAP CLINICAL DATA FEATURES with GTEX AGE######################
row_ha2 = rowAnnotation(Race  = pheno_best$race, 
                       Age = pheno_best$AGE2, 
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
                                  Age = c("30-39" = "#E6D6F5",  
                                          "40-49" = "#C49DDE", 
                                          "50-59" = "#9C6CD3",  
                                          "60-69" = "#713AB6", 
                                          "70-79" = "#4B1C74", 
                                          "NA" = "grey"),
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
#MAKE SMALL HEATMAP WITH GETX AGE################
heatmap_tcga3 <- Heatmap(as.matrix(small_test_matrix), 
                         row_split = group,   
                         show_row_names = F,
                         column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Gene Expression"),
                         right_annotation = row_ha2,
                         cluster_rows = F)
#SAVE SMALL HEATMAP2###################################################################
png("~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics//gtextcga_heatmap_GTEXAGE_test_small2025-04-18.png",
    width = 2000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
heatmap_tcga3 #
dev.off() # Close the PNG device