##THIS IS TCGA-OV-RISK-GENES project script No. 13
#TEST the risk model

# Load packages ##########################################
library(tidyverse)
library(data.table)
library(glmnet)
library(tidyverse)
library(gplots)
library(survminer)
library(survivalROC)
library(survival)
library(gridExtra)
library(grid)  
library(timeROC)
library(RColorBrewer) 
library(ggprism)
library(rstatix) 
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load test data ###################################
gtex_counts_test <- readRDS("test_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_test <- gtex_counts_test[colnames(gtex_counts_test) %in% gtex_genes] 
#filter for TCGA cases
gtex_filtered_counts_test2 <- gtex_filtered_counts_test %>%
  dplyr::filter(grepl("^TCGA", rownames(gtex_filtered_counts_test)))
dim(gtex_filtered_counts_test2)
# Load clinical data ###################################
clin_df <- readRDS("joinedTCGA_XENA_clinical2025.RDS")
#fix survival
clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
#filter clinical data###################################
test_ids <- rownames(gtex_filtered_counts_test2)
clin_df <- clin_df[clin_df$barcode %in% test_ids, ]  #79 samples
rownames(clin_df) <- clin_df$barcode
#JOIN CLINICAL (MOSTLY SURVIVAL DATA) WITH TRAIN DATA ########################
#add all genes to clin_df
colnames(clin_df)
gtex_filtered_counts_test2$barcode <- rownames(gtex_filtered_counts_test2)
clin_df$barcode == gtex_filtered_counts_test2$barcode 
clin_df_joined_test <- left_join(clin_df, gtex_filtered_counts_test2, by = "barcode")
rownames(clin_df_joined_test) <- clin_df_joined_test$barcode
colnames(clin_df_joined_test)

#GTEX vs TCGA, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
snames = rownames(gtex_counts_test)
group = substr(snames, 1, 4)
group = as.factor(group)
levels(group) <- c("GTEx", "TCGA-OV")
gtex_counts_test2 <- gtex_counts_test
gtex_counts_test2$group <- group
#get genes of interest
expression <- c( "EXO1",   "RAD50",  "PPT2",   "LUC7L2", "PKP3",
                 "CDCA5",  "ZFPL1" , "VPS33B", "GRB7",   "TCEAL4")
#get long df
gtcga_table_full <- reshape2::melt(gtex_counts_test2[, colnames(gtex_counts_test2) %in%
                                                        c("group", expression )],
                                   id.vars="group",
                                   measure.vars= expression)
#get t test
t.test_gtex <- gtcga_table_full %>%
  group_by(variable) %>%
  t_test(value ~ group,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.test_gtex
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.test_gtex_tibble <- t.test_gtex %>% 
  select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(6, 8, 8, 8.5, 10, 
                   8, 6, 6, 10, 10) #choose where to plot p values
  )
t.test_gtex_tibble$p_custom <- ifelse(t.test_gtex_tibble$p < 0.001, 
                                      "p < 0.001", 
                                      paste0("p = ", sprintf("%.3f",
                                                             each.vs.ref_sig$pj)))
#get colors 
custom_colors <- c("GTEx" = "darkgreen","TCGA-OV" = "purple") 
#plot
gtex_plot <- ggplot(gtcga_table_full, aes(x=group , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = group )) +
  geom_jitter(aes(color = group ), size=1, alpha=0.5) +
  ylab(label = expression("Normalised expression")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_gtex_tibble, label = "p_custom") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) 

gtex_plot
