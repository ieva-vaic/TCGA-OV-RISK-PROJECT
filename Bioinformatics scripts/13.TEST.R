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

# Load cox model ###################################
cvfit <- readRDS("coxnet_cvfit_2025.RDS")
cox_fitx <- readRDS("coxnet_fit_2025.RDS")

#get coeficients - gene names
coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10


#RISK SCORE ##################################
# - cv_model: Your trained cv.glmnet model
cv_model <- cvfit
# - gene_data: A data frame (or matrix) with the expression values for your genes 
# (rows = samples, columns = genes)
# It should have the same feature names as used in the model (genes).
gene_data_test <- clin_df_joined_test[, res_coef_cox_names]
# - gene_list: A vector with the list of genes you're interested in.
res_coef_cox_names
# coefs
# Extract coefficients at optimal lambda
coefs <- coef_x
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df <- rownames_to_column(coefs_df, var = "Feature")
colnames(coefs_df)[2] <- "Coefficient"
print(coefs_df)
# Calculate the risk score: linear combination of gene expressions and coefficients
# Risk score = sum( gene_expression * coefficient )
risk_scores_test <- rowSums(sweep(gene_data_test, 2, coefs_df$Coefficient, "*"))
# View the risk scores
print(risk_scores_test) # now I have some risk scores

#add risk scores to the clin_df_joined_test
clin_df_joined_test$RiskScore <- risk_scores_test[rownames(clin_df_joined_test)]
#create df wih survival data
surv_df_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in%
        c("deceased", "overall_survival", res_coef_cox_names, "RiskScore")]
surv_df_test <- surv_df_test %>%
  dplyr::rename(censor = deceased, surv_time = overall_survival) 
#create features for timeroc
nobs <- NROW(surv_df_test)
#make surf df but only of my genes!
time <- surv_df_test$surv_time
event <- surv_df_test$censor
#time roc
t_eval <- c(365, 1095, 1825)  # time points
roc_result <- timeROC(
  T = time,       # Survival time from df
  delta = event, # Event indicator from df
  marker = surv_df_test[, "RiskScore"], # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)

#time rocs for separate biomarkers ######################################

coxnet.df <- surv_df_test[, (colnames(surv_df_test) %in% res_coef_cox_names)]
dim(coxnet.df)

rez_list <- apply(coxnet.df, 2, timeROC,
                   T = time,       # Survival time from df
                   delta = event, # Event indicator from df
                   #marker  # Predictor already in the df
                   cause = 1,         # Event of interest
                   times = t_eval,    # Time points for ROC
                   iid = TRUE )        # Compute confidence intervals)

auc_table <- map_dfr(names(rez_list), function(gene) {
  roc <- rez_list[[gene]]
  
  tibble(
    gene = gene,
    time = roc$times,
    cases = roc$cases,
    survivors = roc$survivors,
    censored = roc$censored,
    auc = roc$AUC,
    se = roc$inference$vect_sd_1
  )
})

auc_table
write.csv( auc_table, "~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/AUC.TABLE_test.csv")

##plot at year 1###############
# Choose target time
target_time <- 365
time_index <- which(rez_list[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity (FPR)",
  ylab = "Sensitivity (TPR)",
  main = paste("Time-dependent ROC Curves at", target_time, "days"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)

# Add legend (last color is black for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)
##plot at year 3###############
# Choose target time
target_time <- 1095    
time_index <- which(rez_list[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity (FPR)",
  ylab = "Sensitivity (TPR)",
  main = paste("Time-dependent ROC Curves at", target_time, "days"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)

# Add legend (last color is black for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

##plot at year 5 ###############
# Choose target time
target_time <- 1825        
time_index <- which(rez_list[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity (FPR)",
  ylab = "Sensitivity (TPR)",
  main = paste("Time-dependent ROC Curves at", target_time, "days"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)

# Add legend (last color is black for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

#KM plot with RISK SCORE#################

# Calculate the median risk score
median_risk <- median(clin_df_joined_test$RiskScore, na.rm = TRUE) #--0.07747393
# Create a new factor column based on the median value
clin_df_joined_test$RiskGroup <- ifelse(clin_df_joined_test$RiskScore <= median_risk,
                                   "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = clin_df_joined_test$overall_survival,
                    event = clin_df_joined_test$deceased )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = clin_df_joined_test)
# Plot the Kaplan-Meier curve using ggsurvplot
ggsurvplot(km_fit, data = clin_df_joined_test, 
           pval = TRUE,  # Show p-value of the log-rank test
           risk.table = TRUE,  # Add risk table below the plot
           title = "Kaplan-Meier Plot: Low vs. High Risk based on Risk Score in test cohort",
           xlab = "Overall Survival Time",
           ylab = "Survival Probability",
           palette = c("turquoise", "deeppink"),  # Color palette for groups
           legend.title = "Risk Group", 
           legend.labs = c("Low Risk", "High Risk"))


