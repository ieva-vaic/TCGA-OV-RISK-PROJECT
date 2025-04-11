#THIS IS TCGA-OV-RISK-GENES project script No. 12
#Make risk score, compare to clinical features

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
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 
#filter for TCGA cases
gtex_filtered_counts_train2 <- gtex_filtered_counts_train %>%
  dplyr::filter(grepl("^TCGA", rownames(gtex_filtered_counts_train)))
dim(gtex_filtered_counts_train2)
#gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))
# Load clinical data ###################################
clin_df <- readRDS("joinedTCGA_XENA_clinical2025.RDS")
#fix survival
clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
#filter clinical data###################################
train_ids <- rownames(gtex_filtered_counts_train2)
clin_df <- clin_df[clin_df$barcode %in% train_ids, ]  #336 samples
rownames(clin_df) <- clin_df$barcode
#JOIN CLINICAL (MOSTLY SURVIVAL DATA) WITH TRAIN DATA ########################
#add all genes to clin_df
colnames(clin_df)
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df$barcode == gtex_filtered_counts_train2$barcode 
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
rownames(clin_df_joined) <- clin_df_joined$barcode
colnames(clin_df_joined)
#saveRDS(clin_df_joined, "clin_df_joined_2025.RDS")
# Load cox model ###################################
cvfit <- readRDS("coxnet_cvfit_2025.RDS")
cox_fitx <- readRDS("coxnet_fit_2025.RDS")

#get coeficients - gene names
coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10

#VALUE OF COEFICIENTS PLOT ##################################
cv_model <- cvfit
# Extract coefficients at optimal lambda
coefs <- coef_x
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df <- rownames_to_column(coefs_df, var = "Feature")
colnames(coefs_df)[2] <- "Coefficient"
print(coefs_df)
# Plot coeficients ###########################
# get colors
n_features <- nrow(coefs_df)
colors <- colorRampPalette(brewer.pal(9, "Set3"))(n_features)  # Light pastel rainbow colors
# Pretty ggplot2 
ggplot(coefs_df, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Feature)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +  # Colorful bars
  scale_fill_manual(values = colors) +  # Apply colors to bars
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs( x = "Feature", y = "Value of coefficients") +
  theme(
    axis.text.y = element_text(face = "italic"),  # Italic feature names
    panel.grid.major.y = element_blank(),  # Remove gridlines for cleaner look
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Hide legend for a clean visual
  )

#RISK SCORE ##################################
# - cv_model: Your trained cv.glmnet model
cv_model
# - gene_data: A data frame (or matrix) with the expression values for your genes 
# (rows = samples, columns = genes)
# It should have the same feature names as used in the model (genes).
gene_data <- clin_df_joined[, res_coef_cox_names]
# - gene_list: A vector with the list of genes you're interested in.
res_coef_cox_names
# coefs
coefs_df
# Calculate the risk score: linear combination of gene expressions and coefficients
# Risk score = sum( gene_expression * coefficient )
risk_scores <- rowSums(sweep(gene_data, 2, coefs_df$Coefficient, "*"))
# View the risk scores
print(risk_scores) # now I have some risk scores

#add risk scores to the clin_df_joined
clin_df_joined$RiskScore <- risk_scores[rownames(clin_df_joined)]

#survival time vs coeficient plot ##################################
clin_df_joined$overall_survival
clin_df_joined$RiskScore
clin_df_joined$vital_status

# Create the dot plot
ggplot(clin_df_joined, aes(x = overall_survival, y = RiskScore, color = vital_status)) +
  geom_point(size = 4, alpha = 0.7) +  # Dot plot with a slight transparency
  scale_color_manual(values = c("turquoise", "deeppink")) +  # Blue for Alive, Red for Deceased
  labs(title = "Risk Score vs. Overall Survival",
       x = "Overall Survival",
       y = "Risk Score") +
  geom_hline(yintercept = -0.05914624, linetype = "dotted", color = "black", size = 1) +  # Dotted line at median (calculated below)
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())  

#KM plot ##################################
# Calculate the median risk score
median_risk <- median(clin_df_joined$RiskScore, na.rm = TRUE) #-0.05914624
# Create a new factor column based on the median value
clin_df_joined$RiskGroup <- ifelse(clin_df_joined$RiskScore <= median_risk, "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = clin_df_joined$overall_survival, event = clin_df_joined$deceased )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = clin_df_joined)
# Plot the Kaplan-Meier curve using ggsurvplot
ggsurvplot(km_fit, data = clin_df_joined, 
           pval = TRUE,  # Show p-value of the log-rank test
           risk.table = TRUE,  # Add risk table below the plot
           title = "Kaplan-Meier Plot: Low vs. High Risk based on Risk Score",
           xlab = "Overall Survival Time",
           ylab = "Survival Probability",
           palette = c("turquoise", "deeppink"),  # Color palette for groups
           legend.title = "Risk Group", 
           legend.labs = c("Low Risk", "High Risk"))

#SURVIVAL ROCs for separate markers ##################################
surv_df <- clin_df_joined[, colnames(clin_df_joined) %in%
                          c("deceased", "overall_survival", res_coef_cox_names)]
surv_df <- surv_df %>% dplyr::rename(censor = deceased, surv_time = overall_survival) 
#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365
#make surf df but only of my genes!
coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df) #334  10 #kazkur pametu 2 zmones, kol kas neieskosiu - turbut tie kurie be deceased
time <- surv_df$surv_time
event <- surv_df$censor
rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
plot_list <- list()
# Loop through rez_list and create individual plots
for (i in seq_along(rez_list)) {
  p <- ggplot() +
    geom_line(aes(x = rez_list[[i]]$FP, y = rez_list[[i]]$TP)) +
    xlim(0, 1) + ylim(0, 1) +
    xlab(paste("FP", "\n", "AUC = ", round(rez_list[[i]]$AUC, 3))) +
    ylab("TP") +
    ggtitle(paste(names(rez_list)[i], ", Method = KM, Year = 1")) +
    theme_minimal() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # Add diagonal line
  
  # Store the plot in the plot_list
  plot_list[[i]] <- p
}
# Now, arrange all the plots in a grid and save the combined plot
combined_plot <- grid.arrange(grobs = plot_list, ncol = 2) 
ggsave("~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/SURVROC_SEPARATE20250410.png",
       combined_plot, width = 10, height = 15, dpi = 300)
#SURVIVAL ROC for my risk score ##################################
surv_df2 <- clin_df_joined[, colnames(clin_df_joined) %in%
                            c("deceased", "overall_survival", "RiskScore")]
surv_df2 <- surv_df2 %>% dplyr::rename(censor = deceased, surv_time = overall_survival) 
#need these for survrock
nobs <- NROW(surv_df2)
cutoff <- 365
#make surf df but only of my genes!
time2 <- surv_df2$surv_time
event2 <- surv_df2$censor
roc_result <- survivalROC(Stime = time2, 
                          status = event2, 
                          marker = surv_df2[, "RiskScore"], 
                          predict.time = cutoff, 
                          method = "KM")
# Plot the ROC curve for this single marker
p <- ggplot() +
  geom_line(aes(x = roc_result$FP, y = roc_result$TP)) +
  xlim(0, 1) + ylim(0, 1) +
  xlab(paste("FP", "\n", "AUC = ", round(roc_result$AUC, 3))) +
  ylab("TP") +
  ggtitle(paste("Marker: ", "RiskScore", ", Method = KM, Year = 1")) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # Add diagonal line

# Print the plot
print(p)

#TIME ROC for the risk score ###############################################
t_eval <- c(365, 1095, 1825)  # time points
roc_result <- timeROC(
  T = time2,       # Survival time from df
  delta = event2, # Event indicator from df
  marker = surv_df2[, "RiskScore"], # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)

roc_result #time roc does basically the same as survival roc, but give se 
plot(roc_result, time = 1825) 
plot(roc_result, time = 1095) 
plot(roc_result, time = 365) 

#With stage risk model ###########################
#choose variables
colnames(clin_df_joined)
table(clin_df_joined$clinicalstage2, useNA = "a")
clin_df_joined$RiskScore

#get anova
kruskal_test <- kruskal.test(clin_df_joined$RiskScore,
                             clin_df_joined$clinicalstage2, 
                             var.equal = F,
                             #alternative = "two.sided", 
                             na.rm = TRUE)
kruskal.testp_value <- kruskal_test$p.value
kruskal.testp_value
#get colors 
custom_colors <- c("Stage I" = "pink2", "Stage II" = "lightpink","Stage III" = "deeppink",
                   "Stage IV" = "darkviolet") 
#plot
stage_plot <- ggplot(clin_df_joined, aes(x=clinicalstage2 , y=RiskScore, fill = clinicalstage2)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = clinicalstage2 )) +
  geom_jitter(aes(color = clinicalstage2 ), size=1, alpha=0.5) +
  ylab(label = expression("Risk score")) + 
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
  scale_color_manual(values = custom_colors) +
  ggtitle("Risk score correlation with stage ")+
  annotate("text", x = 2, y = max(clin_df_joined$RiskScore) + 0.11, 
           label = paste("p = ", format(kruskal.testp_value, digits = 3)), 
           size = 5, color = "black")

stage_plot

#With grade risk model###########################
table(clin_df_joined$neoplasmhistologicgrade, useNA = "a")
clin_df_joined$neoplasmhistologicgrade <- recode(clin_df_joined$neoplasmhistologicgrade,
                                                  "GB" = NA_character_,
                                                  "GX" = NA_character_,
                                                  "G4" = NA_character_)
# Remove NA values from both columns using complete.cases()
clin_df_joined_gr <- clin_df_joined[complete.cases(clin_df_joined$RiskScore,
                                                   clin_df_joined$neoplasmhistologicgrade), ]
clin_df_joined_gr$neoplasmhistologicgrade <- as.factor(clin_df_joined_gr$neoplasmhistologicgrade)

ttest_grade <- t.test(clin_df_joined_gr$RiskScore ~ clin_df_joined_gr$neoplasmhistologicgrade ,  
                      var.equal = F,
                      alternative = "two.sided")
ttest_grade_p <- ttest_grade$p.value

#get colors 
custom_colors <- c("G2" = "lightpink","G3" = "deeppink") 
#plot
grade_plot <- ggplot(clin_df_joined, aes(x=neoplasmhistologicgrade ,
                                          y=RiskScore, fill = neoplasmhistologicgrade)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = neoplasmhistologicgrade )) +
  geom_jitter(aes(color = neoplasmhistologicgrade ), size=1, alpha=0.5) +
  ylab(label = expression("Risk score")) + 
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
  scale_color_manual(values = custom_colors) +
  ggtitle("Risk score correlation with grade ")+
  annotate("text", x = 2, y = max(clin_df_joined$RiskScore) + 0.11, 
           label = paste("p = ", format(ttest_grade_p, digits = 3)), 
           size = 5, color = "black")

grade_plot
