###THIS IS TCGA-OV-RISK-GENES project script No. 16
#TEST the model on GSE data

# Load packages#######################################
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer) 
library(circlize)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/GEO/")
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(annotate)
library(hgu133plus2.db)
library(rstatix) 
library(ggprism)
library(car)
library(survminer)
library(survivalROC)
library(timeROC)
library(survival)

###########redo######################################################################
library(curatedOvarianData)
data(package = "curatedOvarianData")  # list datasets
#get GSE9891
data("GSE9891_eset")  # load the ExpressionSet for GSE9891
#get clinical df GSE9891
clinical_data <- pData(GSE9891_eset)
head(clinical_data[, c("days_to_death", "vital_status")])  
#look at clincial data
table(clinical_data$tumorstage, useNA = "a") 
table(clinical_data$summarystage, useNA = "a") 
table(clinical_data$summarygrade, useNA = "a") 
table(clinical_data$grade, useNA = "a") #
table(clinical_data$sample_type, useNA = "a") 
#get expression matrix
exprs_matrix <- exprs(GSE9891_eset)
#find out if my chosen genes are in the matrix
#my 10 gene matrix
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5",
                "ZFPL1","VPS33B", "GRB7","TCEAL4")
expression %in% rownames(exprs_matrix) #no LUC7L2
#fix luc by aggregating the three transciprts
luc7l2 <- exprs_matrix[grep("*LUC7L2", rownames(exprs_matrix)), ]
luc7l2_names <- rownames(luc7l2)
luc_mean <- colMeans(exprs_matrix[luc7l2_names, , drop = FALSE])
# Add the mean as a new row named "luc"
exprs_matrix <- rbind(exprs_matrix, LUC7L2 = luc_mean)
expression %in% rownames(exprs_matrix) #no LUC7L2
#214 genes
gtex_genes <- readRDS("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/gtcga_elastic_2025.RDS")
gtex_genes %in%  rownames(exprs_matrix) #NOT ALL in there

#make expression df of interest genes #########################
small_expr_df <- exprs_matrix[rownames(exprs_matrix) %in% expression,]
#transform df
small_expr_df_t <- t(small_expr_df)
small_expr_df_t <- data.frame(small_expr_df_t)

#make heatmap#######################################
#create clinical annotation
row_ha2 = rowAnnotation(Type = clinical_data$sample_type,
                        Grade = clinical_data$grade,
                        Stage =clinical_data$tumorstage,
                        #choose colors
                        col = list(
                          `Type` = c("borderline" = "#9cd4c4",  
                                     "tumor" = "#a89cd4"),
                          Grade = c("1" = "#9cd4c4",  
                                    "2" = "turquoise", 
                                    "3" = "#a89cd4",
                                    "NA" = "grey"),
                          Stage = c("1" = "#FFB6C1", 
                                    "2" = "#E6E6FA",
                                    "3" = "#F08080", 
                                    "4" = "#9cd4c4",  
                                    "NA" = "grey")
                        ))
#MAKE SMALL HEATMAP WITH GETX AGE################
col_fun <- colorRamp2(c(1, 7, 13), c("#013220", "white", "#8B0000"))
heatmap_geo <- Heatmap(as.matrix(t(small_expr_df)) ,  
                       show_row_names = F,
                       row_split = clinical_data$sample_type, 
                       column_names_gp = gpar(fontsize = 6, fontface = "italic"), 
                       row_names_gp = gpar(fontsize = 2), # 
                       heatmap_legend_param = list(title = "Gene Expression"),
                       right_annotation = row_ha2,
                       cluster_rows = F, col = col_fun)
heatmap_geo

# Load cox model ###################################
cvfit <- readRDS("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/coxnet_cvfit_2025.RDS")
cox_fitx <- readRDS("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/coxnet_fit_2025.RDS")

#get coeficients - gene names
coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10
# - cv_model: Your trained cv.glmnet model
cv_model <- cvfit
# Extract coefficients at optimal lambda
coefs <- coef_x
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df <- rownames_to_column(coefs_df, var = "Feature")
colnames(coefs_df)[2] <- "Coefficient"
print(coefs_df)
# Risk score = sum( gene_expression * coefficient )
risk_scores_test <- rowSums(sweep(small_expr_df_t, 2, coefs_df$Coefficient, "*"))
# View the risk scores
print(risk_scores_test) # now I have some risk scores
#add risk scores to df
# Match by rownames and assign
small_expr_df_t$risk_scores <- unlist(risk_scores_test)[match(rownames(small_expr_df_t),
                                                              names(risk_scores_test))]


#create one df wih survival data #######################
#add clinical data
leave_clincial <- c("sample_type", "histological_type",
                    "tumorstage", "grade", "age_at_initial_pathologic_diagnosis",
                    "recurrence_status","days_to_tumor_recurrence", "days_to_death","vital_status") #choose what to leave
small_clinical <- clinical_data[ colnames(clinical_data) %in% leave_clincial]
str(small_clinical)
str(small_expr_df_t)

# Merge on 'ID' column
small_geo_df <- merge(small_clinical, small_expr_df_t, by = "row.names")
rownames(small_geo_df) <- small_geo_df$Row.names

#make only surv df
surv_df_geo <- small_geo_df[, colnames(small_geo_df) %in%
                              c("vital_status", "days_to_death", expression, "risk_scores")]
surv_df_geo <- surv_df_geo %>%
  dplyr::rename(censor = vital_status, surv_time = days_to_death) 
#create features for timeroc
nobs <- NROW(surv_df_geo)
#make surf df but only of my genes!
time <- surv_df_geo$surv_time
event <- surv_df_geo$censor

event <- ifelse(event == "deceased", 1, 0)

small_geo_df$recurrence_status_f <- ifelse(small_geo_df$recurrence_status == "recurrence", 1, 0)
small_geo_df$vital_status_f <- ifelse(small_geo_df$vital_status == "deceased", 1, 0)
#time roc for risk score########################
t_eval <- c(365, 1095, 1825, 3650)  # time points
roc_result <- timeROC(
  T = time,       # Survival time from df
  delta = event, # Event indicator from df
  marker = surv_df_geo[, "risk_scores"], # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result #for risk score

#time rocs for separate biomarkers ######################################
exprs.df <- small_geo_df[, (colnames(small_geo_df) %in% expression)]
dim(exprs.df)

rez_list <- apply(exprs.df, 2, timeROC,
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
write.csv(auc_table, "~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/geo9891_auc_table.csv")

##plot at year 5 ###############
par(pty = "s")
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
  main = paste("Time-dependent ROC Curves at", target_time, "days, in GSE9891"),
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
# Add risk score ROC line in bold
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)
# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")
# legend names
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)
# Add legend 
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
median_risk <- median(small_geo_df$risk_scores, na.rm = TRUE) #--0.7791149
# Create a new factor column based on the median value
small_geo_df$RiskGroup <- ifelse(small_geo_df$risk_scores <= median_risk,
                                 "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time ,
                    event  )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = small_geo_df)
# Plot the Kaplan-Meier curve using ggsurvplot
ggsurvplot(km_fit, data = small_geo_df, 
           pval = TRUE,  # Show p-value of the log-rank test
           risk.table = TRUE,  # Add risk table below the plot
           title = "Kaplan-Meier Plot: Low vs. High Risk based on Risk Score in gse9891 cohort",
           xlab = "Overall Survival Time",
           ylab = "Survival Probability",
           palette = c("turquoise", "deeppink"),  # Color palette for groups
           legend.title = "Risk Group", 
           legend.labs = c("Low Risk", "High Risk"))

#KM plot, PFS, with RISK SCORE#################
# Calculate the median risk score
median_risk <- median(small_geo_df$risk_scores, na.rm = TRUE) #-0.7791149
# Create a new factor column based on the median value
small_geo_df$RiskGroup <- ifelse(small_geo_df$risk_scores <= median_risk,
                                 "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(small_geo_df$days_to_tumor_recurrence ,
                    small_geo_df$recurrence_status_f  )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = small_geo_df)
# Plot the Kaplan-Meier curve using ggsurvplot
ggsurvplot(km_fit, data = small_geo_df, 
           pval = TRUE,  # Show p-value of the log-rank test
           risk.table = TRUE,  # Add risk table below the plot
           title = "Kaplan-Meier Plot: Low vs. High Risk based on Risk Score in gse9891 cohort",
           xlab = "Progress-free Survival Time",
           ylab = "Progress-free survival Probability",
           palette = c("turquoise", "deeppink"),  # Color palette for groups
           legend.title = "Risk Group", 
           legend.labs = c("Low Risk", "High Risk"))

#time roc for pfs risc score #################
time_pfs <- small_geo_df$days_to_tumor_recurrence
event_pfs <- small_geo_df$recurrence_status_f
t_eval <- c(365, 1095, 1825, 2920)  # time points
roc_result3 <- timeROC(
  T = time_pfs,       # Survival time from df
  delta = event_pfs, # Event indicator from df
  marker = small_geo_df[, "risk_scores"], # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result3 #for risk score
#time roc for pfs for separate biomarkers ##################
exprs.df <- small_geo_df[, (colnames(small_geo_df) %in% expression)]
dim(exprs.df)

rez_list2 <- apply(exprs.df, 2, timeROC,
                  T = time_pfs,       # Survival time from df
                  delta = event_pfs, # Event indicator from df
                  #marker  # Predictor already in the df
                  cause = 1,         # Event of interest
                  times = t_eval,    # Time points for ROC
                  iid = TRUE )        # Compute confidence intervals)

auc_table2 <- map_dfr(names(rez_list2), function(gene) {
  roc <- rez_list2[[gene]]
  
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

auc_table2

write.csv(auc_table2, "~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/geo9891_auc_table_pfs.csv")

##pfs plot at year 5 ###############
par(pty = "s")
# Choose target time
target_time <- 1825        
time_index <- which(rez_list2[[1]]$times == target_time)
# Set up base plot with gene 1
plot(
  rez_list2[[1]]$FP[, time_index],
  rez_list2[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity (FPR)",
  ylab = "Sensitivity (TPR)",
  main = paste("Time-dependent PFS ROC Curves at", target_time, "days, in GSE9891"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)
# Add ROC lines for all genes
for (i in 2:length(rez_list2)) {
  lines(
    rez_list2[[i]]$FP[, time_index],
    rez_list2[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}
# Add risk score ROC line in bold
lines(
  roc_result3$FP[, time_index],
  roc_result3$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)
# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")
# legend names
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)
# Add legend 
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list2), "maroon"),
  lwd = c(rep(2, length(rez_list2)), 3),
  cex = 0.6,
  bty = "n"
)

#boxplot borderline vs malignant ###################
colnames(small_geo_df)
#shapriro test?
shapiro_results <- small_geo_df[, c(11:20, 2)] %>%
  pivot_longer(cols = -sample_type , names_to = "gene", values_to = "value") %>%
  group_by(sample_type , gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results #half is not normal:CDCA5 , GRB7 , PKP3, PPT2  , TCEAL4  , VPS33B 
#var test
var_results <- small_geo_df[, c(11:20, 2)] %>%
  pivot_longer(cols = -sample_type, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[sample_type == unique(sample_type)[1]], 
                       value[sample_type == unique(sample_type)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results #EXO1, RAD50, LUC7L2 normal but unequal variances, zfpl1 is normal and equal  

#make long df of lmp vs malignant
#get long df
geo_table_long <- reshape2::melt(small_geo_df[, colnames(small_geo_df) %in%
                                             c("sample_type", expression )],
                                 id.vars="sample_type",
                                 measure.vars= expression)
#get t test:
t.test_geo <- geo_table_long %>%
  filter(variable %in% c("ZFPL1"))%>%
  group_by(variable) %>%
  t_test(value ~ sample_type,
         #p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=F 
  )%>% 
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance()# Format p-values to remove scientific notation
t.test_geo
#get t test, unequal variance
t.test_geo2 <- geo_table_long %>%
  filter(variable %in% c("EXO1", "RAD50", "LUC7L2" ))%>%
  group_by(variable) %>%
  t_test(value ~ sample_type,
         #p.adjust.method = "BH", 
         var.equal = FALSE, # Welch’s t-test 
         paired = FALSE, 
         detailed=F 
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance()# Format p-values to remove scientific notation
t.test_geo2
#get t test, non normal dist. CDCA5 , GRB7 , PKP3, PPT2  , TCEAL4  , VPS33B
wilcox.test_geo <- geo_table_long %>%
  filter(variable %in% c("CDCA5", "PKP3", "GRB7", "PPT2", "VPS33B", "TCEAL4"))%>%
  group_by(variable) %>%
  summarise(
    p_value = wilcox.test(value ~ sample_type)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  add_significance()  # Keep significant results
wilcox.test_geo

t.test_geo_tibble  <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "borderline",   "tumor", 0.00018 , 9, "ZFPL1", #stjudent's
  
  "borderline",   "tumor", 9.45e-23, 9, "EXO1", #welch's
  "borderline",   "tumor", 6.40e-3, 10, "RAD50", #welch's
  "borderline",   "tumor", 2.53e-1, 11, "LUC7L2", #welch's
  
  "borderline",   "tumor", 2.47e-2, 9, "PPT2",  #mann-whitney's
  "borderline",   "tumor", 2.68e-4, 10, "PKP3", #mann-whitneys
  "borderline",   "tumor", 4.76e-10, 10, "CDCA5",#mann-whitney's
  "borderline",   "tumor", 1.31e-1, 10, "VPS33B", #mann-whitney's
  "borderline",   "tumor", 2.43e-2, 12, "GRB7", #mann-whitney's
  "borderline",   "tumor", 3.43e-7, 13, "TCEAL4" #mann-whitney's
)
#leave 3 digits after .:
t.test_geo_tibble$p_custom <- ifelse(t.test_geo_tibble$p.adj < 0.001, 
                            "p < 0.001", 
                            paste0("p = ", sprintf("%.3f",
                            t.test_geo_tibble$p.adj)))
#get colors 
custom_colors <- c("borderline" = "darkblue","tumor" = "darkred") 

#plot
geo_plot <- ggplot(geo_table_long, aes(x=sample_type , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = sample_type )) +
  geom_jitter(aes(color = sample_type ), size=1, alpha=0.5) +
  ylab(label = expression("Normalised expression")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_geo_tibble, label = "p_custom") + #pvalue
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

geo_plot

#boxplot with grade#############################
small_geo_df$grade <- as.factor(small_geo_df$grade)
#shapriro test with grade
shapiro_results_grade <- small_geo_df[, c(11:20, 5)] %>%
  pivot_longer(cols = -grade, names_to = "gene", values_to = "value") %>%
  group_by(grade , gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
print(shapiro_results_grade, n = 40) #not normal: CDCA5, EXO1,PPT2, GRB7, TCEAL4, VPS33B, 

#get long df
geo_table_long_grade <- reshape2::melt(small_geo_df[, colnames(small_geo_df) %in%
                                                   c("grade", expression )],
                                       id.vars="grade",
                                       measure.vars= expression)
#get t test:
t.test_geo <- geo_table_long_grade %>%
  filter(!variable %in% c("CDCA5", "EXO1","PPT2", "GRB7", "TCEAL4", "VPS33B"  ))%>%
  group_by(variable) %>%
  t_test(value ~ grade,
         #p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=F 
  )%>% 
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance() %>% # Format p-values to remove scientific notation
  filter(p_adj < 0.1)
t.test_geo 

#get t test, non normal dist. multiple group
wilcox.test_geo <- geo_table_long_grade %>%
  filter(variable %in% c("CDCA5", "EXO1","PPT2", "GRB7", "TCEAL4", "VPS33B")) %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ grade,
                       p.adjust.method = "BH")  %>% 
  filter(p.adj  < 0.05)
wilcox.test_geo

#make a tibble
t.test_geo_tibble_grade  <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "2",   "3", 0.032, 10, "RAD50", #stjudent's
  "1",   "2", 0.02  , 9, "PKP3", #stjudent's
  "1",   "3", 0.02  , 10, "PKP3", #stjudent's
  
  "1",   "2", 0.0000555 , 9.5, "EXO1", #mann-whitney's
  "1",   "3", 0.00000232, 10, "EXO1", #mann-whitney's
  "2",   "3", 0.01, 10.5, "EXO1", #mann-whitney's
  "1",   "2", 0.000015, 9.5, "CDCA5", #mann-whitney's
  "1",   "3", 0.00000063, 10, "CDCA5", #mann-whitney's
  "2",   "3", 0.016, 10.5, "CDCA5", #mann-whitney's
  
  "2",   "3", 0.031, 9.5, "VPS33B", #mann-whitney's
  "1",   "2", 0.000402, 12, "TCEAL4", #mann-whitney's
  "1",   "3", 0.000614, 12.5, "TCEAL4", #mann-whitney's
  

  
)
#leave 3 digits after .:
t.test_geo_tibble_grade$p_custom <- ifelse(t.test_geo_tibble_grade$p.adj < 0.001, 
                                           "p < 0.001", 
                                           paste0("p = ", sprintf("%.3f",
                                          t.test_geo_tibble_grade$p.adj)))

#get colors 
custom_colors <- c("1" = "darkblue","2" = "darkred","3" = "deeppink", "NA" = "grey") 

#plot
geo_plot_grade <- ggplot(geo_table_long_grade, aes(x=grade , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = grade )) +
  geom_jitter(aes(color = grade ), size=1, alpha=0.5) +
  ylab(label = expression("Normalised expression")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_geo_tibble_grade, label = "p_custom") + #pvalue
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

geo_plot_grade
#boxplot with stage #########################
small_geo_df$tumorstage <- as.factor(small_geo_df$tumorstage)
#shapriro test with stage 
shapiro_results_stage <- small_geo_df[, c(11:20, 4)] %>%
  pivot_longer(cols = -tumorstage, names_to = "gene", values_to = "value") %>%
  group_by(tumorstage , gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
#not normal: EXO1, GRB7   PPT2  LUC7L2  CDCA5 PKP3  VPS33B ZFPL1  
#EXO1 stage III,GRB7 stage III, IV, PPT2  stage III  
results <- lapply(expression, function(var) {
  formula <- as.formula(paste(var, "~ tumorstage"))
  leveneTest(formula, data = small_geo_df)
})
names(results) <- expression
print(results) # RAD50, TCEAL4 will be stujudents test

#get long df
geo_table_long_stage <- reshape2::melt(small_geo_df[, colnames(small_geo_df) %in%
                                                   c("tumorstage", expression )],
                                       id.vars="tumorstage",
                                       measure.vars= expression)
#get t test:
t.test_geo <- geo_table_long_stage %>%
  filter(variable %in% c( "TCEAL4", "RAD50"  ))%>%
  group_by(variable) %>%
  t_test(value ~ tumorstage,
         #p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=F 
  )%>% 
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance() %>% # Format p-values to remove scientific notation
  filter(p_adj < 0.1)
t.test_geo  #TCEAL4  

#get t test, non normal dist. multiple group
wilcox.test_geo <- geo_table_long_stage %>%
  filter(!variable %in% c( "TCEAL4", "RAD50"  ))%>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ tumorstage,
                       p.adjust.method = "BH")  %>% 
  filter(p.adj  < 0.1)
wilcox.test_geo 

#make a tibble
t.test_geo_tibble_stage  <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "1",   "3", 0.053 , 10, "EXO1", 
  "2",   "3", 0.019  , 10, "PKP3",
  "2",   "4", 0.019  , 11, "PKP3", 
  "1",   "3", 0.003  , 10, "CDCA5", 
  
  "1",   "2", 0.058  , 12.5, "TCEAL4", 
  "1",   "3", 0.007   , 12, "TCEAL4",
  "1",   "4", 0.058  , 13, "TCEAL4"
)

#get colors 
custom_colors <- c("1" = "darkblue","2" = "darkred",
                   "3" = "deeppink", "4" = "darkgreen") 

#plot
geo_plot_stage <- ggplot(geo_table_long_stage, aes(x=tumorstage , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumorstage )) +
  geom_jitter(aes(color = tumorstage ), size=1, alpha=0.5) +
  ylab(label = expression("Normalised expression")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_geo_tibble_stage, label = "p.adj") + #pvalue
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

geo_plot_stage

