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
# load series and platform data from GEO

gset <- getGEO("GSE9891", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#saveRDS(gset, "GSE9891.RDS")
ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

#make df
ex_df <- as.data.frame(ex)

#get names
probe_ids <- rownames(ex_df)  

# Get gene symbols
gene_symbols <- getSYMBOL(probe_ids, "hgu133plus2.db")

# Add gene symbols to your data
ex_df$GeneSymbol <- gene_symbols

#chek if useful
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5",
                "ZFPL1","VPS33B", "GRB7","TCEAL4")
expression %in% ex_df$GeneSymbol #all in there
gtex_genes <- readRDS("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/gtcga_elastic_2025.RDS")
gtex_genes %in% ex_df$GeneSymbol #NOT ALL in there
# Remove NAs
ex_df <- ex_df[!is.na(ex_df$GeneSymbol), ]

# Optional: Aggregate to one row per gene symbol (e.g., by mean)
expression_by_gene <- ex_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean))
expression_by_gene <- as.data.frame(expression_by_gene)
rownames(expression_by_gene) <- expression_by_gene$GeneSymbol
expression_by_gene <- expression_by_gene[, -1]

#saveRDS(expression_by_gene, "GSE9891expression_by_gene_.RDS")
#make smaler df of my genes only
SMALL_GEO_MX <- expression_by_gene[rownames(expression_by_gene) %in% expression, ] 
#saveRDS(SMALL_GEO_MX, "SMALL_GEO_MX_GSE9891.RDS")
#clinical df##############################################
clin_df <- pData(gset)

#fix clinical
colnames(clin_df)
#choose what to keep
table(clin_df$characteristics_ch1.1, useNA = "a") #lmp - borderline
table(clin_df$characteristics_ch1.2, useNA = "a") #1 adeno
table(clin_df$characteristics_ch1.3, useNA = "a") #stage full
table(clin_df$characteristics_ch1.4, useNA = "a") #grade full
table(clin_df$`Consolidated.Grade :ch1`, useNA = "a") #same as grade 
table(clin_df$`Primary.Site :ch1`, useNA = "a") #same as ch1.2 
table(clin_df$`StageCode :ch1`, useNA = "a") #same as ch1.3 
table(clin_df$`Type :ch1`, useNA = "a") #same as ch1.1

leave_clincial <- c("Consolidated.Grade :ch1", "Primary.Site :ch1", 
  "StageCode :ch1", "Type :ch1")
clin_df <- clin_df[, colnames(clin_df) %in% leave_clincial]
#make heatmap based on type ########################
#create clinical annotation
row_ha2 = rowAnnotation(Type = clin_df$`Type :ch1`,
                        Grade = clin_df$`Consolidated.Grade :ch1`,
                        Stage =clin_df$`StageCode :ch1`,
                        #choose colors
                        col = list(
                                   `Type` = c("LMP" = "#9cd4c4",  
                                              "Malignant" = "#a89cd4"),
                                   Grade = c("1" = "#9cd4c4",  
                                             "2" = "turquoise", 
                                             "3" = "#a89cd4",
                                             "NA" = "grey"),
                                   Stage = c("IA" = "#FFB6C1",  
                                             "IB" = "#FF69B4",
                                             "IC" = "#FF1493",  
                                             
                                             "IIA" = "#E6E6FA", 
                                             "IIB" = "#9370DB",  
                                             "IIC" = "#9932CC",
                                             
                                             "III" = "#F08080",  
                                             "IIIA" = "#DC143C",
                                             "IIIB" = "#B22222",  
                                             "IIIC" = "#8B0000",
                                             
                                             "IV" = "#9cd4c4",  
                                             "NK" = "turquoise",
                                             "NA" = "grey")
                        ))
#MAKE SMALL HEATMAP WITH GETX AGE################
heatmap_geo <- Heatmap(as.matrix(t(SMALL_GEO_MX)) ,  
                         show_row_names = F,
                         row_split = clin_df$`Type :ch1`, 
                         column_names_gp = gpar(fontsize = 6, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Gene Expression"),
                         right_annotation = row_ha2,
                         cluster_rows = F)
heatmap_geo

#make boxplots##################################
#first, join clincal with my data
SMALL_GEO_MX_t <- t(SMALL_GEO_MX) #transorm 
df_joined <- cbind(SMALL_GEO_MX_t, clin_df)
colnames(df_joined) <- c(expression, "Grade", "Primary_site", "Stage", "Type")
saveRDS(df_joined, "df_joined_GSE9891.RDS")
#shapriro test?
shapiro_results <- df_joined[, c(1:10, 14)] %>%
  pivot_longer(cols = -Type, names_to = "gene", values_to = "value") %>%
  group_by(Type, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results #half is not normal: CDCA5, EXO1, GRB7, PPT2, VPS33B 
#var test
var_results <- df_joined[, c(1:10, 14)] %>%
  pivot_longer(cols = -Type, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Type == unique(Type)[1]], 
                       value[Type == unique(Type)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results #LUC7L2, PKP3, RAD50, ZFPL1 normal but unequal variances   
#make long df of lmp vs malignant
#get long df
geo_table_long <- reshape2::melt(df_joined[, colnames(df_joined) %in%
                                  c("Type", expression )],
                                   id.vars="Type",
                                   measure.vars= expression)
#get t test:
t.test_geo <- geo_table_long %>%
  filter(!variable %in% c("LUC7L2", "PKP3", "RAD50", "ZFPL1",
                          "CDCA5", "EXO1", "GRB7", "PPT2", "VPS33B"))%>%
  group_by(variable) %>%
  t_test(value ~ Type,
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
  filter(variable %in% c("LUC7L2", "PKP3", "RAD50", "ZFPL1"))%>%
  group_by(variable) %>%
  t_test(value ~ Type,
         #p.adjust.method = "BH", 
         var.equal = FALSE, # Welch’s t-test 
         paired = FALSE, 
         detailed=F 
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance()# Format p-values to remove scientific notation
t.test_geo2
#get t test, non normal dist. CDCA5, EXO1, GRB7, PPT2, VPS33B
wilcox.test_geo <- geo_table_long %>%
  filter(variable %in% c("CDCA5", "EXO1", "GRB7", "PPT2", "VPS33B"))%>%
  group_by(variable) %>%
  summarise(
    p_value = wilcox.test(value ~ Type)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  add_significance()  # Keep significant results
wilcox.test_geo
#make tibble for p values
t.test_geo_tibble  <- tibble::tribble(
    ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
    "LMP",   "Malignant", 0.0000116, 8, "TCEAL4", #stjudent's
    
    "LMP",   "Malignant", 4.52e-12, 10, "RAD50", #welch's
    "LMP",   "Malignant", 7.24e-4, 10, "LUC7L2", #welch's
    "LMP",   "Malignant", 4.21e-3, 9, "PKP3", #welch's
    "LMP",   "Malignant", 1.89e-5 , 9, "ZFPL1",  #welch's
    
    "LMP",   "Malignant", 4.64e-10, 10, "EXO1",  #mann-whitney's
    "LMP",   "Malignant", 4.33e-2, 12, "PPT2", #mann-whitneys
    "LMP",   "Malignant", 4.61e-2, 9, "CDCA5",#mann-whitney's
    "LMP",   "Malignant", 5.60e-7, 12, "VPS33B", #mann-whitney's
    "LMP",   "Malignant", 6.08e-2, 9, "GRB7" #mann-whitney's

  )
#leave 3 digits after .:
t.test_geo_tibble$p_custom <- ifelse(t.test_geo_tibble$p.adj < 0.001, 
                                      "p < 0.001", 
                                      paste0("p = ", sprintf("%.3f",
                                     t.test_geo_tibble$p.adj)))
#get colors 
custom_colors <- c("LMP" = "darkblue","Malignant" = "darkred") 

#plot
geo_plot <- ggplot(geo_table_long, aes(x=Type , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Type )) +
  geom_jitter(aes(color = Type ), size=1, alpha=0.5) +
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

#with grade#############################
#shapriro test with grade
shapiro_results_grade <- df_joined[, c(1:11)] %>%
  pivot_longer(cols = -Grade, names_to = "gene", values_to = "value") %>%
  group_by(Grade , gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
print(shapiro_results_grade, n = 40) #not normal: CDCA5, EXO1, RAD50, GRB7, PPT2, VPS33B   
#var test 
var_results_grade1_2 <- df_joined[, c(1:11)] %>%
  filter(Grade %in% c(1,2)) %>%
  pivot_longer(cols = -Grade, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grade == unique(Grade)[1]], 
                       value[Grade == unique(Grade)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
print(var_results_grade1_2) #same as not normal

var_results_grade1_3 <- df_joined[, c(1:11)] %>%
  filter(Grade %in% c(1,3)) %>%
  pivot_longer(cols = -Grade, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grade == unique(Grade)[1]], 
                       value[Grade == unique(Grade)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
print(var_results_grade1_3) #PKP3     

var_results_grade2_3 <- df_joined[, c(1:11)] %>%
  filter(Grade %in% c(2,3)) %>%
  pivot_longer(cols = -Grade, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grade == unique(Grade)[1]], 
                       value[Grade == unique(Grade)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
print(var_results_grade2_3) #PKP3    

#get long df
geo_table_long_grade <- reshape2::melt(df_joined[, colnames(df_joined) %in%
                                             c("Grade", expression )],
                                 id.vars="Grade",
                                 measure.vars= expression)
#get t test:
t.test_geo <- geo_table_long_grade %>%
  filter(!variable %in% c("CDCA5", "EXO1", "RAD50", "GRB7", "PPT2", "VPS33B", "PKP3" ))%>%
  group_by(variable) %>%
  t_test(value ~ Grade,
         #p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=F 
  )%>% 
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance() %>% # Format p-values to remove scientific notation
  filter(p_adj < 0.05)
t.test_geo 
#get t test, unequal variance
t.test_geo2 <- geo_table_long_grade %>%
  filter(variable %in% c("PKP3"))%>%
  group_by(variable) %>%
  t_test(value ~ Grade,
         #p.adjust.method = "BH", 
         var.equal = FALSE, # Welch’s t-test 
         paired = FALSE, 
         detailed=F 
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance() %>% # Format p-values to remove scientific notation
  filter(p_adj < 0.05)
t.test_geo2
#get t test, non normal dist. multiple group
wilcox.test_geo <- geo_table_long_grade %>%
  filter(variable %in% c("CDCA5", "EXO1", "RAD50", "GRB7", "PPT2", "VPS33B")) %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ Grade,
                       p.adjust.method = "BH")  %>% 
  filter(p.adj  < 0.05)
wilcox.test_geo
#make a tibble
t.test_geo_tibble_grade  <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "1",   "3", 0.036, 10, "LUC7L2", #stjudent's
  "1",   "2", 0.019, 7.5, "TCEAL4", #stjudent's
  "1",   "3", 0.009, 8, "TCEAL4", #stjudent's
  
  "1",   "2", 0.079, 9, "PKP3", #welch's
  "1",   "3", 0.05, 9.5, "PKP3", #welch's
  
  "1",   "2", 0.0000333, 9.5, "EXO1", #mann-whitney's
  "1",   "3", 0.000000205, 10, "EXO1", #mann-whitney's
  "2",   "3", 0.008, 10.5, "EXO1", #mann-whitney's
  "1",   "2", 0.000274, 9, "RAD50", #mann-whitney's
  "1",   "3", 0.000000906, 9.5, "RAD50", #mann-whitney's
  "2",   "3", 0.005, 10, "RAD50", #mann-whitney's
  "1",   "2", 0.000723, 12, "VPS33B", #mann-whitney's
  "1",   "3", 0.000723, 12.5, "VPS33B" #mann-whitney's
  
)
#leave 3 digits after .:
t.test_geo_tibble_grade$p_custom <- ifelse(t.test_geo_tibble_grade$p.adj < 0.001, 
                                     "p < 0.001", 
                                     paste0("p = ", sprintf("%.3f",
                                                            t.test_geo_tibble_grade$p.adj)))

#get colors 
custom_colors <- c("1" = "darkblue","2" = "darkred","3" = "deeppink", "NA" = "grey") 

#plot
geo_plot_grade <- ggplot(geo_table_long_grade, aes(x=Grade , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grade )) +
  geom_jitter(aes(color = Grade ), size=1, alpha=0.5) +
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

#with stage#############################
table(df_joined$Stage, useNA = "a")
#fix up stage
df_joined$stage <- gsub("^(I[A-C])$", "stage I", df_joined$Stage)
df_joined$stage <- gsub("^(II[A-C])$", "stage II", df_joined$stage)
df_joined$stage <- gsub("^(III[A-C]?)$", "stage III", df_joined$stage)
df_joined$stage <- gsub("^IV$", "stage IV", df_joined$stage)
df_joined$stage[df_joined$stage %in% c("NK")] <- NA
df_joined$stage[df_joined$stage %in% c("NA")] <- NA
table(df_joined$stage, useNA = "a")
#shapriro test with stage
shapiro_results_stage <- df_joined[, c(1:10, 15)] %>%
  pivot_longer(cols = -stage, names_to = "gene", values_to = "value") %>%
  group_by(stage , gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
print(shapiro_results_stage, n = 40) #not normal: LUC7L2 stage1, CDCA5 stage III, IV,
#EXO1 stage III,GRB7 stage III, IV, PPT2  stage III  
results <- lapply(expression, function(var) {
  formula <- as.formula(paste(var, "~ stage"))
  leveneTest(formula, data = df_joined)
})
names(results) <- expression
print(results) #EXO1, RAD50, CDCA5


#get long df
geo_table_long_stage <- reshape2::melt(df_joined[, colnames(df_joined) %in%
                                             c("stage", expression )],
                                 id.vars="stage",
                                 measure.vars= expression)
#get t test:
t.test_geo <- geo_table_long_stage %>%
  group_by(variable) %>%
  t_test(value ~ stage,
         #p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=F 
  )%>% 
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance() %>% # Format p-values to remove scientific notation
  filter(p_adj < 0.1)
t.test_geo  #EXO1 stage III not notmal thus, this is not valid
#get t test, unequal variance
t.test_geo2 <- geo_table_long_stage %>%
  group_by(variable) %>%
  t_test(value ~ stage,
         #p.adjust.method = "BH", 
         var.equal = FALSE, # Welch’s t-test 
         paired = FALSE, 
         detailed=F 
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  add_significance() %>% # Format p-values to remove scientific notation
  filter(p_adj < 0.1)
t.test_geo2 #no significant
#get t test, non normal dist. multiple group
wilcox.test_geo <- geo_table_long_stage %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ stage,
                       p.adjust.method = "BH")  %>% 
  filter(p.adj  < 0.1)
wilcox.test_geo ##not normal:  #EXO1 stage III

#make a tibble
t.test_geo_tibble_stage  <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "stage I",   "stage III", 0.011, 12, "VPS33B", #stjudent's
  "stage I",   "stage III", 0.003, 9.5, "EXO1" #mann-whitney's
)

#get colors 
custom_colors <- c("stage I" = "darkblue","stage II" = "darkred","stage III" = "deeppink", "stage IV" = "darkgreen") 

#plot
geo_plot_stage <- ggplot(geo_table_long_stage, aes(x=stage , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = stage )) +
  geom_jitter(aes(color = stage ), size=1, alpha=0.5) +
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
expression %in% rownames(exprs_matrix) #all in there
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
                    "recurrence_status","days_to_death","vital_status") #choose what to leave
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
#time roc for risk score########################
t_eval <- c(365, 1095, 1825)  # time points
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
write.csv(auc_table, "~/rprojects/TCGA-OV-RISK-PROJECT/Outputs bioinfomatics/geo9891_auc_table_year5.csv")

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



