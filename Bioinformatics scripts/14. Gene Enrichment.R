#THIS IS TCGA-OV-RISK-GENES project script No. 13
#Describe main genes

# Load packages ##########################################
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 

# Get annotations for 214 genes ##########################
gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = gtex_genes,
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL"))

head(gene_info)
# Get annotations for main 10 genes ##########################
#get genes of interest
expression <- c( "EXO1",   "RAD50",  "PPT2",   "LUC7L2", "PKP3",
                 "CDCA5",  "ZFPL1" , "VPS33B", "GRB7",   "TCEAL4")
gene_info10 <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = expression,
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL", "GENENAME", "ENTREZID", "ENSEMBL"))

head(gene_info10)

#functional enrichment 214 genes#################################
# Get Entrez IDs
entrez_ids <- gene_info$ENTREZID

# GO enrichment (Biological Process)
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment for LASSO selected genes") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#functional enrichment top 10 genes#################################
# Get Entrez IDs
entrez_ids2 <- gene_info10$ENTREZID

# GO enrichment (Biological Process)
ego2 <- enrichGO(gene = entrez_ids2,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)
#0 enriched terms

#KEGG for 214 LASSO genes ########################
ekegg <- enrichKEGG(gene = entrez_ids, organism = "hsa")
dotplot(ekegg, showCategory = 10) +
  ggtitle("KEGG Enrichment for LASSO selected genes") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Get MSigDB gene sets ##############
# Use Hallmark (category "H"), or you can change this to C2, C5, etc.
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# Prepare TERM2GENE format for enricher()
msig_term2gene <- msig_hallmark[, c("gs_name", "entrez_gene")]

#enrichment
enrich_res <- enricher(gene = entrez_ids,
                       TERM2GENE = msig_term2gene,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
dotplot(enrich_res, showCategory = 10) +
  ggtitle("MSigDB Hallmark Enrichment for LASSO selected genes") +
  theme_minimal()

#Get MSigDB gene sets from other collections ##############
View(msigdbr_collections())
# C4 - 	Oncogenic modules
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "C4")

# Prepare TERM2GENE format for enricher()
msig_term2gene <- msig_hallmark[, c("gs_name", "entrez_gene")]

#enrichment
enrich_res <- enricher(gene = entrez_ids,
                       TERM2GENE = msig_term2gene,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
dotplot(enrich_res, showCategory = 10) +
  ggtitle("MSigDB Cancer Modules Enrichment for LASSO selected genes") +
  theme_minimal()
head(enrich_res@result)

#Get MSigDB gene sets from other collections ##############
# C2 - 	Kegg
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "C2")

# Prepare TERM2GENE format for enricher()
msig_term2gene <- msig_hallmark[, c("gs_name", "entrez_gene")]

#enrichment
enrich_res <- enricher(gene = entrez_ids,
                       TERM2GENE = msig_term2gene,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
dotplot(enrich_res, showCategory = 10) +
  ggtitle("MSigDB Cancer Modules Enrichment for LASSO selected genes") +
  theme_minimal()
head(enrich_res@result)

#Get MSigDB gene sets from other collections ##############
# C2 - 	Kegg
# Get annotations for main 10 genes ##########################
#get genes of interest
expression <- c( "EXO1",   "RAD50",  "PPT2",   "LUC7L2", "PKP3",
                 "CDCA5",  "ZFPL1" , "VPS33B", "GRB7",   "TCEAL4")
entrez_ids10 <- gene_info10$ENTREZID
#enrichment
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "C2")
enrich_res <- enricher(gene = entrez_ids10,
                       TERM2GENE = msig_term2gene,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
dotplot(enrich_res, showCategory = 10) +
  ggtitle("MSigDB KEGG Enrichment for 10 selected genes") +
  theme_minimal()
head(enrich_res@result)