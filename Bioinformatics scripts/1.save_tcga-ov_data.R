#THIS IS OV-RISK-GENES project script No. 1 
#Downloading the TCGA data with TCGAbiolinks

# Load packages ##########################################
library("TCGAbiolinks")
library("SummarizedExperiment")
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/") #use data directory

## Build your query ##########################################
query_TCGA = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", 
  access = "open" )

## Download your data ##########################################
# (do not forget to setwd to data folder, creates a big folder)
GDCdownload(query = query_TCGA, files.per.chunk = 200) 

## Prepare your data ##########################################
tcga_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)

#get clinical data file
pheno <- as.data.frame(colData(tcga_data)) 
# *there will be additional clinical data from XENA database

## Save gene expression and clinical data ##########################################
#final saved TCGA-OC data: full file RangedSummarizedExperiment tcga_data and clinical data (colData)
#you may need to remove all of the downloaded files, because it takes a lot of computer space
#gene expression data:
#saveRDS(object = tcga_data,
#        file = "tcga_data.RDS",
#        compress = FALSE)
#clincial data:
#saveRDS(object = pheno,
#        file = "pheno.RDS",
#        compress = FALSE)

