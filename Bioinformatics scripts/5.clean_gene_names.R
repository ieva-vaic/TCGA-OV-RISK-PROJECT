#THIS IS TCGA-OV-RISK-GENES project script No. 5
#tidy gene names

# Load packages ##########################################
library(tidyverse)
library(biomaRt)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")
# Load TCGA mRNA data ###################################
tcga_data <- readRDS("tcga_selected_counts.RDS") #
#transform
mRNA_counts <- t(tcga_data)
#transform to df
mRNA_counts <- as.data.frame(mRNA_counts) 
dim(mRNA_counts)#60660 transcripts (rows), 416 cases (cols) 
#turn form integer to numeric
mRNA_counts[,1:416] <- lapply(mRNA_counts[,1:416], as.numeric)
str(mRNA_counts) #60660 genes
#gather gene names in a column
mRNA_counts$ensembl <- rownames(mRNA_counts)
#remove ensg id endings to leave just gene names
mRNA_counts$ensembl_gene_id <- gsub("\\..*", "",mRNA_counts$ensembl)
#biomart
listEnsembl() #shows available databases
ensembl <- useEnsembl(biomart = "genes") #choose genes
datasets <- listDatasets(ensembl) #choose human genes
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl') #name of data base ir data set
attr <- listAttributes(ensembl.con) #list atributes
filters <- listFilters(ensembl.con) #list filters

tcga_genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', "gene_biotype"), #values to retreve
                    filters = "ensembl_gene_id", #input on quary
                    values = mRNA_counts$ensembl_gene_id, #my ids
                    mart = ensembl.con) #conection object
dim(tcga_genes) #60419 = 241 lost    

#join dfs
counts_tcga_with_gene_names <- left_join(mRNA_counts, tcga_genes, by= "ensembl_gene_id")
dim(counts_tcga_with_gene_names) #60660   420
#saveRDS(counts_tcga_with_gene_names, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/tcga_with_names_all.RDS")

