#THIS IS TCGA-OV-RISK-GENES project script No. 6
#connect all tcga and gtex mrna data together, filter for protein coding genes

# Load packages ##########################################
library(tidyverse)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load data ###################################
gtex_counts <- readRDS("gtex_counts.RDS") #GTEx ovarian cancer v8 mRNA counts
tcga_counts <- readRDS("tcga_with_names_all.RDS") #TCGA-OV project mRNA counts

#the gtex is now a matrix that looks like this:
#rows: ENSEMBLIDS (transcipt level), all 56200 of them 
head(colnames(gtex_counts))
dim(gtex_counts) #56200 counts and 182 cases
#make a numeric dataframe from gtex
gtex_df <- as.data.frame(gtex_counts)
#add row names (ENSG ids) to gtex_df
rownames(gtex_df) <- gtex_df$Name
#make numeric
gtex_df[,3:182] <- lapply(gtex_df[,3:182], as.numeric)
dim(gtex_df) #56200 genes, 180 cases (1 and 2 cols are gene names)
#equalize colnames with tcga_df
gtex_df$external_gene_name <- gtex_df$Description
gtex_df$ensembl_gene_id <- gsub("\\..*", "",gtex_df$Name)
sum(gtex_df$ensembl_gene_id %in% tcga_counts$ensembl_gene_id) #55661 the same
#find duplicate genes
dup_tcga <- duplicated(tcga_counts$ensembl_gene_id)
dup_tcga <- tcga_counts[dup_tcga, ] 
dup_tcga$ensembl #44 parY
dup_gtex <- duplicated(gtex_df$ensembl_gene_id) 
dup_gtex <- gtex_df[dup_gtex, ]
rownames(dup_gtex) #44 parY
#remove par Y from tcga
ENTG_Y <- tcga_counts[grepl('_PAR_Y', tcga_counts$ensembl), ]
ENTG_Y_names <- ENTG_Y$ensembl
tcga_counts_filt <- tcga_counts[!(tcga_counts$ensembl %in% ENTG_Y_names), ]
rownames(tcga_counts_filt) <- tcga_counts_filt$ensembl_gene_id #
dim(tcga_counts_filt) #lieka  60616   420  # new
#remove par Y from gtex
ENTG_YG <- gtex_df[grepl('_PAR_Y', gtex_df$Name), ]
ENTG_YG_names <- ENTG_YG$Name #turbut tie patys bet anyways
gtex_counts_filt <- gtex_df[!(gtex_df$Name %in% ENTG_YG_names), ]
rownames(gtex_counts_filt) <- gtex_counts_filt$ensembl_gene_id #veik rownames 
dim(gtex_counts_filt) #lieka  56156   184
#join TCGA and GTEx
gtcga <- left_join(gtex_counts_filt, tcga_counts_filt, by = "ensembl_gene_id")
dim(gtcga) #56156   603

#some genes are only in GTEX:
sum(is.na(gtcga$external_gene_name.y)) #this amount of genes are missing from tcga: #687
gtcga_na <- gtcga[is.na(gtcga$external_gene_name.y), ] #df of genes without external_gene_name.y
#remove them:
gtcga_final <- gtcga[!(is.na(gtcga$external_gene_name.y)), ] 
dim(gtcga_final)#final gene count #55469    

#remove non-protein-coding ###########################
#get gene functions:
biotypes <- table(gtcga_final$gene_biotype, useNA = "a")
biotypes <- as.data.frame(biotypes)
biotypes
#filter for protein coding:
gtgca_protein <- gtcga_final[(gtcga_final$gene_biotype == "protein_coding" ), ]
dim(gtgca_protein) #left 19197 genes 
table(gtgca_protein$gene_biotype, useNA = "a") 

#find repeating names
gtgca_protein["external_gene_name.y"][gtgca_protein["external_gene_name.y"] == ''] <- NA
sum(is.na(gtgca_protein$external_gene_name.y)) #29 without names
#remove empty names
gtgca_protein <- gtgca_protein %>% 
  mutate(external_gene_name.y = coalesce(external_gene_name.y, Description))  
any(is.na(gtgca_protein$external_gene_name.y)) #no empty names left
#can't make the gene names rownames due to "NPIPA9", remove
which(gtgca_protein$external_gene_name.y == 'NPIPA9')
gtgca_protein[13803, 602] <- NA
gtgca_protein[13807, 602] <- NA
gtgca_protein <- gtgca_protein %>% 
  mutate(external_gene_name.y = coalesce(external_gene_name.y, Description))  
any(is.na(gtgca_protein$external_gene_name.y))
#make gene names rownames
rownames(gtgca_protein) <- gtgca_protein$external_gene_name.y

#loose the gene descriptions now
colnames(gtgca_protein)
gtgca_final_no_names <- gtgca_protein[, -c(1:2, 183:184, 601:603)] 
dim(gtgca_final_no_names) #19197,   596
#save
#saveRDS(gtgca_final_no_names, "gtcga_final_counts.RDS")
