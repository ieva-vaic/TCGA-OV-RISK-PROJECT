#THIS IS OV-RISK-GENES project script No. 2
#Saving the GTEx data and additional XENA data
#this script will read the gtex to a manageable data format
# Load packages ##########################################
library(phantasus)
library(SummarizedExperiment)
library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/") #use data directory

# READ GTEx FILE ###############################################
#data downloaded form https://www.gtexportal.org/home/datasets; 
#GTEx Analysis V8 release; counts by tissue: file name: gene_reads_2017-06-05_v8_ovary.gct.gz
#first, gz the file on your system
#then, open the file here: 
gtex <- read.gct("gene_reads_2017-06-05_v8_ovary.gct")
#get expression:
gtex_counts <- exprs(gtex)
# save for later:
#saveRDS(gtex_counts, "gtex_counts.RDS")

#READ XENA FILE #############################################
#XENA database has the same samples but more clinical data
#protocol on downloading the data can be found here: https://doi.org/10.1016/j.xpro.2022.101168 
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/XENA/")
data(XenaData)
write.csv(XenaData, "00_tblXenaHubInfo.csv")

GeneExpectedCnt_toil = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")

XenaQuery(GeneExpectedCnt_toil) %>%
  XenaDownload(destdir = "./") #download once! slow!

paraCohort = "TCGA Ovarian Cancer"; #Selecting the OVARIAN Cancer cohort.
paraDatasets = "TCGA.OV.sampleMap/OV_clinicalMatrix"; #Selecting the OVARIAN Cancer clinical matrix.
Clin_TCGA = XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterCohorts = paraCohort) %>%
  XenaFilter(filterDatasets = paraDatasets);
XenaQuery(Clin_TCGA) %>%
  XenaDownload(destdir = "./") #download data

#subset to get only OVARIAN data
filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz");
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01));
paraStudy = "GTEX"; #Setting "GTEx" as the study of interest.
paraPrimarySiteGTEx = "Ovary"; #Setting "Ovary" as the primary site of interest.
paraPrimaryTissueGTEx = "^Ovary"; #Setting "Ovary" as the primary tissue of interest.
filterGTEx02 = subset(filterGTEx01,
                      study == paraStudy &
                        primarysite == paraPrimarySiteGTEx &
                        grepl(paraPrimaryTissueGTEx, filterGTEx01$`primary disease or tissue`))
filterTCGA01 = fread(paraDatasets);
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01));
paraSampleType = "Primary Tumor"; #Setting "Primary Tumor" as the sample type of interest.
paraPrimarySiteTCGA = "Ovary"; #Setting "Ovary" as the primary site of interest.
paraHistologicalType = "Serous Cystadenocarcinoma"; #Setting "Serous Cystadenocarcinoma" as the histological type of interest.
filterTCGA02 = subset(filterTCGA01,
                      sampletype == paraSampleType &
                        primarysite == paraPrimarySiteTCGA &
                        grepl(paraHistologicalType, filterTCGA01$histologicaltype))


#CLEAN XENA CLINICAL DATA ###########################################
#Keep variable "Lymphatic Invasion".
varClinKeep = c("sampleID", "lymphaticinvasion", "clinicalstage", "vitalstatus", "ageatinitialpathologicdiagnosis", "neoplasmhistologicgrade", 
                "additionalpharmaceuticaltherapy", "additionalradiationtherapy",
                "anatomicneoplasmsubdivision", "daystobirth",
                "daystodeath", "daystoinitialpathologicdiagnosis", "daystolastfollowup", 
                "daystonewtumoreventadditionalsurgeryprocedure",
                "daystonewtumoreventafterinitialtreatment",
                "newneoplasmeventtype",
                "easterncanceroncologygroup", 
                "followupcasereportformsubmissionreason", "followuptreatmentsuccess", 
                "historyofneoadjuvanttreatment", #only one yes
                "intermediatedimension", "longestdimension", "shortestdimension", 
                "newtumoreventadditionalsurgeryprocedure", #tik 7 yes
                "newtumoreventafterinitialtreatment" , #"22no, 67 yes"
                "personneoplasmcancerstatus", 
                "postoperativerxtx",
                "primarytherapyoutcomesuccess",
                "radiationtherapy" ,
                "residualtumor", #r0:15  r1:33   r2:5   rx:3
                "sampletype", "tumorresidualdisease", "venousinvasion")
clinDF01 = as.data.frame(do.call(cbind, filterTCGA02))
clinFinal = clinDF01[varClinKeep]
#Identify observations/samples with no values assigned to the kept variables.
colSums(clinFinal == "")
colSums(is.na(clinFinal)) #rcount NAs
#Replace "no values" with "NA"
NA -> clinFinal[clinFinal == ""]
colSums(is.na(clinFinal))
#For the variable "Lymphatic Invasion", check count of YES/NO.
table(clinFinal$lymphaticinvasion)
#Re-code the variable "Lymphatic Invasion" from YES/NO to 1/0.
clinFinal$lymphaticinvasion = if_else(clinFinal$lymphaticinvasion == "YES",
                                      1, 0, missing = NULL)
#Verify that the count of 1/0 mirrors previous count of YES/NO.
table(clinFinal$lymphaticinvasion)
#For the variable "stage", check counts
table(clinFinal$clinicalstage, useNA = "a")
#Re-code the variable "Stage" from Stage IA.. to.. -2 -1 0 +1.
clinFinal$clinicalstage2= gsub("[ABC]$", "", clinFinal$clinicalstage)
table(clinFinal$clinicalstage2, useNA = "a") 
clinFinal$clinicalstage_num <- clinFinal$clinicalstage2
clinFinal$clinicalstage_num <- recode(clinFinal$clinicalstage_num, "Stage I" = -2, "Stage II" = -1, "Stage III" = 0, "Stage IV" = 1)
table(clinFinal$clinicalstage_num, useNA = "a") 

#Save as csv with XENA clinical data data ####################################
#write.csv(clinFinal, "~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/00_ClinTraits.csv")
