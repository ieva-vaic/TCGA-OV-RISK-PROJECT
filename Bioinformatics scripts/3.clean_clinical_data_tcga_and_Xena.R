#THIS IS TCGA-OV-RISK-GENES project script No. 3 
#Tidy clinical data

# Load packages ##########################################
library(tidyverse)
setwd("~/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs")
#Load TCGA downloaded data################################
pheno <- readRDS("pheno.RDS")
#look at clinical data 
table(pheno$shortLetterCode, useNA = "a") #7
table(pheno$preservation_method, useNA = "a") #8
#remove uninformative columns
non_informative <- c("shortLetterCode", "sample_type_id", #same as others
                     "tumor_descriptor", "state", "is_ffpe", "tissue_type",
                     "days_to_diagnosis", "last_known_disease_status",
                     "tissue_or_organ_of_origin","project_id",
                     "classification_of_tumor", "tumor_grade", 
                     "progression_or_recurrence", "alcohol_history", "gender",
                     "name", "releasable","released","days_to_sample_procurement", 
                     "composition","preservation_method",
                     #all the same, or nas
                     "ethnicity","prior_malignancy", "oct_embedded", 
                     #only 9 cases
                     "disease_type", "primary_site","treatments",
                     #separate dfs 
                     "pathology_report_uuid","diagnosis_id", "exposure_id",
                     "demographic_id","sample_submitter_id", "patient_names", #extra ids
                     "year_of_diagnosis" #not informative 
                      )
new_pheno <- pheno[ ,!names(pheno) %in% non_informative] #new pheno data with no non-informative columns
#fix treatments################################################################
#create a new df of just tretments column
treatments_dataframe <- bind_rows(pheno$treatment, .id = "id") 
#add and fix patient names
treatments_dataframe$patient_names<- treatments_dataframe$submitter_id  
treatments_dataframe$patient_names <-gsub("_treatment.*", "", treatments_dataframe$patient_names)
#drop uninformative columns
drop_from_treamtment <- c("days_to_treatment_end", "days_to_treatment_start", "treatment_id",
                          "regimen_or_line_of_therapy", "treatment_effect", "therapeutic_agents", "initial_disease_status", 
                          "treatment_intent_type", "treatment_anatomic_site", "treatment_outcome", "state",
                          "created_datetime", "updated_datetime", "id", "submitter_id")
drop_from_treamtment %in% colnames(treatments_dataframe) 
treatments_dataframe <- treatments_dataframe[ ,!names(treatments_dataframe) %in% drop_from_treamtment]
#every person has 2 rows for treatments except for those who had no treatment
#drop the rows where the treatment was not given (indicated by no/na in treatment_or_therapy)
table(treatments_dataframe$treatment_or_therapy, useNA='a') #429 yes
had_treatment <- subset(treatments_dataframe, treatment_or_therapy == 'yes') #new treatments dataframe
table(had_treatment$patient_names, useNA='a') #this showed that some are double
had_treatment$patient <- had_treatment$patient_names #add patient names for easy merging
#paste treatment types to have 1 row per person
had_treatment <- had_treatment %>% group_by(patient_names) %>% mutate(undergone_treatments=paste(sort(treatment_type), collapse="_"))   
had_treatment_collaped <- had_treatment[!duplicated(had_treatment$patient_names), ] # remove duplications, resulting smaller dataframe
#merge treatments and normal pheno data frame (final pheno dataframe)
pheno_final <- merge(x = new_pheno, y = had_treatment_collaped, by = "patient", all = TRUE) 
dim(pheno_final)#429 patients & 35 clinical features
#clean up dataspace
rm(new_pheno)
rm(treatments_dataframe) 
rm(had_treatment)
rm(had_treatment_collaped)
rm(drop_from_treamtment)
rm(non_informative)
#recikact columns are remove extra ids
pheno_final2 <- pheno_final %>%
  relocate(definition, .after = sample_type)  %>%
  relocate(primary_diagnosis, .after = icd_10_code) %>%
  relocate(treatment_type, .before = undergone_treatments) %>%
  select(-(submitter_id))%>%
  relocate(days_to_last_follow_up, .before = days_to_death) %>%
  relocate(age_at_diagnosis, .before = age_at_index) %>%
  relocate(bcr_patient_barcode, .before = sample) %>%
  relocate(sample.aux, .before = sample) 
dim(pheno_final2) #32 columns left
#saveRDS(pheno_final2, "pheno_no_empty_data.RDS")

#Add clinical from XENA database#####################################
#load XENA data:
XENAclin <- read_csv("~/rprojects/TCGA-OV-data/00_ClinTraits.csv")
XENAclin <- XENAclin[, -1]
XENAclin$sample.aux <- XENAclin$sampleID
dim(XENAclin) #584 genes 36 vars
XENAclin <- XENAclin[(XENAclin$sample.aux %in% pheno_final2$sample.aux), ] 
dim(XENAclin)#cases left: 420
joined_clin <- left_join(pheno_final2, XENAclin, by = "sample.aux" ) 
dim(joined_clin)#9 does not have the data!
saveRDS(joined_clin, "joinedTCGA_XENA_clinical2025.RDS")

