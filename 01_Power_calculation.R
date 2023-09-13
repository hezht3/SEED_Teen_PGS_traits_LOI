#####################################
# Power and sample size calculation #
# Author: Johnathan He              #
# Date: September 5, 2023           #
#####################################

setwd("/dcs05/ladd/NDEpi/data/Projects/SEED/LOI/17/SEED_teen_phenotype_LOI")

require(tidyverse)
require(gtsummary)
require(avengeme)


#########
# Input #
#########
# SEED teen sample list
SEED_teen_participant <- read.csv(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/RDA/TEEN/CSV data", "VIEW_PARTICIPANT.csv"))
SEED_teen_hds_1a <- read.csv(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/RDA/TEEN/CSV data", "HDS_HS_1A.csv"))
SEED_teen_hds_2a <- read.csv(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/RDA/TEEN/CSV data", "HDS_HS_2A.csv"))
SEED_teen_hds_3a <- read.csv(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/RDA/TEEN/CSV data", "HDS_HS_3A.csv"))

# SEED cleaned GWAS fam file
SEED1_fam <- read_table(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/SNPArray_DataFiles/CleanMeasuredGenotypes/cleaned2020", "wave1.fam"),
                        col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"))
SEED2_fam <- read_table(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/SNPArray_DataFiles/CleanMeasuredGenotypes/cleaned2020", "wave2.fam"),
                        col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"))
SEED3_fam <- read_table(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/SNPArray_DataFiles/CleanMeasuredGenotypes/cleaned2020", "wave3.fam"),
                        col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"))

# Merge and limit to individuals with genetic data
SEED_teen <- SEED_teen_participant %>% 
    left_join(SEED_teen_hds_1a, by = "FamilyID") %>% 
    left_join(SEED_teen_hds_2a, by = "FamilyID") %>% 
    left_join(SEED_teen_hds_3a, by = "FamilyID") %>% 
    filter(FamilyID %in% c(SEED1_fam$FID, SEED2_fam$FID, SEED3_fam$FID))

# Wave
SEED_teen <- SEED_teen %>% 
    mutate(wave = case_when(FamilyID %in% SEED1_fam$FID ~ 1,
                            FamilyID %in% SEED2_fam$FID ~ 2,
                            FamilyID %in% SEED3_fam$FID ~ 3)) %>% 
    mutate(wave = factor(wave))

# Final classification
SEED_FC <- read.csv("/dcs05/ladd/NDEpi/data/Projects/SEED/LOI/44/output/phenoFiles/moms_kids.allData.all_pheno_PC_data.csv") %>% 
    select(FID, DR_FC)
SEED_teen <- SEED_teen %>% 
    mutate(FamilyID = as.character(FamilyID)) %>% 
    left_join(SEED_FC %>% distinct(FID, .keep_all = TRUE) %>% rename("FamilyID" = "FID"), by = "FamilyID")

# Clean selected outcome variables
### Binary outcomes
SEED_teen <- SEED_teen %>% 
    # ADHD
    mutate(ADHD = ifelse(HS_A9AA_ADHD == "." | HS_A9AA_ADHD == "99", as.numeric(NA), as.numeric(HS_A9AA_ADHD))) %>% 
    # Anxiety
    mutate(Anxiety = ifelse(HS_A9FA_ANX == "." | HS_A9FA_ANX == "99", as.numeric(NA), as.numeric(HS_A9FA_ANX))) %>% 
    # Depression
    mutate(Depression = ifelse(HS_A9RA_DEP == "." | HS_A9RA_DEP == "99", as.numeric(NA), as.numeric(HS_A9RA_DEP))) %>% 
    # Speech or other language disorder
    mutate(Speech_disorder = ifelse(HS_A9AMA_SPE == "." | HS_A9AMA_SPE == "99", as.numeric(NA), as.numeric(HS_A9AMA_SPE))) %>% 
    # Sensory integration disorder
    mutate(Sensory_disorder = ifelse(HS_A9AEA_SENS == "-1" | HS_A9AEA_SENS == "." | HS_A9AEA_SENS == "99", as.numeric(NA), as.numeric(HS_A9AEA_SENS)))

### BMI
SEED_teen <- SEED_teen %>% 
    # Height
    mutate(height = case_when(HS_A1_HT %in% c("-1", ".", "99") ~ as.numeric(NA),
                              HS_A1_HT == 1 ~ as.numeric(HS_A1_HTTAPE),
                              HS_A1_HT == 2 ~ as.numeric(HS_A1_HTRECALL))) %>% 
    # Weight
    mutate(weight = case_when(HS_A2_WT %in% c("-1", ".", "99") ~ as.numeric(NA),
                              HS_A2_WT == 1 ~ as.numeric(HS_A2_WTHOME),
                              HS_A2_WT == 2 ~ as.numeric(HS_A2_WTRECALL))) %>% 
    # BMI: formula is weight (lb) / [height (in)]2 x 703
    mutate(BMI = (weight/(height^2)) * 703) %>%    # Formula: weight (lb) / [height (in)]2 x 703
    mutate(BMI_cat = case_when(BMI < 18.5 ~ 1,                 # Underweight
                               BMI >= 18.5 & BMI < 25.0 ~ 2,   # Normal weight
                               BMI >= 25.0 ~ 3))               # Overweight or obesity

### Verbal communication
SEED_teen <- SEED_teen %>%  # Verbal communications
    mutate(Verbal_com = ifelse(HS_A7_VERBAL == "-1"|HS_A7_VERBAL == "."|is.na(HS_A7_VERBAL), as.numeric(NA), as.numeric(HS_A7_VERBAL))) %>% 
    mutate(Verbal_com = case_when(is.na(Verbal_com) ~ as.numeric(NA),
                                  Verbal_com == 1 ~ 0,
                                  Verbal_com %in% 2:5 ~ 1))

# Limit to participants with at least one outcome measurements
SEED_teen <- SEED_teen %>% 
    filter(if_any(c(ADHD, Anxiety, Depression, Speech_disorder, Sensory_disorder, BMI_cat, Verbal_com), ~ !is.na(.x)))   # 392 participants
#save(SEED_teen, file = "./Intermediate/SEED_teen.rda")

# Table for sample size
SEED_teen %>% 
    mutate(DR_FC = factor(DR_FC, levels = c("CASE", "POP", "DD"))) %>% 
    select(wave, DR_FC) %>% 
    tbl_summary(by = "DR_FC") %>% 
    add_overall()


#####################
# Power calculation #
#####################
# ASD & ADHD
polygenescore(1000000, 
              n = c(18381+27969, 392), 
              vg1 = 0.025, 
              cov12 = 0.38, 
              binary = TRUE, 
              sampling = c(18381/(8381+27969), 125/392), 
              pupper = c(0, 1), 
              prevalence = 1/44)

# ASD & Anxiety
polygenescore(1000000, 
              n = c(18381+27969, 392), 
              vg1 = 0.025, 
              cov12 = 0.22, 
              binary = TRUE, 
              sampling = c(18381/(8381+27969), 125/392), 
              pupper = c(0, 1), 
              prevalence = 1/44)

# ASD & Depression
polygenescore(1000000, 
              n = c(18381+27969, 392), 
              vg1 = 0.025, 
              cov12 = 0.54, 
              binary = TRUE, 
              sampling = c(18381/(8381+27969), 125/392), 
              pupper = c(0, 1), 
              prevalence = 1/44)
