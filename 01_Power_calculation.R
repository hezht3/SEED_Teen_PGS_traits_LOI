#####################################
# Power and sample size calculation #
# Author: Johnathan He              #
# Date: September 5, 2023           #
#####################################

setwd("/dcs05/ladd/NDEpi/data/Projects/SEED/LOI/17/SEED_teen_phenotype_LOI")

require(tidyverse)
require(gtsummary)
require(pwrss)


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
# Individual associations
power.rslt <- NULL
for(prob in c(0.08, 0.10, 0.12)) {
    for(or in c(1.5, 1.6, 1.7)) {
        power.rslt.current <- tibble(p0 = prob, odds.ratio = or,
                                     power = pwrss.z.logistic(n = 392, p0 = prob, odds.ratio = or, r2.other.x = 0.10, 
                                                              alpha = 0.05, dist = "normal")$power)
        power.rslt <- bind_rows(power.rslt, power.rslt.current)
    }
}

tiff("./Output/Figure 1. Power calculation.tiff",
     width = 2100, height = 1350, pointsize = 5, res = 300)
power.rslt %>% 
    mutate(p0 = factor(p0, levels = c(0.08, 0.10, 0.12), labels = c("8%", "10%", "12%"))) %>% 
    ggplot(aes(x = odds.ratio, y = power, color = p0, group = p0)) +
    geom_point() +
    geom_line() +
    xlab("Detectable odds ratio") +
    ylab("Power") +
    ggsci::scale_color_nejm(name = "Base probability under H0") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()
