setwd("/dcs05/ladd/NDEpi/data/Projects/InProgress/zhe/SEED_teen_phenotype_LOI")

require(tidyverse)
require(gtsummary)


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
SEED4_fam <- read_table(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/SNPArray_DataFiles/CleanMeasuredGenotypes/cleaned2020", "wave4.fam"),
                        col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"))
SEED5_fam <- read_table(file.path("/dcs05/ladd/NDEpi/data/MasterCohortData/SEED/SNPArray_DataFiles/CleanMeasuredGenotypes/cleaned2020", "wave5.fam"),
                        col_names = c("FID", "IID", "PID", "MID", "sex", "phenotype"))

# Overlapping individuals
intersect(SEED_teen_hds_1a[SEED_teen_hds_1a$HS_A9AA_ADHD != ".",]$FamilyID, SEED1_fam$FID) %>% length()   # 258
intersect(SEED_teen_hds_1a[SEED_teen_hds_1a$HS_A9AA_ADHD != ".",]$FamilyID, SEED2_fam$FID) %>% length()   # 40
intersect(SEED_teen_hds_1a[SEED_teen_hds_1a$HS_A9AA_ADHD != ".",]$FamilyID, SEED3_fam$FID) %>% length()   # 94
intersect(SEED_teen_hds_1a[SEED_teen_hds_1a$HS_A9AA_ADHD != ".",]$FamilyID, SEED4_fam$FID) %>% length()   # 0
intersect(SEED_teen_hds_1a[SEED_teen_hds_1a$HS_A9AA_ADHD != ".",]$FamilyID, SEED5_fam$FID) %>% length()   # 0

# Merge and limit to individuals with genetic data
SEED_teen <- SEED_teen_participant %>% 
    left_join(SEED_teen_hds_1a, by = "FamilyID") %>% 
    left_join(SEED_teen_hds_2a, by = "FamilyID") %>% 
    left_join(SEED_teen_hds_3a, by = "FamilyID") %>% 
    filter(FamilyID %in% c(SEED1_fam$FID, SEED2_fam$FID, SEED3_fam$FID))


######################
# Summary statistics #
######################
# Psychiatric phenotypes
table(SEED_teen$HS_A9AA_ADHD)   # ADHD
table(SEED_teen$HS_A9FA_ANX)    # Anxiety
table(SEED_teen$HS_A9IA_AUT)    # ASD
table(SEED_teen$HS_A9KA_BPD)    # Bipolar disorder
table(SEED_teen$HS_A9RA_DEP)    # Depression
table(SEED_teen$HS_A9VA_EAT)    # Eating disorder
table(SEED_teen$HS_A9ACA_OCD)   # OCD
table(SEED_teen$HS_A9AHA_TOU)   # Tourette syndrome

# Neurological phenotypes
table(SEED_teen$HS_A9LA_BRAIN)  # Brain injury, concussion or head injury
table(SEED_teen$HS_A9OA_CP)     # Cerebral palsy
table(SEED_teen$HS_A9WA_SEZ)    # Epilepsy or seizure disorder
table(SEED_teen$HS_A9YA_MIG)    # Frequent or severe headaches, including migraine
table(SEED_teen$HS_A9SA_DEV)    # Developmental delay
table(SEED_teen$HS_A9AAA_MR)    # Intellectual disability
table(SEED_teen$HS_A9AMA_SPE)   # Speech or other language disorder

# Behavioral-cognitive phenotypes
table(SEED_teen$HS_A1_HT)    # Height (measured class)
SEED_teen <- SEED_teen %>% 
    mutate(height = case_when(HS_A1_HT %in% c("-1", ".", "99") ~ as.numeric(NA),
                              HS_A1_HT == 1 ~ as.numeric(HS_A1_HTTAPE),
                              HS_A1_HT == 2 ~ as.numeric(HS_A1_HTRECALL)))
summary(SEED_teen$height)

table(SEED_teen$HS_A2_WT)    # Weight (measured class)
SEED_teen <- SEED_teen %>% 
    mutate(weight = case_when(HS_A2_WT %in% c("-1", ".", "99") ~ as.numeric(NA),
                              HS_A2_WT == 1 ~ as.numeric(HS_A2_WTHOME),
                              HS_A2_WT == 2 ~ as.numeric(HS_A2_WTRECALL)))
summary(SEED_teen$weight)

SEED_teen <- SEED_teen %>%   # BMI
    mutate(BMI = (weight/height^2) * 703)   # Formula: weight (lb) / [height (in)]2 x 703
SEED_teen <- SEED_teen %>% 
    mutate(BMI_cat = case_when(BMI < 18.5 ~ 1,                 # Underweight
                               BMI >= 18.5 & BMI < 25.0 ~ 2,   # Normal weight
                               BMI >= 25.0 & BMI < 30.0 ~ 3,   # Overweight
                               BMI >= 30.0 ~ 4))               # Obesity

SEED_teen <- SEED_teen %>%  # Subjective well-being
    mutate(HS_A3_HEALTH = ifelse(HS_A3_HEALTH == "-1"|is.na(HS_A3_HEALTH), as.numeric(NA), as.numeric(HS_A3_HEALTH)))
table(SEED_teen$HS_A3_HEALTH)

SEED_teen <- SEED_teen %>%  # Verbal communications
    mutate(HS_A7_VERBAL = ifelse(HS_A7_VERBAL == "-1"|HS_A7_VERBAL == "."|is.na(HS_A7_VERBAL), as.numeric(NA), as.numeric(HS_A7_VERBAL))) %>% 
    mutate(HS_A7_VERBAL = case_when(is.na(HS_A7_VERBAL) ~ as.numeric(NA),
                                    HS_A7_VERBAL == 1 ~ 0,
                                    HS_A7_VERBAL %in% 2:5 ~ 1))
table(SEED_teen$HS_A7_VERBAL)

SEED_teen <- SEED_teen %>%  # Sleep hours
    mutate(HS_A17_USUHRS = ifelse(HS_A17_USUHRS == "-1"|HS_A17_USUHRS == "."|is.na(HS_A17_USUHRS), as.numeric(NA), as.numeric(HS_A17_USUHRS)))
table(SEED_teen$HS_A17_USUHRS)

table(SEED_teen$HS_A9JA_BEH)   # Behavioral or conduct problems

table(SEED_teen$HS_A9ADA_INJ)  # Self-injurious behavior

table(SEED_teen$HS_A9AEA_SENS) # Sensory integration disorder


###################
# Chi-square test #
###################
SEED_teen_pheno <- SEED_teen %>% 
    select(FamilyID, HS_A9AA_ADHD, HS_A9FA_ANX, HS_A9IA_AUT, HS_A9KA_BPD, HS_A9RA_DEP, HS_A9VA_EAT, HS_A9ACA_OCD, HS_A9AHA_TOU,
           HS_A9LA_BRAIN, HS_A9OA_CP, HS_A9WA_SEZ, HS_A9YA_MIG, HS_A9SA_DEV, HS_A9AAA_MR, HS_A9AMA_SPE, 
           height, BMI_cat, HS_A3_HEALTH, HS_A7_VERBAL, HS_A17_USUHRS, HS_A9JA_BEH, HS_A9ADA_INJ, HS_A9AEA_SENS)

SEED_teen_pheno <- SEED_teen_pheno %>% 
    filter(if_any(c(HS_A9AA_ADHD:HS_A9AMA_SPE, HS_A9JA_BEH:HS_A9AEA_SENS), ~ .x != ".")) %>% 
    mutate(across(c(HS_A9AA_ADHD:HS_A9AMA_SPE, HS_A9JA_BEH:HS_A9AEA_SENS), 
                  ~ ifelse(.x == "-1"|.x == "99", as.numeric(NA), as.numeric(.x))))

SEED_teen_pheno <- SEED_teen_pheno %>% 
    rename("ADHD" = "HS_A9AA_ADHD", "Anxiety" = "HS_A9FA_ANX", "ASD" = "HS_A9IA_AUT", "Bipolar disorder" = "HS_A9KA_BPD",
           "Depression" = "HS_A9RA_DEP", "Eating disorder" = "HS_A9VA_EAT", "OCD" = "HS_A9ACA_OCD", "Tourette syndrome" = "HS_A9AHA_TOU",
           "Brain injury" = "HS_A9LA_BRAIN", "Cerebral palsy" = "HS_A9OA_CP", "Seizure disorder" = "HS_A9WA_SEZ",
           "Migraine" = "HS_A9YA_MIG", "Developmental delay" = "HS_A9SA_DEV", "Intellectual disability" = "HS_A9AAA_MR", 
           "Language disorder" = "HS_A9AMA_SPE", "Height" = "height", "BMI" = "BMI_cat", "Subjective well-being" = "HS_A3_HEALTH",
           "Verbal communication" = "HS_A7_VERBAL", "Sleep hours" = "HS_A17_USUHRS", "Behavioral problems" = "HS_A9JA_BEH",
           "Self-injurious" = "HS_A9ADA_INJ", "Sensory integration" = "HS_A9AEA_SENS")

require(lsr)
# function to get chi square p value and Cramers V
f = function(x,y) {
    tbl = SEED_teen_pheno %>% select(x,y) %>% table()
    chisq_pval = round(chisq.test(tbl)$p.value, 3)
    cramV = round(cramersV(tbl), 3) 
    data.frame(x, y, chisq_pval, cramV) }

### Psychiatric phenotypes
df_comb = data.frame(t(combn(sort(names(SEED_teen_pheno %>% select(ADHD:`Tourette syndrome`))), 2)), stringsAsFactors = F)

df_res = map2_df(df_comb$X1, df_comb$X2, f)
df_res %>%
    ggplot(aes(x, y, fill = chisq_pval)) +
    geom_tile() +
    geom_text(aes(x, y, label = cramV)) +
    scale_fill_gradient("Chi-square p-value", low = "#f3f4f7", high = "#0099e5", na.value = "white") +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### Neurological phenotypes
df_comb = data.frame(t(combn(sort(names(SEED_teen_pheno %>% select(`Brain injury`:`Language disorder`))), 2)), stringsAsFactors = F)

df_res = map2_df(df_comb$X1, df_comb$X2, f)
df_res %>%
    ggplot(aes(x, y, fill = chisq_pval)) +
    geom_tile() +
    geom_text(aes(x, y, label = cramV)) +
    scale_fill_gradient("Chi-square p-value", low = "#f3f4f7", high = "#0099e5", na.value = "white") +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### Behavioral-cognitive phenotypes
df_comb = data.frame(t(combn(sort(names(SEED_teen_pheno %>% select(Height:`Sensory integration`))), 2)), stringsAsFactors = F)

df_res = map2_df(df_comb$X1, df_comb$X2, f)
df_res %>%
    ggplot(aes(x, y, fill = chisq_pval)) +
    geom_tile() +
    geom_text(aes(x, y, label = cramV)) +
    scale_fill_gradient("Chi-square p-value", low = "#f3f4f7", high = "#0099e5", na.value = "white") +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### Selected phenotypes
df_comb = data.frame(t(combn(sort(names(SEED_teen_pheno %>% 
                                            select(ADHD, Anxiety, `Developmental delay`, `Language disorder`, `Sensory integration`,
                                                   BMI, `Verbal communication`))), 2)), stringsAsFactors = F)

df_res = map2_df(df_comb$X1, df_comb$X2, f)
df_res %>%
    ggplot(aes(x, y, fill = chisq_pval)) +
    geom_tile() +
    geom_text(aes(x, y, label = cramV)) +
    scale_fill_gradient("Chi-square p-value", low = "#f3f4f7", high = "#0099e5", na.value = "white") +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


####################################
# Stratify by final classification #
####################################
SEED_FC <- read.csv("/dcs05/ladd/NDEpi/data/Projects/SEED/LOI/44/output/phenoFiles/moms_kids.allData.all_pheno_PC_data.csv") %>% 
    select(FID, DR_FC)
SEED_teen_pheno <- SEED_teen_pheno %>% 
    mutate(FamilyID = as.character(FamilyID)) %>% 
    left_join(SEED_FC %>% distinct(FID, .keep_all = TRUE) %>% rename("FamilyID" = "FID"), by = "FamilyID")

# All phenotypes
SEED_teen_pheno %>% 
    select(- FamilyID) %>% 
    tbl_summary(by = "DR_FC", missing = "no") %>% 
    add_overall()

# Selected phenotypes
SEED_teen_pheno %>% 
    select(DR_FC, ADHD, Anxiety, `Developmental delay`, `Language disorder`, `Sensory integration`, BMI, `Verbal communication`) %>% 
    tbl_summary(by = "DR_FC", missing = "no") %>% 
    add_overall()
