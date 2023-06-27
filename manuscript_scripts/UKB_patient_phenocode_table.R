setwd("~/OneDrive/Tian_Workfile/BWH/AA_matrix_proj/")
library(openxlsx)
library(tidyverse)
library(pROC)
library(patchwork)
library(ptm)
library(corrplot)
library(Biostrings)
library(hrbrthemes)
source("AA_proj_functions_101421.R")

# supplement table for patient phenocodes
LFS_cancer_types # from TP53 patient LFS analysis (includes breastcancer for BRCA1 analysis)
cancer_types = gsub(pattern = "__", replacement = "/", x = LFS_cancer_types)
cancer_types = gsub(pattern = "_", replacement = " ", x = cancer_types)
cancer_types[cancer_types == 'sarcomafibrosarcoma'] = 'sarcoma/fibrosarcoma'
cancer_types[cancer_types == 'kidneyrenal cell cancer'] = 'kidney/renal cell cancer'
cancer_types[cancer_types == 'colon cancersigmoid cancer'] = 'colon cancer/sigmoid cancer'

df.patient_mapping = read_tsv("UKBB/icd10_mapping.tsv")
df.patient_mapping = df.patient_mapping %>% 
  mutate('Cancer types' = gsub(pattern=" / ", replacement="/", x=meaning))
ind = match(x = cancer_types, table = df.patient_mapping$`Cancer types`)
df.patient_mapping = df.patient_mapping[ind,] %>% 
  mutate('UKB self-reported codes' = `self-reported coding`, 'ICD-10 codes' = `ICD-10 Codes`) %>%
  select('Cancer types', 'UKB self-reported codes', 'ICD-10 codes')
write.xlsx(df.patient_mapping, file = "~/Google Drive/My Drive/manuscript_resub/sup_tables/patient_phenocodes.xlsx")
