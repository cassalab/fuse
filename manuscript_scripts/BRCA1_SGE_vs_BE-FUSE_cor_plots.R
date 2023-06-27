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

include_LOF = T
df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_", ifelse(include_LOF, "", "noLOF_"), "042423.rds"))
ls.norm_score = readRDS("maveDB/norm_score_26genes_042423.rds")

# BRCA1 SGE data
df.BRCA1_SGE = ls.norm_score$BRCA1 

# BRCA1 baseedit data
gene_id = "BRCA1"
df.raw1 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
df.raw2 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
df = rbind(df.raw1, df.raw2)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df$gene = gene_id
df.BRCA1_BE_FUSE = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
ind = match(df.BRCA1_SGE$gene_aa_str, table = df.BRCA1_BE_FUSE$gene_aa_str)
df.BRCA1_SGE$BE_FUSE_score = df.BRCA1_BE_FUSE$final_score[ind]
df.BRCA1_SGE$BE_FUSE_score_lite = df.BRCA1_BE_FUSE$final_score_lite[ind]
df.BRCA1_SGE$BE_norm_raw_score = df.BRCA1_BE_FUSE$norm_raw_score[ind]


wb <- createWorkbook()
for (var_type in c("allVars", "misOnly")){
  df.BRCA1_SGE_2 = df.BRCA1_SGE
  if (var_type == "misOnly"){
    df.BRCA1_SGE_2 = df.BRCA1_SGE %>% filter(aaref != "*" & aaalt != "*" & aaref != aaalt) # missense only
  } 
  
  df = c()
  # BE adjusted score vs. DMS raw score
  res = cor.test(x = df.BRCA1_SGE_2$BE_FUSE_score_lite, y = df.BRCA1_SGE_2$norm_raw_score, method = "pearson")
  res2 = cor.test(x = df.BRCA1_SGE_2$BE_FUSE_score_lite, y = df.BRCA1_SGE_2$norm_raw_score, method = "spearman")
  overlapping_variants = nrow(df.BRCA1_SGE_2 %>% filter(!is.na(BE_FUSE_score_lite)))
  df = rbind(df, data.frame(name = "BE FUSE score vs. DMS raw score", overlapping_variants,
                            pearson_cor = res$estimate, pearson_pval = res$p.value,
                            spearman_cor = res2$estimate, spearman_pval = res2$p.value))
  
  # BE adjusted score (inferred) vs. DMS raw score
  res = cor.test(x = df.BRCA1_SGE_2$BE_FUSE_score, y = df.BRCA1_SGE_2$norm_raw_score, method = "pearson")
  res2 = cor.test(x = df.BRCA1_SGE_2$BE_FUSE_score, y = df.BRCA1_SGE_2$norm_raw_score, method = "spearman")
  overlapping_variants = nrow(df.BRCA1_SGE_2 %>% filter(!is.na(BE_FUSE_score)))
  df = rbind(df, data.frame(name = "BE FUSE score (with prediction) vs. DMS raw score", overlapping_variants,
                            pearson_cor = res$estimate, pearson_pval = res$p.value,
                            spearman_cor = res2$estimate, spearman_pval = res2$p.value))
  
  # BE raw score vs. DMS raw score
  res = cor.test(x = df.BRCA1_SGE_2$BE_norm_raw_score, y = df.BRCA1_SGE_2$norm_raw_score, method = "pearson")
  res2 = cor.test(x = df.BRCA1_SGE_2$BE_norm_raw_score, y = df.BRCA1_SGE_2$norm_raw_score, method = "spearman")
  overlapping_variants = nrow(df.BRCA1_SGE_2 %>% filter(!is.na(BE_norm_raw_score)))
  df = rbind(df, data.frame(name = "BE raw score vs. DMS raw score", overlapping_variants,
                            pearson_cor = res$estimate, pearson_pval = res$p.value,
                            spearman_cor = res2$estimate, spearman_pval = res2$p.value))
  
  addWorksheet(wb, var_type)
  writeData(wb, var_type, df)
}

saveWorkbook(wb, file = "~/Google Drive/My Drive/manuscript_resub/sup_tables/BRCA1_SGE_BE-FUSE_cor.xlsx", overwrite = TRUE)


