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

# make FUNSUM from 26 gene's DMS
ls.norm_score = readRDS("maveDB/norm_score_26genes_042423.rds")
gene_ids.exclude = c("LDLR", "AICDA")

include_LOF = T
df.all_sub_score = c()
ls.funsum = list()
for (gene in names(ls.norm_score)){
  df = ls.norm_score[[gene]]
  df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF, score_type = "norm_raw_score", sd_norm = F)
  df.funsum = get_FUNSUM(df.sub_score, avg_method = "js", include_LOF = include_LOF)
  
  if (!(gene %in% gene_ids.exclude)){
    ls.funsum[[gene]] = df.funsum
    df.all_sub_score = rbind(df.all_sub_score, cbind(gene, df.sub_score))
  }
}
df.funsum_all = get_FUNSUM(df.all_sub_score[,-1], avg_method = "js", include_LOF = include_LOF ) # get combined funsum
df.sub_tb = funsum_to_subTable(df.funsum_all) #convert FUNSUM to tabular format

saveRDS(df.funsum_all, file = "./shiny_app/DMS_denoiser/funsum_maveDB_042423.rds")