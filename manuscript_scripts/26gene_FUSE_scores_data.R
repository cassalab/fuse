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

# de-noise training DMS data for testing base editing data
include_LOF = T
df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_", ifelse(include_LOF, "", "noLOF_"), "042423.rds"))
gene_ids.exclude = c("LDLR", "AICDA")

ls.norm_score = readRDS("maveDB/norm_score_26genes_042423.rds")
all_gene_ids = names(ls.norm_score)
df.out_all = c()
for (gene_id in all_gene_ids){
  df = ls.norm_score[[gene_id]]
  df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
  df.out_all = rbind(df.out_all, df.out)
}
df.out_all = df.out_all %>% filter(!(gene %in% gene_ids.exclude))
write.xlsx(df.out_all, "~/Google Drive/My Drive/manuscript_resub/sup_tables/Table_S1_FUSE_estimated_24genes.xlsx")

# write_csv(df.out_all, file = "./maveDB/FUSE_estimated_26genes.csv")


# number of assayed and predicted variants in the selected 9 genes
selected_gene_ids = c("TP53", "PTEN", "HRAS", "MSH2", "GCK", "MTHFR", "CBS", "HMBS", "SNCA")

temp = df.out_all %>% filter(gene %in% selected_gene_ids)
nrow(temp %>% filter(!is.na(norm_raw_score)))
nrow(temp %>% filter(!is.na(final_score)))

# number of assayed and predicted variants in 26 genes
nrow(df.out_all %>% filter(!is.na(norm_raw_score)))
nrow(df.out_all %>% filter(!is.na(final_score)))

