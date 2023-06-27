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

# normalize original DMS data and selected maveDB scoresets
df.clinvar = read.table("clinvar_lite_GRCh38_081722.txt", header = T, sep = "\t")
df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
ls.norm_score = list()
ls.DMS_score = readRDS("DMS_score_14genes.rds")

for (gene_id in names(ls.DMS_score)){
  df = ls.DMS_score[[gene_id]]
  df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9)
  df.clinvar_gene = df.clinvar %>% filter(gene == gene_id)
  pos_diff = clinvar_align_aapos(df, df.clinvar_gene, diagnose = T)
  if (pos_diff != 0){
    print(paste0(gene_id, " pos_diff=", pos_diff))
  }
  df = df %>% 
    mutate(gene, aapos = aapos+pos_diff, gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt)) %>%
    select(gene, aapos, aaref, aaalt, gene_aa_str, functional_class, raw_score, norm_raw_score)
  ls.norm_score[[gene_id]] = df
  # plt = draw_score_distribution_ggplot(df, score_type = "norm_raw_score")
  # ggsave(filename = paste0("./maveDB/plots/norm_raw_score/original_DMS/", gene, ".png"), plot = plt, device = "png", width = 8, height = 5)
}

selected_scoreset_ids = c("00000050-a-1", "00000096-a-1", "00000049-a-4", "00000005-a-2", "00000108-a-3", "00000057-a-1", 
                          "00000048-a-1", "00000045-h-1", "00000066-a-1", "00000095-b-1", "00000047-a-1", "00000106-b-1",
                          "00000002-a-2")
df.maveDB_4 = read_csv(file = "./maveDB/maveDB_extended_scoreset_with_clinvar.csv")
for (scoreset_id in selected_scoreset_ids){
  ind = which(df.maveDB_4$scoreset_id == scoreset_id)
  gene_id = df.maveDB_4$gene_symbol[ind]
  df = ls.scoreset[[scoreset_id]]
  df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9)
  df.clinvar_gene = df.clinvar %>% filter(gene == gene_id)
  pos_diff = clinvar_align_aapos(df, df.clinvar_gene, diagnose = T)
  if (pos_diff != 0){
    print(paste0(gene_id, " pos_diff=", pos_diff))
  }
  df = df %>% 
    mutate(gene, aapos = aapos+pos_diff, gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt)) %>%
    select(gene, aapos, aaref, aaalt, gene_aa_str, functional_class, raw_score, norm_raw_score)
  ls.norm_score[[gene_id]] = df
  # plt = draw_score_distribution_ggplot(df, score_type = "norm_raw_score")
  # ggsave(filename = paste0("./maveDB/plots/norm_raw_score/maveDB_selected_DMS/", gene, ".png"), plot = plt, device = "png", width = 8, height = 5)
}

saveRDS(ls.norm_score, file = "./maveDB/norm_score_26genes_042423.rds")