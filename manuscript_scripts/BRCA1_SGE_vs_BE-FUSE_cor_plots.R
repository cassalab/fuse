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

include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")
df.clinvar = read.table("./clinvar/clinvar_lite_GRCh38_091323.txt", header = T, sep = "\t")

# BRCA1 SGE data
gene_id = "BRCA1"
df = ls.combined_scoreset[[gene_id]]
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df$gene = df.combined_scoreset_info$gene_symbol[ind]
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df.out$SGE_norm_raw_score = df.out$norm_raw_score
df.out$SGE_FUSE_score = df.out$final_score_ss
df.out$SGE_FUSE_score_lite = df.out$final_score_ss
df.out$SGE_FUSE_score_lite[is.na(df.out$SGE_norm_raw_score)] = NA
df.out = df.out %>% filter(aaref != aaalt & aaref != "*" & aaalt != "*") %>% 
  select(gene_aa_str, SGE_norm_raw_score, SGE_FUSE_score, SGE_FUSE_score_lite)

# BRCA1 baseedit data
gene_id = "BRCA1"
df.raw1 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
df.raw2 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
df = rbind(df.raw1, df.raw2)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df$gene = gene_id
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df.out2= de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df.out2$BE_norm_raw_score = df.out2$norm_raw_score
df.out2$BE_FUSE_score = df.out2$final_score_ss
df.out2$BE_FUSE_score_lite = df.out2$final_score_ss
df.out2$BE_FUSE_score_lite[is.na(df.out2$BE_norm_raw_score)] = NA
df.out2 = df.out2 %>% filter(aaref != aaalt & aaref != "*" & aaalt != "*") %>%
  select(gene_aa_str, BE_norm_raw_score, BE_FUSE_score, BE_FUSE_score_lite)

df.out3 = merge(df.out, df.out2, all=T, by = c("gene", "gene_aa_str"))



df = c()
# BE adjusted score vs. DMS raw score
res = cor.test(x = df.out3$BE_FUSE_score_lite, y = df.out3$SGE_norm_raw_score, method = "pearson")
res2 = cor.test(x = df.out3$BE_FUSE_score_lite, y = df.out3$SGE_norm_raw_score, method = "spearman")
overlapping_variants = nrow(df.out3 %>% filter(!is.na(BE_FUSE_score_lite) & !is.na(SGE_norm_raw_score)))
df = rbind(df, data.frame(name = "BE FUSE score vs. DMS raw score", overlapping_variants,
                          pearson_cor = res$estimate, pearson_pval = res$p.value,
                          spearman_cor = res2$estimate, spearman_pval = res2$p.value))

# BE adjusted score (inferred) vs. DMS raw score
res = cor.test(x = df.out3$BE_FUSE_score, y = df.out3$SGE_norm_raw_score, method = "pearson")
res2 = cor.test(x = df.out3$BE_FUSE_score, y = df.out3$SGE_norm_raw_score, method = "spearman")
overlapping_variants = nrow(df.out3 %>% filter(!is.na(BE_FUSE_score) & !is.na(SGE_norm_raw_score)))
df = rbind(df, data.frame(name = "BE FUSE score (with prediction) vs. DMS raw score", overlapping_variants,
                          pearson_cor = res$estimate, pearson_pval = res$p.value,
                          spearman_cor = res2$estimate, spearman_pval = res2$p.value))

# BE raw score vs. DMS raw score
res = cor.test(x = df.out3$BE_norm_raw_score, y = df.out3$SGE_norm_raw_score, method = "pearson")
res2 = cor.test(x = df.out3$BE_norm_raw_score, y = df.out3$SGE_norm_raw_score, method = "spearman")
overlapping_variants = nrow(df.out3 %>% filter(!is.na(BE_norm_raw_score) & !is.na(SGE_norm_raw_score)))
df = rbind(df, data.frame(name = "BE raw score vs. DMS raw score", overlapping_variants,
                          pearson_cor = res$estimate, pearson_pval = res$p.value,
                          spearman_cor = res2$estimate, spearman_pval = res2$p.value))

write_csv(df, file = "./manuscript/new_figures/BRCA1_SGE_BE-FUSE_cor_2.csv")

