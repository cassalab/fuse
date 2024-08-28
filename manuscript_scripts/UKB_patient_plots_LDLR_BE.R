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

# # extract am scores for LDLR
# df.am_LDLR = read_tsv("~/Downloads/AlphaMissense_aa_substitutions.tsv.gz", skip = 3, lazy = T) %>% filter(uniprot_id == "P01130")
# df.am_LDLR$gene = "LDLR"
# df.am_LDLR = df.am_LDLR %>% mutate(gene_aa_str = paste0(gene, "---", protein_variant))
# write_csv(df.am_LDLR, "./ProteinGym/LDLR_am_scores.csv")


# LDLR base editing data
# load phenotype info
df.patient = read.csv("./UKBB/updated_patient_to_phenotype_042822.csv")
df.cad = read_csv("../Jayoung_proj/cad.csv")

# load LDLR v2p and convert to p2v
df.v2p = read_csv("../Jayoung_proj/variant_list.csv") %>% filter(gene == "LDLR")
temp = str_split(string = df.v2p$variant, pattern = "-", simplify = T)
ind = which((nchar(temp[,3]) == 1) & (nchar(temp[,4]) == 1))
df.v2p = df.v2p[ind,]
df.p2v = v2p_to_p2v(df.v2p)
df.p2v = df.p2v %>% filter(!grepl(pattern="|", fixed = T, x = variant))

# load dfNSFP to convert from nucleotide to amino acid
df.NSFP = read.xlsx("../AA_matrix_proj/dbNSFP_BRCA1_LDLR_TP53_includeLOF.xlsx") %>% filter(genename == "LDLR") %>%
  mutate(variant_str = paste(chr, pos, ref, alt, sep = "-"), gene_aa_str = paste0(genename, "---", aaref, aapos, aaalt))

# merge cad status, attach to df.p2v
df.p2v$cad = 0
ind = match(df.p2v$patient, table = df.cad$eid)
df.p2v$cad = df.cad$disease[ind]
df.p2v$cad = factor(df.p2v$cad)
df.p2v$cad_group = "Yes"
df.p2v$cad_group[df.p2v$cad == 0] = "No"

# attach other info to df.p2v
ind = match(df.p2v$patient, table = df.patient$new_id)
df.p2v$ldl = df.patient$ldl[ind]
df.p2v$normalized_cad_PRS = df.patient$normalized_cad_PRS[ind]
ind = match(df.p2v$variant, table = df.NSFP$variant_str)
df.p2v$gene_aa_str = df.NSFP$gene_aa_str[ind]
df.p2v$aapos = df.NSFP$aapos[ind]
df.p2v$aaref = df.NSFP$aaref[ind]
df.p2v$aaalt = df.NSFP$aaalt[ind]

## attach functional scores
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")

# load LDLR Bean scores
df.raw = read_csv("../Jayoung_proj/functional_score/TableSX_tiling_screen_UKBmerged.20230817.alpha_fixed.csv")
df.raw = df.raw %>% filter(coding == "coding" & !grepl("control", pos, ignore.case=T)) %>% 
  mutate(aaref = ref, aaalt = alt, aapos = as.numeric(gsub(pattern = "A", replacement = "", x = pos)), index=`0`) %>% 
  select(index, aapos, aaref, aaalt, contains("_z_"), mu_sd_adj)

# calculate FUSE score
gene_id = "LDLR"
df = df.raw %>% 
  mutate(gene=gene_id, raw_score = mu_z_adj, gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt), mu_sd_adj) %>% 
  select(gene, aapos, aaref, aaalt, gene_aa_str, raw_score, mu_sd_adj) %>%
  filter(aapos>0)
df$functional_class = NA
ind = which(df$aaref != df$aaalt & df$aaalt != '*')
df$functional_class[ind] = "MIS"
ind = which(df$aaref != "*" & df$aaalt == "*")
df$functional_class[ind] = "LOF"
ind = which(df$aaref == df$aaalt)
df$functional_class[ind] = "SYN"

df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", gene_id, ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df.out$final_score_ss_lite = df.out$final_score_ss
df.out$final_score_ss_lite[is.na(df.out$raw_score)] = NA

# attach am score
df.am_LDLR = read_csv("./ProteinGym/LDLR_am_scores.csv")
ind = match(df.out$gene_aa_str, table = df.am_LDLR$gene_aa_str)
df.out$am_score = df.am$am_pathogenicity[ind]

## draw plots
df.p2v_2 = merge(df.p2v, df.out)

# score distribution between patient groups
ymax = 2.5
bw = 0.1
x_range = c(-2, 2.5)
plt1 = draw_patient_ggplot_v2(df.p2v_2, score_type = "norm_raw_score", score_type_name = "Normalized original score", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.22, 0.85))
plt2 = draw_patient_ggplot_v2(df.p2v_2, score_type = "final_score_ss_lite", score_type_name = "FUSE score\n(without predicted sites)", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot_v2(df.p2v_2, score_type = "final_score_ss", score_type_name = "FUSE score\n(with predicted sites)", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw)

plt = (plt1|plt2|plt3) + plot_annotation(tag_levels = 'A')
ggsave(paste0("./manuscript/new_figures/UKBB_plots/UKBB_allPlots_LDLR_BE.png"), plt, device = "png", width = 12, height = 4)








