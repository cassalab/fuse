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

# barplot to show number of assayed and predicted sites in functional datasets
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_052024.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")


# estimate FUSE score for all human datasets from maveDB and proteinGym
df.out_all = c()
for (i in 1:nrow(df.combined_scoreset_info)){
  scoreset_id = df.combined_scoreset_info$scoreset_id[i]
  df = ls.combined_scoreset[[scoreset_id]]
  df$gene = df.combined_scoreset_info$gene_symbol[i]
  df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                             dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[i], ".dss"), 
                             include_LOF = include_LOF, show_func_class = T)
  df.out = df.out %>% mutate(scoreset_id, source = df.combined_scoreset_info$source[i], 
                             assay_type = ifelse(gene == "BRCA1", "SGE", "DMS"))
  df.out$source = df.combined_scoreset_info$source[i]
  df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
  df.out_all = rbind(df.out_all, df.out)
}


# add DMS datasets
df.obs = data.frame(facet = "DMS", gene = unique(df.out_all$gene), Predicted=NA, Assayed=NA)
# df.obs = df.obs %>% filter(gene != "LDLR" & gene != "AID") # remove LDLR and AID
df.obs$facet[df.obs$gene == "BRCA1"] = "SGE"

df.out_all2 = df.out_all %>% group_by(gene, gene_aa_str) %>% 
  summarise(norm_raw_score = mean(norm_raw_score), final_score_ss = mean(final_score_ss))
for (i in 1:nrow(df.obs)){
  df = df.out_all2 %>% filter(gene == df.obs$gene[i])
  df1 = df %>% filter(!is.na(norm_raw_score))
  df.obs$Predicted[i] = nrow(df) - nrow(df1)
  df.obs$Assayed[i] = nrow(df1)
}
df.obs = df.obs %>% arrange(facet, gene)

# add base editng datasets
gene_id = "BRCA1"
df.raw1 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
df.raw2 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
df = rbind(df.raw1, df.raw2)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df$gene = gene_id
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)

df.out = df.out %>% filter(functional_class == "MIS")
temp = data.frame(facet = "Base editing", gene = gene_id, 
                  Predicted = length(which(is.na(df.out$norm_raw_score))), 
                  Assayed = length(which(!is.na(df.out$norm_raw_score))))
df.obs = rbind(df.obs, temp)

gene_id = "LDLR"
df.raw = read_csv("../Jayoung_proj/functional_score/TableSX_tiling_screen_UKBmerged.20230817.alpha_fixed.csv")
df.raw = df.raw %>% filter(coding == "coding" & !grepl("control", pos, ignore.case=T)) %>% 
  mutate(aaref = ref, aaalt = alt, aapos = as.numeric(gsub(pattern = "A", replacement = "", x = pos)), index=`0`) %>% 
  select(index, aapos, aaref, aaalt, contains("_z_"), mu_sd_adj)

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

df.out = df.out %>% filter(functional_class == "MIS")
temp = data.frame(facet = "Base editing", gene = gene_id, 
                  Predicted = length(which(is.na(df.out$norm_raw_score))), 
                  Assayed = length(which(!is.na(df.out$norm_raw_score))))
df.obs = rbind(df.obs, temp)



# df.raw = read.csv("./shiny_app/DMS_denoiser/DMS_data/DDR_genes_BE.csv", header = T)
# df = normalize_scoreset(df.raw, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
# df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
# temp = data.frame(facet = "Base editing", gene = "DDR genes", 
#                   Predicted = length(which(is.na(df.out$norm_raw_score))), 
#                   Assayed = length(which(!is.na(df.out$norm_raw_score))))
# df.obs = rbind(df.obs, temp)
# df.obs$facet = factor(df.obs$facet, levels = c("DMS", "Base editing", "SGE"))

df.obs$facet = factor(df.obs$facet, levels = c("DMS", "SGE", "Base editing"))
df.obs = df.obs %>% arrange(facet, gene)
df.obs_A = df.obs[1:27,] 
df.obs_B = df.obs[28:50,]

df.obs_A_long = df.obs_A %>% pivot_longer(cols = c("Predicted", "Assayed"), names_to = "Score_type", values_to = "Variant_count")
df.obs_A_long$Score_type = factor(df.obs_A_long$Score_type, levels = c("Predicted", "Assayed"))
df.obs_A_long$facet = factor(df.obs_A_long$facet, levels = c("DMS", "SGE", "Base editing"))
df.obs_B_long = df.obs_B %>% pivot_longer(cols = c("Predicted", "Assayed"), names_to = "Score_type", values_to = "Variant_count")
df.obs_B_long$Score_type = factor(df.obs_B_long$Score_type, levels = c("Predicted", "Assayed"))
df.obs_B_long$facet = factor(df.obs_B_long$facet, levels = c("DMS", "SGE", "Base editing"))

plt_A <- ggplot(df.obs_A_long, aes(fill=Score_type, y=gene, x=Variant_count)) + 
  geom_col(position="stack", show.legend = T, color="black") + xlim(0, 40000) + 
  scale_fill_manual(values = c("#FD8972", "#8EDAF5"), name = "") + 
  ggtitle("Variants with FUSE estimated functional score\n") + xlab(label = "") + ylab(label = "Count")  + 
  facet_grid(cols = vars(facet), scales='free', space = "free", switch = "y") +
  guides(fill="none") +
  theme_ipsum(base_family = "Arial", plot_title_size = 15, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust=0.5, vjust=2), 
        # strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

plt_B <- ggplot(df.obs_B_long, aes(fill=Score_type, y=gene, x=Variant_count)) + 
  geom_col(position="stack", show.legend = T, color="black") + xlim(0, 40000) + 
  scale_fill_manual(values = c("#FD8972", "#8EDAF5"), name = "") + 
  ggtitle("\n") + xlab(label = "") + ylab(label = "Count")  + 
  facet_grid(facet~., scales='free', space = "free", switch = "y") +
  theme_ipsum(base_family = "Arial", plot_title_size = 15, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust=0.5, vjust=2), 
        legend.position = c(0.85, 0.95), 
        legend.background = element_rect(fill="white", linetype=0),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plt = plt_A | plt_B
ggsave(paste0("./manuscript/new_figures/variant_count.png"), plt, device = "png", width = 10, height = 10) 



# generate datatables for all functional scoresets
# estimate FUSE score for all human datasets from maveDB and proteinGym
df.out_all = c()
for (i in 1:nrow(df.combined_scoreset_info)){
  scoreset_id = df.combined_scoreset_info$scoreset_id[i]
  df = ls.combined_scoreset[[scoreset_id]]
  df$gene = df.combined_scoreset_info$gene_symbol[i]
  df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                             dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[i], ".dss"), 
                             include_LOF = include_LOF, show_func_class = T)
  df.out = df.out %>% mutate(scoreset_id, source = df.combined_scoreset_info$source[i], 
                             assay_type = ifelse(gene == "BRCA1", "SGE", "DMS"))
  df.out$source = df.combined_scoreset_info$source[i]
  df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
  df.out_all = rbind(df.out_all, df.out)
}

gene_id = "BRCA1"
df.raw1 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
df.raw2 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
df = rbind(df.raw1, df.raw2)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df$gene = gene_id
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                          dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                          include_LOF = include_LOF, show_func_class = T)
df.out = df.out %>% filter(functional_class == "MIS")
df.out$scoreset_id = NA
df.out$source = "Individual studies"
df.out$assay_type = "Base editing"
df.final = rbind(df.out_all, df.out)
df.final = df.final %>% select(gene, scoreset_id, source, assay_type, aapos, aaref, aaalt, functional_class, 
                               DSSP_ss = ss, DSSP_acc = acc, raw_score, norm_raw_score, pos_score, sub_score_ss, FUSE_score = final_score_ss)

gene_id = "LDLR"
df.raw = read_csv("../Jayoung_proj/functional_score/TableSX_tiling_screen_UKBmerged.20230817.alpha_fixed.csv")
df.raw = df.raw %>% filter(coding == "coding" & !grepl("control", pos, ignore.case=T)) %>% 
  mutate(aaref = ref, aaalt = alt, aapos = as.numeric(gsub(pattern = "A", replacement = "", x = pos)), index=`0`) %>% 
  select(index, aapos, aaref, aaalt, contains("_z_"), mu_sd_adj)

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
df.out = df.out %>% filter(functional_class == "MIS")
df.out$scoreset_id = NA
df.out$source = "Individual studies"
df.out$assay_type = "Base editing"
df.final = rbind(df.out_all, df.out)
df.final = df.final %>% select(gene, scoreset_id, source, assay_type, aapos, aaref, aaalt, functional_class, 
                               DSSP_ss = ss, DSSP_acc = acc, raw_score, norm_raw_score, pos_score, sub_score_ss, FUSE_score = final_score_ss)

write_csv(df.final, "./manuscript/new_figures/supplementary_table/Table_S1__all_FUSE_scores.csv")

