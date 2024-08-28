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


# compare FUSE with Imputation Tool on selected genes
df.NSFP = read_csv("./weile_tool/weile_tool_genes_dbNSFP_parsed.csv") 
df.NSFP$aaref[df.NSFP$aaref == "X"] = "*"
df.NSFP$aaalt[df.NSFP$aaalt == "X"] = "*"
df.NSFP = df.NSFP %>% 
  mutate(gene_aa_str = paste0(genename, "---", aaref, aapos, aaalt)) %>%
  mutate(gene_aapos = paste0(genename, "---", aapos))

include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")
df.clinvar = read.table("./clinvar/clinvar_lite_GRCh38_091323.txt", header = T, sep = "\t")


# de-noise functional scores
gene_id = "BRCA1"
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df = read_csv("weile_tool/_P38398[BRCA1_BRCT]_imputation.csv") %>% 
  mutate(gene = "BRCA1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9)
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df = merge(df, df.out, all=T)

gene_id = "BRCA1"
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df2 = read_csv("weile_tool/_P38398[BRCA1_RING]_imputation.csv") %>% 
  mutate(gene = "BRCA1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df2 = normalize_scoreset(df2, lower_bound = 0.1, upper_bound = 0.9)
df.out = de_noise_ss_1gene(df = df2, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df2 = merge(df2, df.out, all=T)

gene_id = "CALM1"
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df3 = read_csv("weile_tool/CALM1_imputation.csv") %>% 
  mutate(gene = "CALM1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df3 = normalize_scoreset(df3, lower_bound = 0.1, upper_bound = 0.9)
df.out = de_noise_ss_1gene(df = df3, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df3 = merge(df3, df.out, all=T)

gene_id = "TPK1"
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df4 = read_csv("weile_tool/TPK1_imputation.csv") %>% 
  mutate(gene = "TPK1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df4 = normalize_scoreset(df4, lower_bound = 0.1, upper_bound = 0.9)
df.out = de_noise_ss_1gene(df = df4, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df4 = merge(df4, df.out, all=T)

df.out = rbind(df, df2, df3, df4)


# refined score was already normalized to [0-1]
# ensure refine score have the same scale and direction as FUSE score
df.out$imputed_score = -df.out$imputed_score+1
df.out$refined_score = -df.out$refined_score+1

# # attach SIFT and PROVEAN score
# ind = match(x = df.out$gene_aa_str, table = df.NSFP$gene_aa_str)
# df.out$SIFT_rankscore = df.NSFP$SIFT_converted_rankscore[ind]
# df.out$PROVEAN_rankscore = df.NSFP$PROVEAN_converted_rankscore[ind]

# load alphaMissense scores and attach to df.out
df.am = read_csv("./ProteinGym/combined_48genes_am_scores.csv")
ind = match(df.out$gene_aa_str, table = df.am$gene_aa_str)
df.out$am_score = df.am$am_pathogenicity[ind]

star_level = 1
df.clinvar2 = df.clinvar %>% 
  mutate(gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt)) %>%
  filter(gold_star >= star_level)

for (data_selected in c("BRCA1", "4genes")){
  df.out2 = df.out
  
  # attach clinvar
  ind = match(df.out2$gene_aa_str, table = df.clinvar2$gene_aa_str)
  df.out2$cs_simple = df.clinvar2$cs_simple[ind]
  df.out2$cs_numeric = NA
  df.out2$cs_numeric[df.out2$cs_simple == "p"] = 1
  df.out2$cs_numeric[df.out2$cs_simple == "b"] = 0
  
  # logit regression
  # leave-one-out: 1 for testing, rest for training
  if (data_selected == "BRCA1"){
    df.out2 = df.out2 %>% filter(gene == "BRCA1")
  }
  
  df.out2 = df.out2 %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
  df.out2 = df.out2 %>% filter(!is.na(cs_numeric) & !is.na(imputed_score) & !is.na(refined_score) & !is.na(sub_score_ss) & !is.na(pos_score) & !is.na(acc) & !is.na(am_score))
  # df.out2 = df.out2 %>% filter(!is.na(cs_numeric) & !is.na(imputed_score) & !is.na(refined_score) & !is.na(sub_score_ss) & !is.na(pos_score) & !is.na(acc)  & !is.na(am_score))
  
  df.test_all = c()
  for (i in 1:nrow(df.out2)){
    df.train = df.out2[-i,]
    df.test = df.out2[i,]
    
    model <- glm(cs_numeric~imputed_score, family="binomial", data=df.train)
    df.test$model1_predicted = predict(model, df.test, type="response")
    
    model <- glm(cs_numeric~refined_score, family="binomial", data=df.train)
    df.test$model2_predicted = predict(model, df.test, type="response")
    
    model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc, family="binomial", data=df.train)
    df.test$model3_predicted = predict(model, df.test, type="response")
    
    # model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc+SIFT_rankscore+PROVEAN_rankscore, family="binomial", data=df.train)
    # df.test$model4_predicted = predict(model, df.test, type="response")
    
    model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc+am_score, family="binomial", data=df.train)
    df.test$model4_predicted = predict(model, df.test, type="response")
    
    model <- glm(cs_numeric~norm_raw_score, family="binomial", data=df.train)
    df.test$model5_predicted = predict(model, df.test, type="response")
    
    model <- glm(cs_numeric~am_score, family="binomial", data=df.train)
    df.test$model6_predicted = predict(model, df.test, type="response")
    
    df.test_all = rbind(df.test_all, df.test)
  }
  
  # re-create Figure 2
  # Figure 2A-E
  # ymax = 2
  # bw = 0.2
  # x_range = c(-3, 3)
  # plt1 = draw_variant_ggplot_v2(df.test_all, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw)
  # plt2 = draw_variant_ggplot_v2(df.test_all, score_type = "final_score", score_type_name = "FUSE score", ymax = ymax, x_range = x_range, bw = bw)
  # plt3 = draw_variant_ggplot_v2(df.test_all, score_type = "imputed_score", score_type_name = "Weile, et al. imputed score", ymax = ymax, x_range = x_range, bw = bw)
  # plt4 = draw_variant_ggplot_v2(df.test_all, score_type = "refined_score", score_type_name = "Weile, et al. refined score", ymax = ymax, x_range = x_range, bw = bw)
  
  roc1 <- roc(df.test_all$cs_numeric, df.test_all$model1_predicted, auc=T, direction = "<")
  roc2 <- roc(df.test_all$cs_numeric, df.test_all$model2_predicted, auc=T, direction = "<")
  roc3 <- roc(df.test_all$cs_numeric, df.test_all$model3_predicted, auc=T, direction = "<")
  roc4 <- roc(df.test_all$cs_numeric, df.test_all$model4_predicted, auc=T, direction = "<")
  roc6 <- roc(df.test_all$cs_numeric, df.test_all$model6_predicted, auc=T, direction = "<")
  
  plot_margin = rep(10, 4)
  
  roc_labels = c(paste0("Weile imputed, AUC=", round(roc1$auc, digits = 2)),
                 paste0("Weile refined, AUC=", round(roc2$auc, digits = 2)))
  list.roc = list(roc1,roc2)
  names(list.roc) = roc_labels
  plt2 = ggroc(list.roc) +
    annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
    scale_color_discrete(name = NULL) + 
    ggtitle("Weile et al.") + xlab("Specificity") + ylab("Sensitivity") + 
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"),
          plot.margin = unit(plot_margin, "points"),
          legend.position = c(0.65, 0.2), 
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))
  
  # roc_labels = c(paste0("FUSE score, AUC=", round(roc3$auc, digits = 2)),
  #                paste0("FUSE+PROVEAN+SIFT,\nAUC=", round(roc4$auc, digits = 2)))
  roc_labels = c(paste0("FUSE score, AUC=", round(roc3$auc, digits = 2)),
                 paste0("AM score, AUC=", round(roc6$auc, digits = 2)),
                 paste0("FUSE + AM,AUC=", round(roc4$auc, digits = 2)))
  list.roc = list(roc3, roc6, roc4)
  names(list.roc) = roc_labels
  plt3 = ggroc(list.roc) +
    annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
    scale_color_discrete(name = NULL) + 
    ggtitle("FUSE pipeline") + xlab("Specificity") + ylab("Sensitivity") + 
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = unit(plot_margin, "points"),
          legend.position = c(0.65, 0.2), 
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))
  
  plt = plt2|plt3 + plot_annotation(tag_levels = 'A')
  ggsave(paste0("./manuscript/new_figures/FUSE_vs_Weile_", star_level, "starPlus_", data_selected, "_2.png"), plt, device = "png", width = 8, height = 4)
}


