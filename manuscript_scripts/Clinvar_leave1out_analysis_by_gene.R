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

# leave-one-out analysis 
ss_list = c("G" = "Helices", "H" = "Helices", "I" = "Helices", 
            "E" = "Strands", "B" = "Strands", 
            "T" = "Loops", "S" = "Loops", "C" = "Loops")
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_052024.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")



# denoise functional data
df.out_all = c()
all_scoreset_ids = unique(df.combined_scoreset_info$scoreset_id)
for (scoreset_id in all_scoreset_ids){
  df = ls.combined_scoreset[[scoreset_id]]
  ind = which(df.combined_scoreset_info$scoreset_id == scoreset_id)
  df$gene = df.combined_scoreset_info$gene_symbol[ind]
  
  # de-noise DMS data
  df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                             dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                             include_LOF = include_LOF, show_func_class = T)
  df.out$scoreset_id = scoreset_id
  df.out_all = rbind(df.out_all, df.out)
}

# attach alphaMissense scores
df.am = read_csv("./ProteinGym/combined_48genes_am_scores.csv")
ind = match(df.out_all$gene_aa_str, table = df.am$gene_aa_str)
df.out_all$am_score = df.am$am_pathogenicity[ind]

# attach clinvar annotations
df.clinvar = read.table("./clinvar/clinvar_lite_GRCh38_091323.txt", header = T, sep = "\t")
df.clinvar = df.clinvar %>% mutate(gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt))
ind = match(df.out_all$gene_aa_str, table = df.clinvar$gene_aa_str)
df.out_all$cs_simple = df.clinvar$cs_simple[ind]
df.out_all$cs_numeric = NA
df.out_all$cs_numeric[df.out_all$cs_simple == "p"] = 1
df.out_all$cs_numeric[df.out_all$cs_simple == "b"] = 0
df.out_all$clinvar_star_level = df.clinvar$gold_star[ind]


### clinvar composite plot (Figure 2)
for (star_level in c(0,1,2)){
  df.out_all2 = df.out_all %>% 
    filter(!is.na(cs_numeric) & !is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(am_score) & clinvar_star_level >= star_level) %>%
    arrange(gene, scoreset_id, aapos, aaref, aaalt)
  
  df.test_all = c()
  df.out_all3 = c()
  
  # df.combined_scoreset_info2 = df.combined_scoreset_info %>% arrange(gene_symbol, desc(norm_raw_score_AUC), desc(pos_mean_coverage), desc(MIS_ct)) %>%
  #   group_by(gene_symbol) %>% summarise(scoreset_id = scoreset_id[1])
  df.combined_scoreset_info2 = df.combined_scoreset_info %>% arrange(gene_symbol, desc(pos_mean_coverage), desc(MIS_ct)) %>%
    group_by(gene_symbol) %>% summarise(scoreset_id = scoreset_id[1])
  
  for (scoreset_id in unique(df.out_all2$scoreset_id)){
    gene = df.out_all2$gene[df.out_all2$scoreset_id == df.out_all2][1]
    
    # for density plot
    ind = which(df.out_all$scoreset_id == scoreset_id)
    df.out_all3 = rbind(df.out_all3, df.out_all[ind,])
    
    # for ROC plot
    ind = which(df.out_all2$scoreset_id == scoreset_id)
    ind2 = which(df.out_all2$gene == gene)
    
    if (length(ind)>0){
      # separate training and testing data
      df.train = df.out_all2[-ind2,]
      df.test = df.out_all2[ind,]
      
      model <- glm(cs_numeric~norm_raw_score, family="binomial", data=df.train)
      df.test$model1_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~final_score, family="binomial", data=df.train)
      df.test$model2_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~final_score_ss, family="binomial", data=df.train)
      df.test$model3_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc, family="binomial", data=df.train)
      df.test$model4_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc+am_score, family="binomial", data=df.train)
      df.test$model5_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~am_score, family="binomial", data=df.train)
      df.test$model6_predicted = predict(model, df.test, type="response")
      
      df.test_all = rbind(df.test_all, df.test)
    }
  }
  
  ### draw Figure 2
  dir.create("./manuscript/new_figures/clinvar_composite/", showWarnings = F, recursive = T)
  
  df.out_all3 = df.out_all3 %>% filter(!is.na(cs_numeric) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(am_score) & clinvar_star_level >= star_level)
  df.out_all3$final_score_ss_lite = df.out_all3$final_score_ss
  df.out_all3$final_score_ss_lite[is.na(df.out_all3$raw_score)] = NA
  
  ymax = 2
  bw = 0.1
  x_range = c(-2, 3)
  
  # Figure 2A-E
  plt1 = draw_variant_ggplot_v2(df.out_all3, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
  plt2 = draw_variant_ggplot_v2(df.out_all3, score_type = "pos_score", score_type_name = "Positional component score", ymax = ymax, x_range = x_range, bw = bw)
  plt3 = draw_variant_ggplot_v2(df.out_all3, score_type = "sub_score_ss", score_type_name = "Substitution component score", ymax = ymax*2, x_range = x_range, bw = bw)
  plt4 = draw_variant_ggplot_v2(df.out_all3, score_type = "final_score_ss_lite", score_type_name = "FUSE score\n(without predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
  plt5 = draw_variant_ggplot_v2(df.out_all3, score_type = "final_score_ss", score_type_name = "FUSE score\n(with predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
  
  df.test_all = df.test_all %>% filter(!is.na(norm_raw_score))
  roc1 <- roc(df.test_all$cs_numeric, df.test_all$model1_predicted, auc=T, direction = "<")
  roc2 <- roc(df.test_all$cs_numeric, df.test_all$model2_predicted, auc=T, direction = "<")
  roc3 <- roc(df.test_all$cs_numeric, df.test_all$model3_predicted, auc=T, direction = "<")
  roc4 <- roc(df.test_all$cs_numeric, df.test_all$model4_predicted, auc=T, direction = "<")
  roc5 <- roc(df.test_all$cs_numeric, df.test_all$model5_predicted, auc=T, direction = "<")
  roc6 <- roc(df.test_all$cs_numeric, df.test_all$model6_predicted, auc=T, direction = "<")
  roc_all = list(roc1,roc4,roc5,roc6)
  
  roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
                 paste0("FUSE, AUC=", round(roc4$auc, digits = 2)),
                 paste0("FUSE + AM, AUC=", round(roc5$auc, digits = 2)),
                 paste0("AM, AUC=", round(roc6$auc, digits = 2)))
  
  plt6 = ggroc(roc_all) + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
    scale_color_discrete(name = NULL, labels = roc_labels) + 
    ggtitle("Variant classification") + xlab("Specificity") + ylab("Sensitivity") + 
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.72, 0.28), 
          legend.text = element_text(size = 8, margin = margin(t = 0, r = 1, b = 0, l = 0)),
          legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2),
          legend.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "pt"))
  
  plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
  ggsave(paste0("./manuscript/new_figures/clinvar_composite/DMS_leaveoneout_", star_level, "starPlus_3.png"), plt, device = "png", width = 12, height = 8) 
  
  
  # generate data table for figure 2
  df.final = df.test_all %>% 
    mutate(ClinVar_ClnSig = ifelse(cs_simple == "b", "benign/likely benign", "pathogenic/likely pathogenic")) %>%
    select(gene, scoreset_id, aapos, aaref, aaalt, functional_class, gene_aa_str, DSSP_ss = ss, DSSP_acc = acc, 
           raw_score, norm_raw_score, pos_score, sub_score_ss, FUSE_score = final_score_ss, 
           Alpha_Missense_score = am_score, ClinVar_ClnSig, ClinVar_gold_star = clinvar_star_level,
           raw_score_pathogenic_prob = model1_predicted, 
           FUSE_pathogenic_prob = model4_predicted, 
           FUSE_AM_pathogenic_prob = model5_predicted,
           AM_pathogenic_prob = model6_predicted) %>%
    arrange(gene, scoreset_id, aapos, aaref, aaalt)
  write.csv(df.final, file = paste0("./manuscript/new_figures/supplementary_table/Supp_table3__DMS_leave1out_", star_level, "starPlus.csv"))
}








### draw supplement figure S4 (ROC plot for individual DMS scoresets)
# filter for scoresets that have at least 3 pathogenic and benign variants

df.roc_info = c()
for (star_level in c(0,1,2)){
  df.out_all2 = df.out_all %>% filter(!is.na(cs_numeric) & !is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(am_score) & clinvar_star_level >= star_level)
  
  df.combined_scoreset_info2 = df.combined_scoreset_info %>% filter(P_LP >= 1 & B_LB >= 1)
  # df.combined_scoreset_info2 = df.combined_scoreset_info2 %>% arrange(gene_symbol, desc(norm_raw_score_AUC), desc(pos_mean_coverage), desc(MIS_ct)) %>%
  #   group_by(gene_symbol) %>% summarise(scoreset_id = scoreset_id[1])
  df.combined_scoreset_info2 = df.combined_scoreset_info2 %>% arrange(gene_symbol, desc(pos_mean_coverage), desc(MIS_ct)) %>%
    group_by(gene_symbol) %>% summarise(scoreset_id = scoreset_id[1])
  
  df.test_all = c()
  for (scoreset_id in unique(df.combined_scoreset_info2$scoreset_id)){
    gene = df.out_all2$gene[df.out_all2$scoreset_id == scoreset_id][1]
    
    # for ROC plot
    ind = which(df.out_all2$scoreset_id == scoreset_id)
    ind2 = which(df.out_all2$gene == gene)
    
    if (length(ind)>0){
      # separate training and testing data
      df.train = df.out_all2[-ind2,]
      df.test = df.out_all2[ind,]
      
      model <- glm(cs_numeric~norm_raw_score, family="binomial", data=df.train)
      df.test$model1_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~final_score, family="binomial", data=df.train)
      df.test$model2_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~final_score_ss, family="binomial", data=df.train)
      df.test$model3_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc, family="binomial", data=df.train)
      df.test$model4_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~pos_score+sub_score+sub_score_ss+acc+am_score, family="binomial", data=df.train)
      df.test$model5_predicted = predict(model, df.test, type="response")
      
      model <- glm(cs_numeric~am_score, family="binomial", data=df.train)
      df.test$model6_predicted = predict(model, df.test, type="response")
      
      df.test_all = rbind(df.test_all, df.test)
    }
  }
  
  
  ls.plt = list()
  for (scoreset_id2 in unique(df.out_all2$scoreset_id)){
    df.test = df.test_all %>% filter(scoreset_id == scoreset_id2)
    gene = df.test$gene[df.test$scoreset_id == scoreset_id2][1]
    
    p_ct = length(which(df.test$cs_simple == "p"))
    b_ct = length(which(df.test$cs_simple == "b"))
    
    if (!(scoreset_id2 %in% df.combined_scoreset_info$scoreset_id) | scoreset_id2 %in% df.combined_scoreset_info2$scoreset_id){
      if (p_ct >= 1 & b_ct >= 1){
        roc1 <- roc(df.test$cs_numeric, df.test$norm_raw_score, auc=T, direction = "<")
        # roc1 <- roc(df.test$cs_numeric, df.test$model1_predicted, auc=T, direction = "<")
        roc4 <- roc(df.test$cs_numeric, df.test$model4_predicted, auc=T, direction = "<")
        roc5 <- roc(df.test$cs_numeric, df.test$model5_predicted, auc=T, direction = "<")
        roc6 <- roc(df.test$cs_numeric, df.test$model6_predicted, auc=T, direction = "<")
        roc_all = list(roc1,roc4,roc5, roc6)
        
        temp = data.frame(star_level, gene, p_ct, b_ct,
                          raw_score_auc = round(roc1$auc, digits = 2),
                          FUSE_auc = round(roc4$auc, digits = 2),
                          FUSE_AM_auc = round(roc5$auc, digits = 2),
                          AM_auc = round(roc6$auc, digits = 2))
        df.roc_info = rbind(df.roc_info, temp)
        
        roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
                       paste0("FUSE, AUC=", round(roc4$auc, digits = 2)),
                       paste0("FUSE + AM, AUC=", round(roc5$auc, digits = 2)),
                       paste0("AM, AUC=", round(roc6$auc, digits = 2)))
        
        plt = ggroc(roc_all) +
          annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
          scale_color_discrete(name = NULL, labels = roc_labels) + 
          ggtitle(gene) + xlab("Specificity") + ylab("Sensitivity") + 
          theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
          theme(plot.title = element_text(hjust = 0.5),
                legend.position = c(0.7, 0.25), 
                legend.text = element_text(size = 8, margin = margin(t = 0, r = 1, b = 0, l = 0)),
                legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2),
                legend.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "pt"))
        ls.plt[[scoreset_id2]] = plt
      }
    }
  }
  
  row_ct = ceiling(length(ls.plt) / 4)
  plt = wrap_plots(ls.plt, ncol = 4, nrow = row_ct)
  ggsave(paste0("./manuscript/new_figures/clinvar_composite/individual_roc_", star_level, "starPlus_4.png"), plt, device = "png", width = 15, height = row_ct*4, dpi=640)
}



# comapre AUC
df.roc_info_long = df.roc_info %>% 
  pivot_longer(cols = contains("auc"), names_to = "id", values_to = "AUC") %>%
  mutate(id = gsub(pattern="_auc", replacement="", id))
df.roc_info_long$id[df.roc_info_long$id == "raw_score"] = "Raw score"
df.roc_info_long$id[df.roc_info_long$id == "FUSE_AM"] = "FUSE + AM"
df.roc_info_long$star_level[df.roc_info_long$star_level == 0] = "0 star or more"
df.roc_info_long$star_level[df.roc_info_long$star_level == 1] = "1 star or more"
df.roc_info_long$star_level[df.roc_info_long$star_level == 2] = "2 star or more"
df.roc_info_long$id = factor(df.roc_info_long$id, levels = c("Raw score", "FUSE", "FUSE + AM", "AM"))

plt = df.roc_info_long %>% ggplot(aes(x=id, y=AUC, color=id, fill = id)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0, alpha=0.2) + 
  geom_jitter(width=0.25, size=1, alpha=0.6) + ylab("AUROC") + 
  facet_wrap(~star_level) + ylim(0.5,1) + xlab("") + theme_light() +
  theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plt

ggsave(paste0("./manuscript/new_figures/clinvar_composite/AUC_comparison.png"), plt, device = "png", width = 8, height = 4)




