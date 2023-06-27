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

# re-create figure 2 format
for (star_level in c(0,1,2)){
  df.test_all = read_csv(paste0("./maveDB/combined_26genes_leaveoneout_", star_level,"starPlus.csv"))
  df.test_all = df.test_all %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
  
  selected_gene_ids = c("TP53", "PTEN", "HRAS", "MSH2", "GCK", "MTHFR", "CBS", "HMBS", "SNCA") # scoresets that have at least 3 variants in b or p group
  df.test_all = df.test_all[(df.test_all$gene %in% selected_gene_ids),]
  
  ymax = 2
  bw = 0.2
  x_range = c(-3, 3)
  
  # Figure 2A-E
  plt1 = draw_variant_ggplot_v2(df.test_all, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
  plt2 = draw_variant_ggplot_v2(df.test_all, score_type = "pos_score", score_type_name = "Positional component score", ymax = ymax, x_range = x_range, bw = bw)
  plt3 = draw_variant_ggplot_v2(df.test_all, score_type = "sub_score", score_type_name = "Substitution component score", ymax = ymax, x_range = x_range, bw = bw)
  plt4 = draw_variant_ggplot_v2(df.test_all, score_type = "final_score_lite", score_type_name = "FUSE score\n(without predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
  plt5 = draw_variant_ggplot_v2(df.test_all, score_type = "final_score", score_type_name = "FUSE score\n(with predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
  
  df.test_all = df.test_all %>% filter(!is.na(norm_raw_score))
  roc1 <- roc(df.test_all$cs_numeric, df.test_all$model1_predicted, auc=T, direction = "<")
  roc2 <- roc(df.test_all$cs_numeric, df.test_all$model2_predicted, auc=T, direction = "<")
  roc3 <- roc(df.test_all$cs_numeric, df.test_all$model3_predicted, auc=T, direction = "<")
  
  roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)), 
                 paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)))
  plt6 = ggroc(list(roc1,roc2)) +
    annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
    scale_color_discrete(name = NULL, labels = roc_labels) + 
    ggtitle("Variant classification") + xlab("Specificity") + ylab("Sensitivity") + 
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.7, 0.2), 
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))
  
  plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
  ggsave(paste0("~/Google Drive/My Drive/manuscript_resub/ClinVar_plots/",star_level, "starPlus/9genes_DMS_leaveoneout_", star_level, "starPlus.png"), plt, device = "png", width = 12, height = 8) 
}

