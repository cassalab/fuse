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

include_LOF = T
df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_", ifelse(include_LOF, "", "noLOF_"), "042423.rds"))
df.clinvar = read.table("clinvar_lite_GRCh38_081722.txt", header = T, sep = "\t")

df = read_csv("weile_tool/_P38398[BRCA1_BRCT]_imputation.csv") %>% 
  mutate(gene = "BRCA1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9)

df2 = read_csv("weile_tool/_P38398[BRCA1_RING]_imputation.csv") %>% 
  mutate(gene = "BRCA1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df2 = normalize_scoreset(df2, lower_bound = 0.1, upper_bound = 0.9)

df3 = read_csv("weile_tool/CALM1_imputation.csv") %>% 
  mutate(gene = "CALM1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df3 = normalize_scoreset(df3, lower_bound = 0.1, upper_bound = 0.9)

df4 = read_csv("weile_tool/TPK1_imputation.csv") %>% 
  mutate(gene = "TPK1", aapos = aa_pos, aaref = aa_ref, aaalt = aa_alt, raw_score = fitness_input, imputed_score = fitness_imputed, refined_score = fitness_refine) %>%
  select(gene, aapos, aaref, aaalt, raw_score, imputed_score, refined_score)
df4 = normalize_scoreset(df4, lower_bound = 0.1, upper_bound = 0.9)

df.refined_all = rbind(df, df2, df3, df4) %>% mutate(gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt))
df.out = de_noise(df = df.refined_all, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF, show_func_class = T)

# attach inputed and refined score
ind = match(x = df.out$gene_aa_str, table = df.refined_all$gene_aa_str)
# refined score was already normalized to [0-1]
# ensure refine score have the same scale and direction as FUSE score
df.out$imputed_score = -df.refined_all$imputed_score[ind]+1
df.out$refined_score = -df.refined_all$refined_score[ind]+1

# attach SIFT and PROVEAN score
ind = match(x = df.out$gene_aa_str, table = df.NSFP$gene_aa_str)
df.out$SIFT_rankscore = df.NSFP$SIFT_converted_rankscore[ind]
df.out$PROVEAN_rankscore = df.NSFP$PROVEAN_converted_rankscore[ind]

for (star_level in c(0,1,2)){
  df.clinvar2 = df.clinvar %>% 
    mutate(gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt)) %>%
    filter(gold_star >= star_level)
  
  for (data_selected in c("BRCA1", "4genes")){
    for (vars_included in c("misOnly", "allVars")){
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
      if (vars_included == "misOnly"){
        df.out2 = df.out2 %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
      }
      df.out2 = df.out2 %>% filter(!is.na(cs_numeric))
      
      df.test_all = c()
      for (i in 1:nrow(df.out2)){
        df.train = df.out2[-i,]
        df.test = df.out2[i,]
        
        data = df.train %>% filter(!is.na(norm_raw_score))
        model <- glm(cs_numeric~norm_raw_score, family="binomial", data=data)
        df.test$model1_predicted = predict(model, df.test, type="response")
        
        data = df.train %>% filter(!is.na(final_score))
        model <- glm(cs_numeric~final_score, family="binomial", data=data)
        df.test$model2_predicted = predict(model, df.test, type="response")
        
        data = df.train %>% filter(!is.na(sub_score))
        model <- glm(cs_numeric~sub_score, family="binomial", data=data)
        df.test$model3_predicted = predict(model, df.test, type="response")
        
        data = df.train %>% filter(!is.na(imputed_score))
        model <- glm(cs_numeric~imputed_score, family="binomial", data=data)
        df.test$model4_predicted = predict(model, df.test, type="response")
        
        data = df.train %>% filter(!is.na(refined_score))
        model <- glm(cs_numeric~refined_score, family="binomial", data=data)
        df.test$model5_predicted = predict(model, df.test, type="response")
        
        data = df.train %>% filter(!is.na(final_score) & !is.na(SIFT_rankscore) & !is.na(PROVEAN_rankscore))
        model <- glm(cs_numeric~final_score+SIFT_rankscore+PROVEAN_rankscore, family="binomial", data=data)
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
      roc5 <- roc(df.test_all$cs_numeric, df.test_all$model5_predicted, auc=T, direction = "<")
      roc6 <- roc(df.test_all$cs_numeric, df.test_all$model6_predicted, auc=T, direction = "<")
      
      plot_margin = rep(10, 4)
      roc_labels = c(paste0("Original score, AUC=", round(roc1$auc, digits = 2)))
      list.roc = list(roc1)
      names(list.roc) = roc_labels
      plt1 = ggroc(list.roc) +
        annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
        scale_color_discrete(name = NULL) + 
        ggtitle("Original assay") + xlab("Specificity") + ylab("Sensitivity") + 
        theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.margin = unit(plot_margin, "points"),
              legend.position = c(0.65, 0.15), 
              legend.text = element_text(size = 10),
              legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))
      
      roc_labels = c(paste0("Weile imputed, AUC=", round(roc4$auc, digits = 2)),
                     paste0("Weile refined, AUC=", round(roc5$auc, digits = 2)))
      list.roc = list(roc4,roc5)
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
      
      roc_labels = c(paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)),
                     paste0("FUSE+PROVEAN+SIFT,\nAUC=", round(roc6$auc, digits = 2)))
      list.roc = list(roc2,roc6)
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
      
      plt = plt1|plt2|plt3 + plot_annotation(tag_levels = 'A')
      ggsave(paste0("~/Google Drive/My Drive/manuscript_resub/ClinVar_plots/", star_level, "starPlus/", "FUSE_vs_Weile_", star_level, "starPlus_", data_selected, "_", vars_included, ".png"), plt, device = "png", width = 11, height = 4)
      
      # list.plt = list()
      # for (i in 1:length(list.roc)){
      #   plt = ggroc(list.roc[i]) +
      #     annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
      #     scale_color_discrete(name = NULL, label = sub(pattern = ", ", replacement = "\n", x = names(list.roc)[i])) + 
      #     ggtitle("Variant classification") + xlab("Specificity") + ylab("Sensitivity") + 
      #     theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
      #     theme(plot.title = element_text(hjust = 0.5),
      #           legend.position = c(0.7, 0.2), 
      #           legend.text = element_text(size = 10),
      #           legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))
      #   list.plt[[i]] = plt
      # }
      
      # plt = (plt1|list.plt[[1]])/(list.plt[[2]]|list.plt[[3]])/(list.plt[[4]]|list.plt[[5]])
      # plt = (list.plt[[1]]|list.plt[[2]]|list.plt[[3]])/(list.plt[[4]]|list.plt[[5]]|plot_spacer()) + plot_annotation(tag_levels = 'A')
      # ggsave(paste0("~/Google Drive/My Drive/manuscript_resub/ClinVar_plots/", star_level, "starPlus/", "FUSE_vs_Weile_", star_level, "starPlus_", data_selected, "_", vars_included, ".png"), plt, device = "png", width = 12, height = 8) 
    }
  }
}

