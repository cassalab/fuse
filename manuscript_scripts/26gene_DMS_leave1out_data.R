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

# leave-one-out analysis for 26 gene's DMS
include_LOF = T
df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_", ifelse(include_LOF, "", "noLOF_"), "042423.rds"))
ls.norm_score = readRDS("maveDB/norm_score_26genes_042423.rds")
df.clinvar = read.table("clinvar_lite_GRCh38_081722.txt", header = T, sep = "\t")
data(BLOSUM62)
aa_list = colnames(df.funsum_all)
df.b62 = -BLOSUM62[aa_list, aa_list]
# scale df.b62 so it has the same mean and SD as df.funsum_all
df.b62 = (df.b62 - mean(df.b62))/sd(df.b62)*sd(df.funsum_all, na.rm = T) + mean(df.funsum_all, na.rm = T)
df.b62_sub_tb = funsum_to_subTable(df.b62)

for (star_level in c(0,1,2)){
  df.clinvar2 = df.clinvar %>% 
    mutate(gene_aa_str = paste0(gene, "---", aaref, aapos, aaalt)) %>%
    filter(gold_star >= star_level)
  
  df.test_all = c()
  all_gene_ids = names(ls.norm_score)
  for (gene_id in all_gene_ids){
    funsum_gene_ids = all_gene_ids[!(all_gene_ids %in% gene_id)]
    
    # calculate FUNSUM from all other scoresets
    df.all_sub_score = c()
    for (funsum_gene_id in funsum_gene_ids){
      df2 = ls.norm_score[[funsum_gene_id]]
      df.sub_score = get_norm_score(df2, pos_mean_method = "js", include_LOF = include_LOF, score_type = "norm_raw_score", sd_norm = F)
      df.all_sub_score = rbind(df.all_sub_score, cbind(funsum_gene_id, df.sub_score))
    }
    df.funsum_all = get_FUNSUM(df.all_sub_score[,-1], avg_method = "js", include_LOF = include_LOF ) # get combined funsum
    
    # de-noise DMS data
    df.out_all = c()
    for (gene_id2 in all_gene_ids){
      df = ls.norm_score[[gene_id2]]
      df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
      df.out_all = rbind(df.out_all, df.out)
    }
    
    # attach B262 scores
    ind = match(df.out_all$aa_pair, table = df.b62_sub_tb$aa_pair)
    df.out_all$b62_score = df.b62_sub_tb$score[ind]
    df.out_all$final_score_b62 = df.out_all$pos_score + df.out_all$b62_score
    df.out_all$final_score_b62_lite = df.out_all$final_score_b62
    ind = which(is.na(df.out_all$norm_raw_score))
    df.out_all$final_score_b62_lite = NA
    
    # attach clinVar
    ind = match(df.out_all$gene_aa_str, table = df.clinvar2$gene_aa_str)
    df.out_all$cs_simple = df.clinvar2$cs_simple[ind]
    df.out_all$cs_numeric = NA
    df.out_all$cs_numeric[df.out_all$cs_simple == "p"] = 1
    df.out_all$cs_numeric[df.out_all$cs_simple == "b"] = 0
    
    ind = which(!is.na(df.out_all$cs_numeric))
    df.out_all = df.out_all[ind,]
    
    # separate training and testing data
    ind = which(df.out_all$gene %in% gene_id)
    if (length(ind)>0){
      df.train = df.out_all[-ind,]
      df.test = df.out_all[ind,]
      
      data = df.train %>% filter(!is.na(norm_raw_score))
      model <- glm(cs_numeric~norm_raw_score, family="binomial", data=data)
      df.test$model1_predicted = predict(model, df.test, type="response")
      
      data = df.train %>% filter(!is.na(final_score))
      model <- glm(cs_numeric~final_score, family="binomial", data=data)
      df.test$model2_predicted = predict(model, df.test, type="response")
      
      data = df.train %>% filter(!is.na(sub_score) & !is.na(pos_score))
      model <- glm(cs_numeric~sub_score+pos_score, family="binomial", data=data)
      df.test$model3_predicted = predict(model, df.test, type="response")
      
      data = df.train %>% filter(!is.na(sub_score))
      model <- glm(cs_numeric~sub_score, family="binomial", data=data)
      df.test$model4_predicted = predict(model, df.test, type="response")
      
      data = df.train %>% filter(!is.na(b62_score))
      model <- glm(cs_numeric~b62_score, family="binomial", data=data)
      df.test$model5_predicted = predict(model, df.test, type="response")
      
      data = df.train %>% filter(!is.na(b62_score) & !is.na(sub_score))
      model <- glm(cs_numeric~b62_score+sub_score, family="binomial", data=data)
      df.test$model6_predicted = predict(model, df.test, type="response")
      
      data = df.train %>% filter(!is.na(final_score_b62) & !is.na(sub_score))
      model <- glm(cs_numeric~final_score_b62, family="binomial", data=data)
      df.test$model7_predicted = predict(model, df.test, type="response")
      
      df.test_all = rbind(df.test_all, df.test)
    }
  }
  write_csv(x = df.test_all, file = paste0("./maveDB/combined_26genes_leaveoneout_", star_level, "starPlus.csv"))
}
