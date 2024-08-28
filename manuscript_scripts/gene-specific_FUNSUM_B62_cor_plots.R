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

# construct correlation matrix
# make gene-specific FUNSUM correlation plot +  BLOSUM
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")

ls.funsum = list()
for (gene in unique(df.combined_scoreset_info$gene_symbol)){
  if (gene %in% names(ls.combined_scoreset)){
    df = ls.combined_scoreset[[gene]]
    df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF, score_type = "norm_raw_score", sd_norm = F)
    df.funsum = get_FUNSUM(df.sub_score, avg_method = "js", include_LOF = include_LOF)
  } else {
    ind = which(df.combined_scoreset_info$gene_symbol == gene)
    df.sub_score_all = c()
    for (i in ind){
      df = ls.combined_scoreset[[df.combined_scoreset_info$scoreset_id[i]]]
      df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF, score_type = "norm_raw_score", sd_norm = F)
      df.sub_score_all = rbind(df.sub_score_all, df.sub_score)
    }
    df.funsum = get_FUNSUM(df.sub_score_all, avg_method = "js", include_LOF = include_LOF)
  }
  
  ls.funsum[[gene]] = df.funsum
}

ls.funsum[["Overall FUNSUM"]] = df.funsum_all

data(BLOSUM62)
aa_list = colnames(df.funsum_all)
df.b62 = -BLOSUM62[aa_list, aa_list]
ls.funsum[["BLOSUM62"]] = df.b62

df.cor = matrix(NA, nrow = length(ls.funsum), ncol = length(ls.funsum))
for (i in 1:length(ls.funsum)){
  for (j in 1:length(ls.funsum)){
    temp = cor.test(as.numeric(ls.funsum[[i]]), as.numeric(ls.funsum[[j]]), method = "spearman")
    df.cor[i,j] = temp$estimate
  }
}
rownames(df.cor) = names(ls.funsum)
colnames(df.cor) = names(ls.funsum)

# correlation in longer format
df.cor_long = data.frame(df.cor) %>% rownames_to_column("data1") %>% 
  pivot_longer(cols = !data1, names_to = "data2", values_to = "pearson_cor") 

# filter for gene-specific FUNSUm and BLOSUM correaltion < 0.5
temp = df.cor_long %>% filter(data2 == "BLOSUM62" & pearson_cor < 0.5)

# correlation between overall FUNSUma nd BLOSUM62
df.cor['Overall FUNSUM', 'BLOSUM62']


# draw correlation plot
png("./manuscript/new_figures/FUNSUM_correlation_2.png", width = 8000, height = 8000, res = 800)
corrplot(df.cor, method = "color", order = "hclust", outline=T, col = colorRampPalette(c("blue","white","red"))(200),
         tl.col="black")
dev.off()

# png("~/Google Drive/My Drive/manuscript_resub/sup_figures/FUNSUM_correlation_2.png", width = 8000, height = 8000, res = 1000)
# corrplot(df.cor, method = "circle", order = "hclust", col = colorRampPalette(c("blue","white","red"))(200),
#          tl.col="black")
# dev.off()