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
include_LOF = T
df.funsum_all = readRDS("./shiny_app/DMS_denoiser/funsum_maveDB_noLOF_042423.rds")

# add DMS datasets
df.out_all = read_csv("./maveDB/FUSE_estimated_26genes.csv")
df.obs = data.frame(facet = "DMS", gene = unique(df.out_all$gene), Predicted=NA, Assayed=NA)
df.obs = df.obs %>% filter(gene != "LDLR" & gene != "AID") # remove LDLR and AID
df.obs$facet[df.obs$gene == "BRCA1"] = "SGE"

for (i in 1:nrow(df.obs)){
  df = df.out_all %>% filter(gene == df.obs$gene[i])
  df1 = df %>% filter(!is.na(norm_raw_score))
  df.obs$Predicted[i] = nrow(df) - nrow(df1)
  df.obs$Assayed[i] = nrow(df1)
}

# add base editng datasets
gene_id = "BRCA1"
df.raw1 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
df.raw2 = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
df = rbind(df.raw1, df.raw2)
df = normalize_scoreset(df, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df$gene = gene_id
df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
temp = data.frame(facet = "Base editing", gene = gene_id, 
                  Predicted = length(which(is.na(df.out$norm_raw_score))), 
                  Assayed = length(which(!is.na(df.out$norm_raw_score))))
df.obs = rbind(df.obs, temp)

df.raw = read.csv("./shiny_app/DMS_denoiser/DMS_data/DDR_genes_BE.csv", header = T)
df = normalize_scoreset(df.raw, lower_bound = 0.1, upper_bound = 0.9, force_quantile = T)
df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
temp = data.frame(facet = "Base editing", gene = "DDR genes", 
                  Predicted = length(which(is.na(df.out$norm_raw_score))), 
                  Assayed = length(which(!is.na(df.out$norm_raw_score))))
df.obs = rbind(df.obs, temp)
df.obs$facet = factor(df.obs$facet, levels = c("DMS", "Base editing", "SGE"))

df.obs_long = df.obs %>% pivot_longer(cols = c("Predicted", "Assayed"), names_to = "Score_type", values_to = "Variant_count")
df.obs_long$Score_type = factor(df.obs_long$Score_type, levels = c("Predicted", "Assayed"))

plt <- ggplot(df.obs_long, aes(fill=Score_type, y=gene, x=Variant_count)) + 
  geom_col(position="stack", show.legend = T, color="black") +
  scale_fill_manual(values = c("#FD8972", "#8EDAF5"), name = "") + 
  ggtitle("Variants with functional estimation\n") + xlab(label = "") + ylab(label = "Count")  + 
  facet_grid(rows = vars(facet), scales='free', space = "free", switch = "x") +
  theme_ipsum(base_family = "Arial", plot_title_size = 15, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust=0.5, vjust=2), 
        legend.position = c(0.9, 0.9), 
        legend.background = element_rect(fill="white", linetype=0),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plt
ggsave(paste0("~/Google Drive/My Drive/manuscript_resub/variant_count.png"), plt, device = "png", width = 6, height = 10) 
