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


# UKB patient anlysis on TP53 DMS data 
# prepare df.p2v
gene_id = "TP53"
df.NSFP_all = read.xlsx("dbNSFP_BRCA1_LDLR_TP53_includeLOF.xlsx")
df.NSFP = df.NSFP_all %>% 
  mutate(variant_str = paste(chr, pos, ref, alt, sep = "-"), 
         gene_aa_str = paste0(genename, "---", aaref, aapos, aaalt)) %>%
  filter(aaref != "X" & aaalt != "X" & genename == gene_id) %>%
  select(genename, variant_str, gene_aa_str)
df.p2v = read_ssv(paste0("./UKBB/maveDB_genes_v2p/", gene_id, ".ssv"))
ind = match(df.p2v$variant, table = df.NSFP$variant_str)
ind2 = which(!is.na(ind))
df.p2v = df.p2v[ind2,]

df.p2v$gene_aa_str = df.NSFP$gene_aa_str[ind[ind2]]
df.p2v$gene = df.NSFP$genename[ind[ind2]]

# attach functional scores
include_LOF = T
df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_", ifelse(include_LOF, "", "noLOF_"), "042423.rds"))
ls.norm_score = readRDS("maveDB/norm_score_26genes_042423.rds")
df = ls.norm_score$TP53
df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
ind = match(df.p2v$gene_aa_str, table = df.out$gene_aa_str)
temp = df.out[ind, c("norm_raw_score", "pos_score", "sub_score", "final_score", "final_score_lite")]
df.p2v = cbind(df.p2v, temp)
df.p2v = df.p2v %>%
  filter(!is.na(norm_raw_score) | !is.na(final_score)) %>% 
  group_by(patient) %>% 
  summarise(n_variants = length(variant), 
            variant_str = paste0(variant, collapse = ','), 
            AF = mean(AF),
            gene_aa_str = paste0(gene_aa_str, collapse = ','),
            gene = paste0(gene, collapse = ','),
            norm_raw_score = mean(norm_raw_score, na.rm=T),
            pos_score = mean(pos_score, na.rm=T),
            sub_score = mean(sub_score, na.rm=T),
            final_score = mean(final_score, na.rm=T),
            final_score_lite = mean(final_score_lite, na.rm=T)) %>%
  filter(n_variants == 1)

# attach cancer phenotypes
LFS_cancer_types = c("sarcomafibrosarcoma",
                     "acute_myeloid_leukaemia", 
                     "adrenal_cancer", 
                     "bone_metastases__bony_secondaries", 
                     "brain_cancer__primary_malignant_brain_tumour", 
                     "breast_cancer", 
                     "primary_bone_cancer",
                     "leukaemia",
                     "hodgkins_lymphoma__hodgkins_disease", "non-hodgkins_lymphoma",
                     "lung_cancer",
                     "skin_cancer",
                     "kidneyrenal_cell_cancer",
                     "colon_cancersigmoid_cancer", "pancreas_cancer",
                     "thyroid_cancer",
                     "ovarian_cancer", "testicular_cancer", "prostate_cancer"
)
df.cancer = readRDS("UKBB/all_cancer_combined.rds")
df.cancer$cancer = 0
for (cancer_type in LFS_cancer_types){
  ind = which(df.cancer[[cancer_type]] == 1)
  df.cancer$cancer[ind] = 1
}
df.cancer = df.cancer %>% 
  select(eid, sex, age, cancer) %>%
  mutate(eid = as.character(eid))
ind = match(df.p2v$patient, table = df.cancer$eid)
df.p2v = cbind(df.p2v, df.cancer[ind,])
df.p2v$cancer = factor(df.p2v$cancer)

# score distribution between patient groups
ymax = 2
bw = 0.2
x_range = c(-3.5, 4.5)
plt1 = draw_patient_ggplot(df.p2v, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
plt2 = draw_patient_ggplot(df.p2v, score_type = "pos_score", score_type_name = "Positional component score", ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot(df.p2v, score_type = "sub_score", score_type_name = "Substitution component score", ymax = ymax, x_range = x_range, bw = bw)
plt4 = draw_patient_ggplot(df.p2v, score_type = "final_score_lite", score_type_name = "FUSE score\n(without predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
plt5 = draw_patient_ggplot(df.p2v, score_type = "final_score", score_type_name = "FUSE score\n(with predicted sites)", ymax = ymax, x_range = x_range, bw = bw)

# ROC analysis with logistic regression model
age_thres1 = 100
age_thres2 = 0

data = df.p2v %>% 
  mutate(func_score = norm_raw_score) %>% 
  filter(is.finite(func_score) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))
cancer_ct = length(which(data$cancer == 1))
if (cancer_ct>0){
  model1 <- glm(cancer~func_score, family="binomial", data=data)
  data$predicted = predict(model1, data, type="response")
  roc1 <- roc(data$cancer, data$predicted, auc=T, direction = "<")
}

data = df.p2v %>% 
  mutate(func_score = final_score) %>% 
  filter(is.finite(func_score) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))
cancer_ct = length(which(data$cancer == 1))
if (cancer_ct>0){
  model2 <- glm(cancer~func_score, family="binomial", data=data)
  data$predicted = predict(model2, data, type="response")
  roc2 <- roc(data$cancer, data$predicted, auc=T, direction = "<")
}

# plot ROC curves
roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)), 
               paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)))

plt6 = ggroc(list('Raw score' = roc1, 'Estimated score' = roc2)) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("LFS-associated cancer classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.2), 
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))

plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
ggsave("~/Google Drive/My Drive/manuscript_resub/UKB_plots/LFScancer_allPlots_TP53_DMS.png", plt, device = "png", width = 12, height = 8)
# ggsave(filename = "./maveDB/plots/LFScancer_ROC_on_TP53_DMS.png", plot = plt, device = "png", width = 7, height = 5)



# UKB patient anlysis on BRCA1 SGE data
# prepare df.p2v
gene_id = "BRCA1"
df.NSFP_all = read.xlsx("dbNSFP_BRCA1_LDLR_TP53_includeLOF.xlsx")
df.NSFP = df.NSFP_all %>% 
  mutate(variant_str = paste(chr, pos, ref, alt, sep = "-"), 
         gene_aa_str = paste0(genename, "---", aaref, aapos, aaalt)) %>%
  filter(aaref != "X" & aaalt != "X" & genename == gene_id) %>%
  select(genename, variant_str, gene_aa_str)
df.p2v = read_ssv(paste0("./UKBB/maveDB_genes_v2p/", gene_id, ".ssv"))
ind = match(df.p2v$variant, table = df.NSFP$variant_str)
ind2 = which(!is.na(ind))
df.p2v = df.p2v[ind2,]

df.p2v$gene_aa_str = df.NSFP$gene_aa_str[ind[ind2]]
df.p2v$gene = df.NSFP$genename[ind[ind2]]

# attach functional scores
include_LOF = T
df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_", ifelse(include_LOF, "", "noLOF_"), "042423.rds"))
ls.norm_score = readRDS("maveDB/norm_score_26genes_042423.rds")
df = ls.norm_score$BRCA1
df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)
ind = match(df.p2v$gene_aa_str, table = df.out$gene_aa_str)
temp = df.out[ind, c("norm_raw_score", "pos_score", "sub_score", "final_score", "final_score_lite")]
df.p2v = cbind(df.p2v, temp)
df.p2v = df.p2v %>%
  filter(!is.na(norm_raw_score) | !is.na(final_score)) %>% 
  group_by(patient) %>% 
  summarise(n_variants = length(variant), 
            variant_str = paste0(variant, collapse = ','), 
            AF = mean(AF),
            gene_aa_str = paste0(gene_aa_str, collapse = ','),
            gene = paste0(gene, collapse = ','),
            norm_raw_score = mean(norm_raw_score, na.rm=T),
            pos_score = mean(pos_score, na.rm=T),
            sub_score = mean(sub_score, na.rm=T),
            final_score = mean(final_score, na.rm=T),
            final_score_lite = mean(final_score_lite, na.rm=T)) %>%
  filter(n_variants == 1)

# attach cancer phenotypes
df.cancer = readRDS("UKBB/all_cancer_combined.rds")
df.cancer = df.cancer %>% 
  select(eid, sex, age, cancer=breast_cancer) %>%
  mutate(eid = as.character(eid))
ind = match(df.p2v$patient, table = df.cancer$eid)
df.p2v = cbind(df.p2v, df.cancer[ind,])
df.p2v = df.p2v %>% filter(sex == 0)
df.p2v$cancer = factor(df.p2v$cancer)

# score distribution between patient groups
ymax = 2
bw = 0.2
x_range = c(-3.5, 4.5)
plt1 = draw_patient_ggplot(df.p2v, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
plt2 = draw_patient_ggplot(df.p2v, score_type = "pos_score", score_type_name = "Positional component score", ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot(df.p2v, score_type = "sub_score", score_type_name = "Substitution component score", ymax = ymax, x_range = x_range, bw = bw)
plt4 = draw_patient_ggplot(df.p2v, score_type = "final_score_lite", score_type_name = "FUSE score\n(without predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
plt5 = draw_patient_ggplot(df.p2v, score_type = "final_score", score_type_name = "FUSE score\n(with predicted sites)", ymax = ymax, x_range = x_range, bw = bw)

# logit models
age_thres1 = 100
age_thres2 = 0

data = df.p2v %>% 
  mutate(func_score = norm_raw_score) %>% 
  filter(is.finite(func_score) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))
cancer_ct = length(which(data$cancer == 1))
if (cancer_ct>0){
  model1 <- glm(cancer~func_score, family="binomial", data=data)
  data$predicted = predict(model1, data, type="response")
  roc1 <- roc(data$cancer, data$predicted, auc=T, direction = "<")
}

data = df.p2v %>% 
  mutate(func_score = final_score) %>% 
  filter(is.finite(func_score) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))
cancer_ct = length(which(data$cancer == 1))
if (cancer_ct>0){
  model2 <- glm(cancer~func_score, family="binomial", data=data)
  data$predicted = predict(model2, data, type="response")
  roc2 <- roc(data$cancer, data$predicted, auc=T, direction = "<")
}

# plot ROC curves
roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)), 
               paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)))

plt6 = ggroc(list('Raw score' = roc1, 'Estimated score' = roc2)) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("Breast cancer classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.2), 
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))

plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
ggsave(filename = "~/Google Drive/My Drive/manuscript_resub/UKB_plots/breastcancer_allPlots_BRCA1_SGE.png", plot = plt, device = "png", width = 12, height = 8)

