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


# ROC plot parameters
legend_pos = c(0.75, 0.25)
legend_margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "pt")
legend_text_margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "pt")

### final version of the UKBB figure
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
df.p2v$gene_aa_str = df.NSFP$gene_aa_str[ind]
df.p2v$gene = gene_id


# attach cancer phenotypes
df.cancer = readRDS("UKBB/all_cancer_combined.rds")
df.cancer = df.cancer %>% 
  select(eid, sex, age, cancer=breast_cancer) %>%
  mutate(eid = as.character(eid))

df.cancer2 = read_csv("~/Downloads/breast_cancer.csv")
df.cancer2 = df.cancer2 %>% select(eid, sex, age, cancer=disease)
df.cancer3 = read_csv("~/Downloads/ovarian_cancer.csv")
df.cancer3 = df.cancer3 %>% select(eid, sex, age, cancer=disease)

df.cancer = rbind(df.cancer, df.cancer2, df.cancer3)
df.cancer = df.cancer %>% group_by(eid) %>%
  summarise(sex = unique(sex), age = mean(age), cancer = sum(cancer))
df.cancer$cancer[df.cancer$cancer>0] = 1


ind = match(df.p2v$patient, table = df.cancer$eid)
df.p2v = cbind(df.p2v, df.cancer[ind,])
df.p2v = df.p2v %>% filter(sex == 0)
df.p2v$cancer = factor(df.p2v$cancer)


# attach functional scores
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")

# de-noise DMS data
df = ls.combined_scoreset[[gene_id]]
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df$gene = gene_id
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
df.out$final_score_ss_lite = df.out$final_score_ss
df.out$final_score_ss_lite[is.na(df.out$raw_score)] = NA

# attach am score
df.am = read_csv("./ProteinGym/combined_48genes_am_scores.csv")
ind = match(df.out$gene_aa_str, table = df.am$gene_aa_str)
df.out$am_score = df.am$am_pathogenicity[ind]

## draw plots
df.p2v_2 = merge(df.p2v, df.out)

# score distribution between patient groups
ymax = 3
bw = 0.1
x_range = c(-2, 2.5)
plt1 = draw_patient_ggplot_v2(df.p2v_2, score_type = "norm_raw_score", score_type_name = "Normalized original score", legend_labels = c("No cancer", "Breast cancer"), ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.26, 0.87))
plt2 = draw_patient_ggplot_v2(df.p2v_2, score_type = "final_score_ss_lite", score_type_name = "FUSE score\n(without predicted sites)", legend_labels = c("No cancer", "Breast cancer"), ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot_v2(df.p2v_2, score_type = "final_score_ss", score_type_name = "FUSE score\n(with predicted sites)", legend_labels = c("No cancer", "Breast cancer"), ymax = ymax, x_range = x_range, bw = bw)

# ROC analysis with logistic regression model
age_thres1 = 100
age_thres2 = 0

df.p2v_3 = df.p2v_2 %>% 
  filter(!is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(am_score) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))

df.train = df.p2v_3
df.test = df.p2v_3

cancer_ct = length(which(df.train$cancer == 1))
if (cancer_ct>0){
  model <- glm(cancer~norm_raw_score, family="binomial", data=df.train)
  df.test$model1_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~pos_score+sub_score_ss+acc, family="binomial", data=df.train)
  df.test$model2_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~pos_score+sub_score_ss+acc+am_score, family="binomial", data=df.train)
  df.test$model3_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~am_score, family="binomial", data=df.train)
  df.test$model4_predicted = predict(model, df.test, type="response")
}

roc1 <- roc(df.test$cancer, df.test$model1_predicted, auc=T, direction = "<")
roc2 <- roc(df.test$cancer, df.test$model2_predicted, auc=T, direction = "<")
roc3 <- roc(df.test$cancer, df.test$model3_predicted, auc=T, direction = "<")
roc4 <- roc(df.test$cancer, df.test$model4_predicted, auc=T, direction = "<")
roc_all = list(roc1,roc2,roc3,roc4)

roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
               paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)),
               paste0("FUSE + AM, AUC=", round(roc3$auc, digits = 2)),
               paste0("AM, AUC=", round(roc4$auc, digits = 2)))

# plot ROC curves
plt4 = ggroc(roc_all) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("HBOC cancer classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = legend_pos, 
        legend.text = element_text(size = 8, margin = legend_text_margin),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2),
        legend.margin = legend_margin)


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
df.p2v = df.p2v %>% group_by(patient) %>% 
  summarise(n_variants = length(variant), variant_str = paste0(variant, collapse = ',')) %>%
  filter(n_variants == 1)
df.p2v = merge(df.p2v, df.NSFP)
df.p2v$gene = df.p2v$genename

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

# attach functional scores
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")

# de-noise DMS data
df = ls.combined_scoreset[[gene_id]]
ind = which(df.combined_scoreset_info$gene_symbol == gene_id)[1]
df$gene = gene_id
df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                           dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                           include_LOF = include_LOF, show_func_class = T)
df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
df.out$final_score_ss_lite = df.out$final_score_ss
df.out$final_score_ss_lite[is.na(df.out$raw_score)] = NA

# attach am score
df.am = read_csv("./ProteinGym/combined_48genes_am_scores.csv")
ind = match(df.out$gene_aa_str, table = df.am$gene_aa_str)
df.out$am_score = df.am$am_pathogenicity[ind]


## draw plots
ind = match(df.p2v$gene_aa_str, table = df.out$gene_aa_str)
temp = df.out[ind, c("norm_raw_score", "pos_score", "sub_score_ss", "final_score", "final_score_ss", "final_score_ss_lite")]
df.p2v_2 = merge(df.p2v, df.out)

# ROC analysis with logistic regression model
age_thres1 = 100
age_thres2 = 0

df.p2v_3 = df.p2v_2 %>% 
  filter(!is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(am_score) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))
df.train = df.p2v_3
df.test = df.p2v_3

cancer_ct = length(which(df.train$cancer == 1))
if (cancer_ct>0){
  model <- glm(cancer~norm_raw_score, family="binomial", data=df.train)
  df.test$model1_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~pos_score+sub_score+sub_score_ss+acc, family="binomial", data=df.train)
  df.test$model2_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~pos_score+sub_score+sub_score_ss+acc+am_score, family="binomial", data=df.train)
  df.test$model3_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~am_score, family="binomial", data=df.train)
  df.test$model4_predicted = predict(model, df.test, type="response")
}

roc1 <- roc(df.test$cancer, df.test$model1_predicted, auc=T, direction = "<")
roc2 <- roc(df.test$cancer, df.test$model2_predicted, auc=T, direction = "<")
roc3 <- roc(df.test$cancer, df.test$model3_predicted, auc=T, direction = "<")
roc4 <- roc(df.test$cancer, df.test$model4_predicted, auc=T, direction = "<")
roc_all = list(roc1,roc2,roc3,roc4)

roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
               paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)),
               paste0("FUSE + AM, AUC=", round(roc3$auc, digits = 2)),
               paste0("AM, AUC=", round(roc4$auc, digits = 2)))
plt5 = ggroc(roc_all) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("LFS-related cancer classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = legend_pos, 
        legend.text = element_text(size = 8, margin = legend_text_margin),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2),
        legend.margin = legend_margin)




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

df.p2v_2 = merge(df.p2v, df.out)

# ROC analysis with logistic regression model
df.p2v_3 = df.p2v_2 %>% 
  filter(!is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(am_score) & !is.na(cad))

df.train = df.p2v_3
df.test = df.p2v_3

cad_ct = length(which(df.train$cad == 1))
if (cad_ct>0){
  model <- glm(cad~norm_raw_score, family="binomial", data=df.train)
  df.test$model1_predicted = predict(model, df.test, type="response")
  
  model <- glm(cad~pos_score+sub_score_ss+acc, family="binomial", data=df.train)
  df.test$model2_predicted = predict(model, df.test, type="response")
  
  model <- glm(cad~pos_score+sub_score_ss+acc+am_score, family="binomial", data=df.train)
  df.test$model3_predicted = predict(model, df.test, type="response")
  
  model <- glm(cad~am_score, family="binomial", data=df.train)
  df.test$model4_predicted = predict(model, df.test, type="response")
}

roc1 <- roc(df.test$cad, df.test$model1_predicted, auc=T, direction = "<")
roc2 <- roc(df.test$cad, df.test$model2_predicted, auc=T, direction = "<")
roc3 <- roc(df.test$cad, df.test$model3_predicted, auc=T, direction = "<")
roc4 <- roc(df.test$cad, df.test$model4_predicted, auc=T, direction = "<")
roc_all = list(roc1,roc2,roc3,roc4)

roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
               paste0("FUSE score, AUC=", round(roc2$auc, digits = 2)),
               paste0("FUSE + AM, AUC=", round(roc3$auc, digits = 2)),
               paste0("AM, AUC=", round(roc4$auc, digits = 2)))

# plot ROC curves
plt6 = ggroc(roc_all) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("CAD classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = legend_pos, 
        legend.text = element_text(size = 8, margin = legend_text_margin),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2),
        legend.margin = legend_margin)

plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
ggsave(paste0("./manuscript/new_figures/UKBB_plots/UKBB_allPlots_addVineel3.png"), plt, device = "png", width = 12, height = 8)



















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
df.p2v = df.p2v %>% group_by(patient) %>% 
  summarise(n_variants = length(variant), variant_str = paste0(variant, collapse = ',')) %>%
  filter(n_variants == 1)
df.p2v = merge(df.p2v, df.NSFP)
df.p2v$gene = df.p2v$genename

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

# attach functional scores
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_120523.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_120523.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_112923.rds")

df.out_all = c()
for (ind in which(df.combined_scoreset_info$gene_symbol == gene_id)){
  scoreset_id = df.combined_scoreset_info$scoreset_id[ind]
  df = ls.combined_scoreset[[scoreset_id]]
  df$gene = df.combined_scoreset_info$gene_symbol[ind]
  df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                             dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                             include_LOF = include_LOF, show_func_class = T)
  df.out$scoreset_id = scoreset_id
  df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
  df.out_all = rbind(df.out_all, df.out)
}
df.out_all$final_score_ss_lite = df.out_all$final_score_ss
df.out_all$final_score_ss_lite[is.na(df.out_all$raw_score)] = NA

## draw plots
scoreset_id = "P53_HUMAN_Giacomelli_NULL_Nutlin_2018"
ind = which(df.out_all$scoreset_id == scoreset_id)
df.out = df.out_all[ind,]

ind = match(df.p2v$gene_aa_str, table = df.out$gene_aa_str)
temp = df.out[ind, c("norm_raw_score", "pos_score", "sub_score_ss", "final_score", "final_score_ss", "final_score_ss_lite")]
df.p2v_2 = merge(df.p2v, df.out)

# score distribution between patient groups
ymax = 2
bw = 0.2
x_range = c(-3.5, 4.5)
plt1 = draw_patient_ggplot(df.p2v_2, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
plt2 = draw_patient_ggplot(df.p2v_2, score_type = "pos_score", score_type_name = "Positional component score", ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot(df.p2v_2, score_type = "sub_score_ss", score_type_name = "Substitution component score", ymax = ymax, x_range = x_range, bw = bw)
plt4 = draw_patient_ggplot(df.p2v_2, score_type = "final_score_ss_lite", score_type_name = "SS-FUSE score\n(without predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
plt5 = draw_patient_ggplot(df.p2v_2, score_type = "final_score_ss", score_type_name = "SS-FUSE score\n(with predicted sites)", ymax = ymax, x_range = x_range, bw = bw)

# ROC analysis with logistic regression model
age_thres1 = 100
age_thres2 = 0

df.p2v_3 = df.p2v_2 %>% 
  filter(!is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))
df.train = df.p2v_3
df.test = df.p2v_3

cancer_ct = length(which(df.train$cancer == 1))
if (cancer_ct>0){
  model <- glm(cancer~norm_raw_score, family="binomial", data=df.train)
  df.test$model1_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~pos_score+sub_score+sub_score_ss+acc, family="binomial", data=df.train)
  df.test$model2_predicted = predict(model, df.test, type="response")
}

# plot ROC curves
roc1 <- roc(df.test$cancer, df.test$model1_predicted, auc=T, direction = "<")
roc2 <- roc(df.test$cancer, df.test$model2_predicted, auc=T, direction = "<")
roc_all = list(roc1,roc2)

roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
               paste0("SS-FUSE, AUC=", round(roc2$auc, digits = 2)))

plt6 = ggroc(roc_all) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("LFS-associated cancer classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = legend_pos, 
        legend.text = element_text(size = 8, margin = legend_text_margin),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))

plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
ggsave(paste0("./manuscript/new_figures/LFScancer_allPlots_TP53_DMS_", scoreset_id, ".png"), plt, device = "png", width = 12, height = 8)




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
df.p2v$gene_aa_str = df.NSFP$gene_aa_str[ind]
df.p2v$gene = gene_id


# attach cancer phenotypes
df.cancer = readRDS("UKBB/all_cancer_combined.rds")
df.cancer = df.cancer %>% 
  select(eid, sex, age, cancer=breast_cancer) %>%
  mutate(eid = as.character(eid))
ind = match(df.p2v$patient, table = df.cancer$eid)
df.p2v = cbind(df.p2v, df.cancer[ind,])
df.p2v = df.p2v %>% filter(sex == 0)
df.p2v$cancer = factor(df.p2v$cancer)


# attach functional scores
include_LOF = F
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_120523.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_120523.rds")
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv")
df.combined_scoreset_info = df.combined_scoreset_info %>% filter(!is.na(gene_symbol) & !is.na(uniprot_id))
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_112923.rds")

df.out_all = c()
for (ind in which(df.combined_scoreset_info$gene_symbol == gene_id)){
  scoreset_id = df.combined_scoreset_info$scoreset_id[ind]
  df = ls.combined_scoreset[[scoreset_id]]
  df$gene = df.combined_scoreset_info$gene_symbol[ind]
  df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                             dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss"), 
                             include_LOF = include_LOF, show_func_class = T)
  df.out$scoreset_id = scoreset_id
  df.out = df.out %>% filter((aaref != aaalt) & (aaref != "*") & (aaalt != "*"))
  df.out_all = rbind(df.out_all, df.out)
}
df.out_all$final_score_ss_lite = df.out_all$final_score_ss
df.out_all$final_score_ss_lite[is.na(df.out_all$raw_score)] = NA


## draw plots
scoreset_id = "BRCA1_HUMAN_Findlay_2018"
ind = which(df.out_all$scoreset_id == scoreset_id)
df.out = df.out_all[ind,]

df.p2v_2 = merge(df.p2v, df.out)

# score distribution between patient groups
ymax = 2
bw = 0.2
x_range = c(-3.5, 4.5)
plt1 = draw_patient_ggplot(df.p2v_2, score_type = "norm_raw_score", score_type_name = "Normalized original score", ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
plt2 = draw_patient_ggplot(df.p2v_2, score_type = "pos_score", score_type_name = "Positional component score", ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot(df.p2v_2, score_type = "sub_score_ss", score_type_name = "Substitution component score", ymax = ymax, x_range = x_range, bw = bw)
plt4 = draw_patient_ggplot(df.p2v_2, score_type = "final_score_ss_lite", score_type_name = "SS-FUSE score\n(without predicted sites)", ymax = ymax, x_range = x_range, bw = bw)
plt5 = draw_patient_ggplot(df.p2v_2, score_type = "final_score_ss", score_type_name = "SS-FUSE score\n(with predicted sites)", ymax = ymax, x_range = x_range, bw = bw)

# ROC analysis with logistic regression model
age_thres1 = 100
age_thres2 = 0

df.p2v_3 = df.p2v_2 %>% 
  filter(!is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(cancer)) %>%
  filter((age<age_thres1 & cancer==1) | (age>age_thres2 & cancer==0))

df.train = df.p2v_3
df.test = df.p2v_3

cancer_ct = length(which(df.train$cancer == 1))
if (cancer_ct>0){
  model <- glm(cancer~norm_raw_score, family="binomial", data=df.train)
  df.test$model1_predicted = predict(model, df.test, type="response")
  
  model <- glm(cancer~pos_score+sub_score_ss+acc, family="binomial", data=df.train)
  df.test$model2_predicted = predict(model, df.test, type="response")
}

roc1 <- roc(df.test$cancer, df.test$model1_predicted, auc=T, direction = "<")
roc2 <- roc(df.test$cancer, df.test$model2_predicted, auc=T, direction = "<")
roc_all = list(roc1,roc2)

roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
               paste0("SS-FUSE, AUC=", round(roc2$auc, digits = 2)))

# plot ROC curves
plt6 = ggroc(roc_all) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("Breast cancer classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = legend_pos, 
        legend.text = element_text(size = 8, margin = legend_text_margin),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2))

plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
ggsave(paste0("./manuscript/new_figures/Breastcancer_allPlots_BRCA1_SGE_", scoreset_id, ".png"), plt, device = "png", width = 12, height = 8)






### LDLR base editing data
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
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_120523.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_120523.rds")

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
df.p2v_2 = merge(df.p2v, df.out)

# score distribution between patient groups
ymax = 2
bw = 0.2
x_range = c(-3.5, 4.5)
plt1 = draw_patient_ggplot_v2(df.p2v_2, score_type = "norm_raw_score", score_type_name = "Normalized original score", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw, show_legend = T, legend_pos = c(0.18, 0.85))
plt2 = draw_patient_ggplot_v2(df.p2v_2, score_type = "pos_score", score_type_name = "Positional component score", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw)
plt3 = draw_patient_ggplot_v2(df.p2v_2, score_type = "sub_score_ss", score_type_name = "Substitution component score", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw)
plt4 = draw_patient_ggplot_v2(df.p2v_2, score_type = "final_score_ss_lite", score_type_name = "FUSE score\n(without predicted sites)", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw)
plt5 = draw_patient_ggplot_v2(df.p2v_2, score_type = "final_score_ss", score_type_name = "FUSE score\n(with predicted sites)", phenotype_col = "cad", legend_labels = c("No CAD", "CAD"), ymax = ymax, x_range = x_range, bw = bw)

# ROC analysis with logistic regression model
df.p2v_3 = df.p2v_2 %>% 
  filter(!is.na(norm_raw_score) & !is.na(final_score) & !is.na(pos_score) & !is.na(sub_score) & !is.na(sub_score_ss) & !is.na(acc) & !is.na(cad))

df.train = df.p2v_3
df.test = df.p2v_3

cad_ct = length(which(df.train$cad == 1))
if (cad_ct>0){
  model <- glm(cad~norm_raw_score, family="binomial", data=df.train)
  df.test$model1_predicted = predict(model, df.test, type="response")
  
  model <- glm(cad~pos_score+sub_score_ss+acc, family="binomial", data=df.train)
  df.test$model2_predicted = predict(model, df.test, type="response")
}

roc1 <- roc(df.test$cad, df.test$model1_predicted, auc=T, direction = "<")
roc2 <- roc(df.test$cad, df.test$model2_predicted, auc=T, direction = "<")
roc_all = list(roc1,roc2)

roc_labels = c(paste0("Raw score, AUC=", round(roc1$auc, digits = 2)),
               paste0("FUSE, AUC=", round(roc2$auc, digits = 2)))

# plot ROC curves
plt6 = ggroc(roc_all) + 
  annotate("segment",x = 1, xend = 0, y = 0, yend = 1, color="lightgrey", linetype="dashed") +
  scale_color_discrete(name = NULL, labels = roc_labels) + 
  ggtitle("CAD cases classification") + xlab("Specificity") + ylab("Sensitivity") + 
  theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = legend_pos, 
        legend.text = element_text(size = 8, margin = legend_text_margin),
        legend.background = element_rect(fill="white", linetype=1, color="lightgrey", linewidth=0.2),
        legend.margin = legend_margin)

plt = (plt1|plt2|plt3)/(plt4|plt5|plt6) + plot_annotation(tag_levels = 'A')
ggsave(paste0("./manuscript/new_figures/CAD_allPlots_LDLR_BE.png"), plt, device = "png", width = 12, height = 8)


