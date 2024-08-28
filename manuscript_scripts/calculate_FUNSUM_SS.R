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

include_LOF = F

### generate SS-FUNSUM with proteinGym + maveDB scoresets
## generate overall FUNSUM
ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds")
df.all_sub_score = c()
for (scoreset_id in names(ls.combined_scoreset)){
  df = ls.combined_scoreset[[scoreset_id]]
  df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF, score_type = "norm_raw_score", sd_norm = F)
  df.all_sub_score = rbind(df.all_sub_score, cbind(scoreset_id, df.sub_score))
}
df.funsum_all = get_FUNSUM(df.all_sub_score, avg_method = "js", include_LOF = include_LOF) # get combined funsum
saveRDS(df.funsum_all, file = "./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")

## separate FUNSUM for each secondary structure in DSSP
ss_list = c("G" = "Helices", "H" = "Helices", "I" = "Helices", 
            "E" = "Strands", "B" = "Strands", 
            "T" = "Loops", "S" = "Loops", "C" = "Loops")
ls.sub_score = list()

ls.combined_scoreset = readRDS("./ProteinGym/combined_scoreset_051624.rds") # load raw scores in standard format
df.combined_scoreset_info = read_csv("./ProteinGym/combined_scoreset_info_111523.csv") # load protein gym dataset info sheet

for (scoreset_id in names(ls.combined_scoreset)){
  df = ls.combined_scoreset[[scoreset_id]]
  df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF, score_type = "norm_raw_score", sd_norm = F)
  df.sub_score$aapos_aaref = paste0(df.sub_score$aapos, df.sub_score$aaref)
  
  if (scoreset_id %in% df.combined_scoreset_info$scoreset_id){
    ind = which(df.combined_scoreset_info$scoreset_id == scoreset_id)
  } else{
    ind = which(df.combined_scoreset_info$gene_symbol == scoreset_id)[1]
  }
  
  if (length(ind)>0){
    dss_path = paste0("./alphaFold/dssp_out/", df.combined_scoreset_info$uniprot_id[ind], ".dss")
    if (file.exists(dss_path)){
      df.dssp = parse.dssp(file = dss_path, keepfiles = T)
      df.dssp$aapos_aaref = paste0(df.dssp$respdb, df.dssp$aa)
      ind = match(df.sub_score$aapos_aaref, table = df.dssp$aapos_aaref)
      
      match_rate = length(ind)/nrow(df.sub_score)
      if (match_rate < 0.9){
        print(paste0(scoreset_id, " DSSP match rate is ", match_rate))
      } else {
        df.sub_score$ss = df.dssp$ss[ind]
        
        for (ss in names(ss_list)){
          ind = which(df.sub_score$ss == ss)
          if (length(ind)>0){
            ls.sub_score[[ss]] = rbind(ls.sub_score[[ss]], cbind(scoreset_id, df.sub_score[ind,]))
          }
        }
      }
    } else {
      print(paste0("No .dss file can be found: ", scoreset_id))
    }
  } else {
    print(paste0("No gene symbol in df.combined_scoreset_info: ", scoreset_id))
  }
}

for (ss2 in unique(ss_list)){
  ss1 = names(ss_list[ss_list == ss2])
  for (ss in ss1){
    ls.sub_score[[ss2]] = rbind(ls.sub_score[[ss2]], ls.sub_score[[ss]])
  }
}

ls.funsum_ss = list()
for (ss in c(names(ss_list), unique(ss_list))){
  ls.funsum_ss[[ss]] = get_FUNSUM(ls.sub_score[[ss]], avg_method = "js", include_LOF = include_LOF ) # get combined funsum
}
saveRDS(ls.funsum_ss, file = "./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")




# FUNSUM-SS and BLOSUM62 heatmaps
df.funsum_all = readRDS("./ProteinGym/FUNSUM/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./ProteinGym/FUNSUM/funsum_SS_combined_noLOF_051624.rds")
data(BLOSUM62)
aa_list = colnames(df.funsum_all)
df.b62 = -BLOSUM62[aa_list, aa_list]

plt1 = draw_FUNSUM_heatmap(df.funsum_all, title_str = "FUNSUM (Overall)", show_diag = F)
plt2 = draw_FUNSUM_heatmap(ls.funsum_ss$Helices, title_str = "FUNSUM (Helices)", show_diag = F)
plt3 = draw_FUNSUM_heatmap(ls.funsum_ss$Strands, title_str = "FUNSUM (Strands)", show_diag = F)
plt4 = draw_FUNSUM_heatmap(ls.funsum_ss$Loops, title_str = "FUNSUM (Loops)", show_diag = F)
plt5 = draw_FUNSUM_heatmap(df.b62, title_str = "BLOSUM62", show_diag = F)

plt = (plt1 | plt2) / (plt3 | plt4) / (plt5 | plot_spacer())
ggsave(filename = "./manuscript/new_figures/Figure1_combined_heatmap.png", plot = plt, device = "png", 
       width = 10, height = 14)


# additional table to show secondary structure availability for all scoresets
df.all = read_csv("./manuscript/new_figures/supplementary_table/Table_S1__all_FUSE_scores.csv")
df.all_vars = unique(df.all %>% mutate(gene_aapos = paste0(gene, "---", aapos)) %>% select(gene, gene_aapos, DSSP_ss))

ss_list = c("G" = "Helices", "H" = "Helices", "I" = "Helices", 
            "E" = "Strands", "B" = "Strands", 
            "T" = "Loops", "S" = "Loops", "C" = "Loops")
df.all_vars$DSSP_ss2 = ss_list[df.all_vars$DSSP_ss]
df.all_vars$DSSP_ss2[is.na(df.all_vars$DSSP_ss2)] = "NA"
df.SS_ct = df.all_vars %>% group_by(DSSP_ss2) %>% summarise(SS_ct = length(DSSP_ss2))
write_csv(df.SS_ct, "./manuscript/new_figures/supplementary_table/AlphaFold_DSSP_secondary_structure_avaliability.csv")

