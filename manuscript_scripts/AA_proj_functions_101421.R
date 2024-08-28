library(openxlsx)
library(tidyverse)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

get_KS_pval <-  function(x, y, alternative = "greater"){
  require(reticulate)
  use_condaenv("base")
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  
  if (length(x)>1 & length(y)>1){
    py_run_string(code = "from scipy.stats import ks_2samp")
    res = py$ks_2samp(r_to_py(x), r_to_py(y), alternative = alternative)
    return(res$pvalue)
  } else {
    return(NA)
  }
}

get_MW_pval <-  function(x, y, alternative = "less"){
  require(reticulate)
  use_condaenv("base")
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  
  if (length(x)>1 & length(y)>1){
    py_run_string(code = "from scipy.stats import mannwhitneyu")
    res = py$mannwhitneyu(r_to_py(x), r_to_py(y), alternative = alternative)
    return(res$pvalue)
  } else {
    return(NA)
  }
}

# JS estimator
get_js <- function(m){
  mbar = mean(colMeans(m, na.rm = T), na.rm = T) # global mean
  mu0 = colMeans(m, na.rm = T) # column means
  s2 = var(as.vector(m), na.rm = T)/(nrow(m)*ncol(m)-length(which(is.na(as.vector(m)) == T))) # global variance
  cval = 1 - (ncol(m)-2)*s2/sum((mu0 - mbar)^2, na.rm = T) # adjusted value
  js_est = mbar + cval*(mu0 - mbar)
  return(js_est)
}

get_score_matrix <- function(df, sd_norm = T, invert_score=F, include_LOF=F, infer_score=F, na.replace=NA){
  aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  }
  df.out = data.frame(aapos = unique(df$aapos), aaref = NA)
  temp = matrix(NA, nrow = nrow(df.out), ncol = length(aa_list))
  colnames(temp) <- aa_list
  
  for (i in 1:nrow(df.out)){
    ind = which(df$aapos == df.out$aapos[i])
    df.out$aaref[i] = unique(df$aaref[ind])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$raw_score[ind[ind2]]
  }
  
  if(invert_score){
    temp = -temp
  }
  
  # normalize the data by sd
  if (sd_norm){
    temp = (temp-median(temp, na.rm = T))/sd(temp, na.rm = T)
  }
  
  if (infer_score){
    df.funsum_all = read.csv("FUNSUM_all_111221.csv")
    if (include_LOF){
      df.funsum_all = read.csv("FUNSUM_all_LOF_122021.csv")
    }
    rownames(df.funsum_all) = df.funsum_all[,1]
    df.funsum_all = df.funsum_all[,-1]
    if (include_LOF){
      colnames(df.funsum_all) = aa_list
    }
    
    for (i in 1:nrow(temp)){
      aaref = df.out$aaref[i]
      ind = which(is.na(temp[i,]))
      ind2 = which(!is.na(temp[i,]))
      if (length(ind)>0 & length(ind2)>0){
        temp2 = temp[i,ind2] - df.funsum_all[aaref, aa_list[ind2]]
        pos_mean = mean(as.numeric(temp2), na.rm=T)
        temp[i,ind] = as.numeric(pos_mean + df.funsum_all[aaref, aa_list[ind]])
      }
    }
  }
  
  temp[is.na(temp)] = na.replace
  df.out = cbind(df.out, temp)
  return(df.out)
}

pos_norm <- function(df.score_mat, pos_mean_method, include_LOF=F){
  aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  }
  
  df.out = data.frame(aapos = df.score_mat$aapos, aaref = df.score_mat$aaref, pos_mean = NA)
  temp = df.score_mat[,aa_list]
  
  if (pos_mean_method == "mean"){
    df.out$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.out$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.out$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "FUNSUM-HMM"){
    df.funsum_low = read.csv("FUNSUM_low_011422.csv")
    rownames(df.funsum_low) = df.funsum_low[,1]
    df.funsum_low = df.funsum_low[,-1]
    df.funsum_high = read.csv("FUNSUM_high_011422.csv")
    rownames(df.funsum_high) = df.funsum_high[,1]
    df.funsum_high = df.funsum_high[,-1]
    
    pos_mean_all = c()
    for (i in 1:nrow(temp)){
      aaref = df.out$aaref[i]
      pos_mean = NA
      ind = which(!is.na(temp[i,]))
      if (df.score_mat$state[i] == 1){
        temp2 = temp[i,ind] - df.funsum_low[aaref, aa_list[ind]]
      } else {
        temp2 = temp[i,ind] - df.funsum_high[aaref, aa_list[ind]]
      }
      
      pos_mean = mean(as.numeric(temp2))
      pos_mean_all = c(pos_mean_all, pos_mean)
    }
    
    df.out$pos_mean = pos_mean_all
  }
  
  temp = temp - df.out$pos_mean
  df.out = cbind(df.out, temp)
  return(df.out)
}

check_state <- function(df.score_mat, est.states){
  ind = which(est.states$state == 1)
  ind2 = which(est.states$state == 2)
  aa_list = colnames(df.score_mat)[-c(1:2)]
  if (mean(unlist(df.score_mat[ind, aa_list]), na.rm=T) > mean(unlist(df.score_mat[ind2, aa_list]), na.rm=T)){
    est.states$state[ind] = 2
    est.states$state[ind2] = 1
    temp = est.states$S1
    est.states$S1 = est.states$S2
    est.states$S2 = temp
  }
  
  return(est.states)
}


get_norm_score <-  function(df, pos_mean_method, pos_norm = T, sd_norm = T, invert_score=F, include_LOF=F, infer_score=F, score_type="raw_score"){
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
  }
  
  ind = which((df$aaref %in% aa_list) & (df$aaalt %in% aa_list))
  df = df[ind,]
  
  df.out = data.frame(aapos = unique(df$aapos), aaref = NA, pos_mean = NA)
  temp = matrix(NA, nrow = nrow(df.out), ncol = length(aa_list))
  colnames(temp) <- aa_list
  
  for (i in 1:nrow(df.out)){
    ind = which(df$aapos == df.out$aapos[i])
    df.out$aaref[i] = unique(df$aaref[ind])[1]
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df[[score_type]][ind[ind2]]
  }
  
  if(invert_score){
    temp = -temp
  }
  
  # normalize the data by sd
  if (sd_norm){
    temp = (temp-median(temp, na.rm = T))/sd(temp, na.rm = T)
  }
  
  if (pos_mean_method == "mean"){
    df.out$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.out$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.out$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "funsum"){
    df.funsum_all = readRDS("./shiny_app/DMS_denoiser/funsum_081722.rds")
    if (!include_LOF){
      df.funsum_all = df.funsum_all[aa_list, aa_list]
    }
    
    pos_mean_all = c()
    for (i in 1:nrow(temp)){
      aaref = df.out$aaref[i]
      pos_mean = NA
      ind = which(!is.na(temp[i,]))
      temp2 = temp[i,ind] - df.funsum_all[aaref, aa_list[ind]]
      pos_mean = mean(as.numeric(temp2))
      pos_mean_all = c(pos_mean_all, pos_mean)
    }
    
    df.out$pos_mean = pos_mean_all
  }
  
  if (pos_norm){
    temp = temp - df.out$pos_mean
  }
  
  if (infer_score){
    temp = df.out$pos_mean + df.funsum_all[df.out$aaref,aa_list]
  }
  
  df.out = cbind(df.out, temp)
  return(df.out)
}

pivot_norm_score <-  function(df, pos_mean_method, include_LOF=F, score_type="norm_raw_score", df.funsum_all=NULL){
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
  }
  
  ind = which((df$aaref %in% aa_list) & (df$aaalt %in% aa_list))
  df = df[ind,]
  
  df.out = data.frame(aapos = unique(df$aapos), aaref = NA, pos_mean = NA)
  temp = matrix(NA, nrow = nrow(df.out), ncol = length(aa_list))
  colnames(temp) <- aa_list
  
  for (i in 1:nrow(df.out)){
    ind = which(df$aapos == df.out$aapos[i])
    df.out$aaref[i] = unique(df$aaref[ind])[1]
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df[[score_type]][ind[ind2]]
  }
  
  if (pos_mean_method == "mean"){
    df.out$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.out$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.out$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "funsum"){
    if (is.null(df.funsum_all)){
      df.funsum_all = readRDS(paste0("./shiny_app/DMS_denoiser/funsum_maveDB_042423.rds"))
    }
    if (!include_LOF){
      df.funsum_all = df.funsum_all[aa_list, aa_list]
    }
    
    pos_mean_all = c()
    for (i in 1:nrow(temp)){
      aaref = df.out$aaref[i]
      pos_mean = NA
      ind = which(!is.na(temp[i,]))
      temp2 = temp[i,ind] - df.funsum_all[aaref, aa_list[ind]]
      pos_mean = mean(as.numeric(temp2))
      pos_mean_all = c(pos_mean_all, pos_mean)
    }
    
    df.out$pos_mean = pos_mean_all
  }
  
  df.out = cbind(df.out, temp)
  return(df.out)
}



get_norm_score_v2 <-  function(df, pos_mean_method, pos_norm = T){
  df = df[!is.na(df$raw_score),] # remove rows with empty scores
  
  ## normalization with synonymous and LOF substitutions
  ind = which(df$aaref == df$aaalt)
  ind2 = which(df$aaalt == "*")
  if (length(ind)>0 & length(ind2)>0){
    m1 = median(df$raw_score[ind])
    m2 = median(df$raw_score[ind2])
  
    if (m2 < m1){
      df$raw_score = -df$raw_score
    }
    m1 = median(df$raw_score[ind])
    m2 = median(df$raw_score[ind2])
    
    df$raw_score = (df$raw_score - m1)/(m2-m1)
  }
  
  aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  df.out = data.frame(aapos = unique(df$aapos), aaref = NA, pos_mean = NA)
  temp = matrix(NA, nrow = nrow(df.out), ncol = length(aa_list))
  colnames(temp) <- aa_list
  
  for (i in 1:nrow(df.out)){
    ind = which(df$aapos == df.out$aapos[i])
    df.out$aaref[i] = unique(df$aaref[ind])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$raw_score[ind[ind2]]
  }
  
  if (pos_mean_method == "mean"){
    df.out$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.out$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.out$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "FUNSUM"){
    df.funsum_all = readRDS("./shiny_app/DMS_denoiser/funsum_031622.rds")
    
    pos_mean_all = c()
    for (i in 1:nrow(temp)){
      aaref = df.out$aaref[i]
      pos_mean = NA
      ind = which(!is.na(temp[i,]))
      temp2 = temp[i,ind] - df.funsum_all[aaref, aa_list[ind]]
      pos_mean = mean(as.numeric(temp2))
      pos_mean_all = c(pos_mean_all, pos_mean)
    }
    
    df.out$pos_mean = pos_mean_all
  }
  
  if (pos_norm){
    temp = temp - df.out$pos_mean
  }
  
  df.out = cbind(df.out, temp)
  return(df.out)
}

sample_rand <- function(df.raw, edit_num, include_LOF=F){
  if (!include_LOF){
    df.raw = df.raw[(!df.raw$aaref %in% "*"),]
    df.raw = df.raw[(!df.raw$aaalt %in% "*"),]
  }
  
  if (edit_num > nrow(df.raw)){
    stop(paste0("edit_num cannot exceed ", nrow(df.raw)))
  }
  df.raw$gene_aapos = paste0(df.raw$gene, "---", df.raw$aapos)
  uniq_gene_aapos = unique(df.raw$gene_aapos)
  gene_aapos_len = length(uniq_gene_aapos)
  df.raw$iteration = 1
  
  df.out = c()
  while (edit_num > gene_aapos_len){
    ind_selected = c()
    past_aapos = c()
    for (gene_aapos in uniq_gene_aapos){
      ind = which(df.raw$gene_aapos == gene_aapos)
      if (length(ind)>1){
        ind2 = sample(ind, size = 1)
      } else {
        ind2 = ind
      }
      ind_selected = c(ind_selected, ind2)
      past_aapos = c(past_aapos, gene_aapos)
    }
    
    df.out = rbind(df.out, df.raw[ind_selected,])
    df.raw = df.raw[-ind_selected,]
    df.raw$iteration = df.raw$iteration+1
    edit_num = edit_num - gene_aapos_len
    uniq_gene_aapos = unique(df.raw$gene_aapos)
    gene_aapos_len = length(uniq_gene_aapos)
  }
  
  uniq_gene_aapos = sample(uniq_gene_aapos, edit_num, replace = F)
  ind_selected = c()
  for (gene_aapos in uniq_gene_aapos){
    ind = which(df.raw$gene_aapos == gene_aapos)
    if (length(ind)>1){
      ind2 = sample(ind, size = 1)
    } else {
      ind2 = ind
    }
    ind_selected = c(ind_selected, ind2)
  }
  df.out = rbind(df.out, df.raw[ind_selected,])
  return(df.out)
}

# sample_rank_v2 <- function(df.raw, edit_num, df.aa_rank_linear, include_LOF=F){
#   if (!include_LOF){
#     df.raw = df.raw[(!df.raw$aaref %in% "*"),]
#     df.raw = df.raw[(!df.raw$aaalt %in% "*"),]
#   }
#   
#   if (edit_num > nrow(df.raw)){
#     stop(paste0("edit_num cannot exceed ", nrow(df.raw)))
#   }
#   df.raw$gene_aapos = paste0(df.raw$gene, "---", df.raw$aapos)
#   uniq_gene_aapos = unique(df.raw$gene_aapos)
#   gene_aapos_len = length(uniq_gene_aapos)
#   df.raw$iteration = 1
#   
#   df.out = c()
#   while (edit_num > gene_aapos_len){
#     ind_selected = c()
#     past_aapos = c()
#     for (gene_aapos in uniq_gene_aapos){
#       ind = which(df.raw$gene_aapos == gene_aapos)
#       if (length(ind)>1){
#         aaref = unique(df.raw$aaref[ind])
#         aaalts = df.raw$aaalt[ind]
#         ind2 = which(df.aa_rank_linear$aaref == aaref & df.aa_rank_linear$aaalt %in% aaalts)
#         ind3 = order(df.aa_rank_linear$rank[ind2])
#         ind4 = ind[which(df.raw$aaalt[ind] == df.aa_rank_linear$aaalt[ind2][ind3][1])]
#       } else {
#         ind4 = ind
#       }
#       ind_selected = c(ind_selected, ind4)
#       past_aapos = c(past_aapos, gene_aapos)
#     }
#     
#     df.out = rbind(df.out, df.raw[ind_selected,])
#     df.raw = df.raw[-ind_selected,]
#     df.raw$iteration = df.raw$iteration+1
#     edit_num = edit_num - gene_aapos_len
#     uniq_gene_aapos = unique(df.raw$gene_aapos)
#     gene_aapos_len = length(uniq_gene_aapos)
#   }
#   
#   temp = data.frame(aapos = uniq_gene_aapos, aaref = NA, best_aaalt = NA, best_rank = NA, best_aaalt_ind = NA)
#   for (i in 1:nrow(temp)){
#     ind = which(df.raw$gene_aapos == temp$aapos[i])
#     aaref = unique(df.raw$aaref[ind])
#     aaalts = df.raw$aaalt[ind]
#     ind2 = which(df.aa_rank_linear$aaref == aaref & df.aa_rank_linear$aaalt %in% aaalts)
#     ind3 = order(df.aa_rank_linear$rank[ind2])
#     
#     temp$aaref[i] = aaref
#     temp$best_aaalt[i] = df.aa_rank_linear$aaalt[ind2][ind3][1]
#     temp$best_rank[i] = df.aa_rank_linear$rank[ind2][ind3][1]
#     temp$best_aaalt_ind[i] = ind[which(df.raw$aaalt[ind] == df.aa_rank_linear$aaalt[ind2][ind3][1])]
#   }
#   
#   ind_selected = temp$best_aaalt_ind[order(temp$best_rank)[1:edit_num]]
#   df.out = rbind(df.out, df.raw[ind_selected,])
#   return(df.out)
# }

sample_rank <- function(df.raw, edit_num, df.aa_rank, include_LOF=F){
  if (!include_LOF){
    df.raw = df.raw[(!df.raw$aaref %in% "*"),]
    df.raw = df.raw[(!df.raw$aaalt %in% "*"),]
  }
  
  if (edit_num > nrow(df.raw)){
    stop(paste0("edit_num cannot exceed ", nrow(df.raw)))
  }
  df.raw$gene_aapos = paste0(df.raw$gene, "---", df.raw$aapos)
  uniq_gene_aapos = unique(df.raw$gene_aapos)
  gene_aapos_len = length(uniq_gene_aapos)
  df.raw$iteration = 1
  
  df.out = c()
  while (edit_num > gene_aapos_len){
    ind_selected = c()
    past_aapos = c()
    for (gene_aapos in uniq_gene_aapos){
      ind = which(df.raw$gene_aapos == gene_aapos)
      if (length(ind)>1){
        aaref = unique(df.raw$aaref[ind])
        aaalts = df.raw$aaalt[ind]
        temp = df.aa_rank[aaref,]
        temp = temp[temp %in% aaalts]
        ind2 = ind[which(df.raw$aaalt[ind] == temp[1])]
      } else {
        ind2 = ind
      }
      ind_selected = c(ind_selected, ind2)
      past_aapos = c(past_aapos, gene_aapos)
    }
    
    df.out = rbind(df.out, df.raw[ind_selected,])
    df.raw = df.raw[-ind_selected,]
    df.raw$iteration = df.raw$iteration+1
    edit_num = edit_num - gene_aapos_len
    uniq_gene_aapos = unique(df.raw$gene_aapos)
    gene_aapos_len = length(uniq_gene_aapos)
  }
  
  uniq_gene_aapos = sample(uniq_gene_aapos, edit_num, replace = F)
  ind_selected = c()
  for (gene_aapos in uniq_gene_aapos){
    ind = which(df.raw$gene_aapos == gene_aapos)
    if (length(ind)>1){
      aaref = unique(df.raw$aaref[ind])
      aaalts = df.raw$aaalt[ind]
      temp = df.aa_rank[aaref,]
      temp = temp[temp %in% aaalts]
      ind2 = ind[which(df.raw$aaalt[ind] == temp[1])]
    } else {
      ind2 = ind
    }
    ind_selected = c(ind_selected, ind2)
  }
  df.out = rbind(df.out, df.raw[ind_selected,])
  return(df.out)
}

compare_mse_wrapper <- function(df.all, permu_num = 10, permu_size = 100, edit_step = 10, invert_score = F){
  df.func_score.all = get_norm_score(df.all, pos_mean_method = "funsum", pos_norm = F, sd_norm = T, invert_score = invert_score)
  
  ls.mse = list()
  # pb = txtProgressBar(min = 0, max = permu_num, initial = 1, style = 3) 
  for (i in 1:permu_num){
    permu_ind = sample(x = 1:nrow(df.func_score.all), size = permu_size, replace = F)
    df.all_permu = df.all[df.all$aapos %in% df.func_score.all$aapos[permu_ind],]
    min_edit_num = 10
    max_edit_num = nrow(df.all_permu)
    
    for (edit_num in seq(min_edit_num, max_edit_num, by = edit_step)){
      df.rand_permu = sample_rand(df.all_permu, edit_num = edit_num)
      df.rank_permu = sample_rank(df.all_permu, edit_num = edit_num, df.aa_rank)
      # df.rank_permu = sample_rank_v2(df.all_permu, edit_num = edit_num, df.aa_rank_linear)
      
      df.func_score.all_permu = get_norm_score(df.all_permu, pos_mean_method = "funsum", pos_norm = F, sd_norm = T, invert_score = invert_score)
      df.func_score.rand_permu = get_norm_score(df.rand_permu, pos_mean_method = "funsum", pos_norm = F, sd_norm = T, invert_score = invert_score)
      df.func_score.rank_permu = get_norm_score(df.rank_permu, pos_mean_method = "funsum", pos_norm = F, sd_norm = T, invert_score = invert_score)
      
      ind = match(x = df.func_score.rand_permu$aapos, table = df.func_score.all_permu$aapos)
      mse_rand = sse(df.func_score.all_permu$pos_mean[ind], df.func_score.rand_permu$pos_mean)/permu_size
      ind = match(x = df.func_score.rank_permu$aapos, table = df.func_score.all_permu$aapos)
      mse_rank = sse(df.func_score.all_permu$pos_mean[ind], df.func_score.rank_permu$pos_mean)/permu_size
      
      ls.mse[["rand"]][[paste0("edits: ", edit_num)]] = c(ls.mse[["rand"]][[paste0("edits: ", edit_num)]], mse_rand)
      ls.mse[["rank"]][[paste0("edits: ", edit_num)]] = c(ls.mse[["rank"]][[paste0("edits: ", edit_num)]], mse_rank)
    }
    
    # setTxtProgressBar(pb,i)
  }
  # close(pb)
  
  return(ls.mse)
}

get_FUNSUM <- function(df.func_score, avg_method, include_LOF=F){
  aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  }
  df.funsum = matrix(NA, nrow = length(aa_list), ncol = length(aa_list))
  rownames(df.funsum) = aa_list
  colnames(df.funsum) = aa_list
  
  for (aa in aa_list){
    ind = which(df.func_score$aaref == aa)
    
    if (length(ind)>0){
      temp = as.matrix(df.func_score[ind, aa_list])
      
      if (avg_method == "mean"){
        df.funsum[aa,] = colMeans(temp, na.rm = T)
      } else if (avg_method == "median"){
        df.funsum[aa,] = apply(temp, 2, FUN = median, na.rm = T)
      } else if (avg_method == "js"){
        df.funsum[aa,] = get_js(temp)
      } else {
        stop()
      }
    }
  }
  
  return(df.funsum)
}

funsum_to_subTable <- function(df.funsum){
  df.sub_tb = c()
  for (i in 1:nrow(df.funsum)){
    temp = data.frame(aa_pair = paste0(rownames(df.funsum)[i], colnames(df.funsum)), score = as.vector(df.funsum[i,]))
    df.sub_tb = rbind(df.sub_tb, temp)
  }
  
  return(df.sub_tb)
}

add_clinvar <- function(df, gene_id, df.clinvar){
  df$aa_str = paste0(df$aaref, df$aapos, df$aaalt)
  df.clinvar_gene = df.clinvar[which(df.clinvar$genename == gene_id),]
  ind2 = match(x = df$aa_str, table = df.clinvar_gene$aa_pair)
  df$cs_merged = df.clinvar_gene$cs_merged[ind2]
  return(df)
}

plot_func_score <- function(df2, gene_id){
  temp = df2
  temp2 = data.frame(pos = unique(temp$aapos), score_mean = NA)
  for (j in 1:nrow(temp2)){
    ind = which(temp$aapos == temp2$pos[j])
    temp2$score_mean[j] = mean(temp$raw_score[ind], na.rm = T)
  }
  
  ind = order(temp2$pos)
  temp2 = temp2[ind,]
  
  plot(x = temp$aapos, y = temp$raw_score, pch = 16, col = add.alpha("black", 0.25), main = gene_id)
  ind = which(temp$cs_merged == "p")
  points(x = temp$aapos[ind], y = temp$raw_score[ind], pch = 16, col = add.alpha("red", 1))
  ind = which(temp$cs_merged == "b")
  points(x = temp$aapos[ind], y = temp$raw_score[ind], pch = 16, col = add.alpha("blue", 1))
  points(x = temp2$pos, y = temp2$score_mean, type = "l")
}

add_phyloP <- function(df, gene_id, df.phyloP){
  df$aa_str = paste0(df$aaref, df$aapos, df$aaalt)
  df.phyloP_gene = df.phyloP[which(df.phyloP$genename == gene_id),]
  ind2 = match(x = df$aa_str, table = df.phyloP_gene$aa_pair)
  df$phyloP = df.phyloP_gene$phyloP[ind2]
  return(df)
}

plot_ldl_by_variants <-  function(variants_g1, variants_g2, plot_title, df.v2p, df.patient){
  ind = match(x = variants_g1, table = df.v2p$variant)
  patient_g1 = unique(as.numeric(unlist(strsplit(x = df.v2p$patient[ind], split = "\\|"))))
  ind = match(x = variants_g2, table = df.v2p$variant)
  patient_g2 = unique(as.numeric(unlist(strsplit(x = df.v2p$patient[ind], split = "\\|"))))
  
  ind = match(x = patient_g1, table = df.patient$new_id)
  ldl_g1 = df.patient$ldl[ind]
  ldl_g1 = ldl_g1[which(ldl_g1 > 0)]
  ind = match(x = patient_g2, table = df.patient$new_id)
  ldl_g2 = df.patient$ldl[ind]
  ldl_g2 = ldl_g2[which(ldl_g2 > 0)]
  
  ldl_g1.density = density(ldl_g1, na.rm = T, bw = 10)
  ldl_g2.density = density(ldl_g2, na.rm = T, bw = 10)
  
  res = ks.test(x = ldl_g1, y = ldl_g2) # KS-test
  
  ymax = max(c(ldl_g1.density$y, ldl_g2.density$y))
  xmax = max(abs(c(ldl_g1.density$x, ldl_g2.density$x)))
  plot(ldl_g1.density, xlim = c(0, xmax*1.1), ylim = c(0, ymax*1.1), col = "red", 
       main = paste0(plot_title, "\nKS-test p-val = ", res$p.value))
  lines(ldl_g2.density, col = "blue")
  legend("topright", legend = c("Bottom 25%", "Top 25%"), col = c("red", "blue"), lty = 1)
}


v2p_to_p2v <- function(df.v2p){
  temp2 = c()
  ind = which(df.v2p$carriers != "")
  pb = txtProgressBar(min = 0, max = nrow(df.v2p), initial = 0) 
  for (i in ind){
    patient_ids = unlist(strsplit(df.v2p$carriers[i], split = "\\|"))
    temp = data.frame(patient = patient_ids, variant = df.v2p$variant[i])
    temp2 = rbind(temp2, temp)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  patient_ids = unique(temp2$patient)
  df.p2v = data.frame(patient = patient_ids, variant = NA)
  pb = txtProgressBar(min = 0, max = nrow(df.p2v), initial = 0) 
  for (i in 1:nrow(df.p2v)){
    patient_id = df.p2v$patient[i]
    ind = which(temp2$patient == patient_id)
    variant_str = paste(unique(temp2$variant[ind]), collapse = "\\|")
    df.p2v$variant[i] = variant_str
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(df.p2v)
}


read_ssv <- function(fname, AF_thres=NA){
  require(tidyverse)
  df.v2p = read.table(fname, sep = " ", header = F)
  df.p2v = df.v2p %>%
    filter(nchar(df.v2p$V4) == 1 & nchar(df.v2p$V4) == 1) %>%
    mutate(variant = paste(sub(pattern = "chr", replacement = "", V1), V3, V4, V5, sep = "-"),
           patient = V7, AF = V6) %>%
    separate_rows(patient, sep = "\\|") %>%
    filter(patient != "") %>%
    select(patient, variant, AF)
  
  if (is.numeric(AF_thres)){
    print(paste0("Filter variants by AF < ", AF_thres, " ..."))
    df.p2v = df.p2v %>% filter(AF < AF_thres)
  }
  
  return(df.p2v)
}

sse <- function(score1, score2){
  return(sum((score1-score2)^2, na.rm = T))
}

de_noise <- function(df, pos_mean_method, df.funsum, include_LOF = T, show_func_class=F){
  # filter out row without amino acid changes
  ind = which(!is.na(df$aaref) & !is.na(df$aaalt))
  df = df[ind,]
  
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
    ind = which((df$aaref != "*") & (df$aaalt != "*"))
    df = df[ind,]
  }
  df.sub_tb = funsum_to_subTable(df.funsum) #convert FUNSUM to tabular format
  
  if (is.null(df[["gene"]])){
    df[["gene"]] = "gene"
  }
  df$gene_aa_str = paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)
  
  # collapse rows with the same amino acid substitutions
  df = df %>% group_by(gene_aa_str) %>%
    summarise(gene = unique(gene), aapos = unique(aapos), aaref = unique(aaref), aaalt = unique(aaalt), 
              raw_score = mean(raw_score), norm_raw_score = mean(norm_raw_score))
  
  ## calculate positional component
  df.pos_score = df %>% group_by(gene, aapos) %>% summarise(aaref = unique(aaref), pos_mean = NA)
  df.pos_score$gene_aapos = paste0(df.pos_score$gene, "---", df.pos_score$aapos)
  
  temp = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
  colnames(temp) <- aa_list
  if (pos_mean_method == "funsum"){
    temp2 = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
    colnames(temp2) <- aa_list
  }
  
  # pb = txtProgressBar(min = 1, max = nrow(df.pos_score), initial = 1, style = 3) 
  for (i in 1:nrow(df.pos_score)){
    ind = which(df$gene == df.pos_score$gene[i] & df$aapos == df.pos_score$aapos[i])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$norm_raw_score[ind[ind2]]
    
    # calculate pos_mean by funsum method
    if (pos_mean_method == "funsum"){
      ind = which(!is.na(temp[i,]))
      temp2[i,ind] = temp[i,ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
    }
    
    # setTxtProgressBar(pb,i)
  }
  # close(pb)
  
  # calculate pos_mean by other methods
  if (pos_mean_method == "mean"){
    df.pos_score$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.pos_score$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.pos_score$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "funsum"){
    df.pos_score$pos_mean = get_js(t(temp2))
  }
  
  ## construct a new df with all possible substitutions
  df.out = df.pos_score %>% dplyr::select(gene, aapos, aaref) %>% 
    dplyr::slice(rep(1:n(), each = length(aa_list)))
  df.out$aaalt = rep(aa_list, nrow(df.pos_score))
  
  # assign functional class
  if (show_func_class){
    df.out$functional_class = NA
    ind = which(df.out$aaref != df.out$aaalt & df.out$aaalt != '*')
    df.out$functional_class[ind] = "MIS"
    ind = which(df.out$aaref != "*" & df.out$aaalt == "*")
    df.out$functional_class[ind] = "LOF"
    ind = which(df.out$aaref == df.out$aaalt)
    df.out$functional_class[ind] = "SYN"
  }
  
  # assign norm_score, pos_score, sub_score to df.out
  df.out$gene_aa_str = paste0(df.out$gene, "---", df.out$aaref, df.out$aapos, df.out$aaalt)
  df.out$gene_aapos = paste0(df.out$gene, "---", df.out$aapos)
  df.out$aa_pair = paste0(df.out$aaref, df.out$aaalt)
  ind = match(df.out$gene_aa_str, table = df$gene_aa_str)
  df.out$raw_score = df$raw_score[ind]
  df.out$norm_raw_score = df$norm_raw_score[ind]
  ind = match(df.out$gene_aapos, table = df.pos_score$gene_aapos)
  df.out$pos_score = df.pos_score$pos_mean[ind]
  ind = match(df.out$aa_pair, table = df.sub_tb$aa_pair)
  df.out$sub_score = df.sub_tb$score[ind]
  df.out$final_score = df.out$pos_score + df.out$sub_score
  df.out$final_score_lite = df.out$final_score
  ind = which(is.na(df.out$norm_raw_score))
  df.out$final_score_lite[ind] = NA
  
  return(df.out)
}

add_blosum_funsum <- function(df.NSFP, df.blosum, df.funsum){
  df.sub_tb = funsum_to_subTable(df.blosum) #convert BLOSUM to tabular format
  df.sub_tb2 = funsum_to_subTable(df.funsum) #convert FUNSUM to tabular format
  ind = match(x = df.NSFP$aa_pair, table = df.sub_tb$aa_pair)
  df.NSFP$blosum = df.sub_tb$score[ind]
  df.NSFP$funsum = df.sub_tb2$score[ind]
  return(df.NSFP)
}

de_noise_ss_1gene <- function(df, pos_mean_method, df.funsum, ls.funsum_ss, dss_path, include_LOF = T, show_func_class=F){
  # filter out row without amino acid changes
  ind = which(!is.na(df$aaref) & !is.na(df$aaalt))
  df = df[ind,]
  
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
    ind = which((df$aaref != "*") & (df$aaalt != "*"))
    df = df[ind,]
  }
  df.sub_tb = funsum_to_subTable(df.funsum) #convert FUNSUM to tabular format
  if (is.null(df[["gene"]])){
    df[["gene"]] = "gene"
  }
  df$gene_aa_str = paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)
  
  # collapse rows with the same amino acid substitutions
  df = df %>% group_by(gene_aa_str) %>%
    summarise(gene = unique(gene), aapos = unique(aapos), aaref = unique(aaref), aaalt = unique(aaalt), 
              raw_score = mean(raw_score), norm_raw_score = mean(norm_raw_score))
  
  ## calculate positional component
  df.pos_score = df %>% group_by(gene, aapos) %>% summarise(aaref = unique(aaref), pos_mean = NA)
  df.pos_score$gene_aapos = paste0(df.pos_score$gene, "---", df.pos_score$aapos)
  
  temp = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
  colnames(temp) <- aa_list
  if (pos_mean_method == "funsum"){
    temp2 = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
    colnames(temp2) <- aa_list
  }
  
  pb = txtProgressBar(min = 1, max = nrow(df.pos_score), initial = 1, style = 3) 
  for (i in 1:nrow(df.pos_score)){
    ind = which(df$gene == df.pos_score$gene[i] & df$aapos == df.pos_score$aapos[i])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$norm_raw_score[ind[ind2]]
    
    # calculate pos_mean by funsum method
    if (pos_mean_method == "funsum"){
      ind = which(!is.na(temp[i,]))
      temp2[i,ind] = temp[i,ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
    }
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # calculate pos_mean by other methods
  if (pos_mean_method == "mean"){
    df.pos_score$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.pos_score$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.pos_score$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "funsum"){
    df.pos_score$pos_mean = get_js(t(temp2))
  }
  
  ## construct a new df with all possible substitutions
  df.out = df.pos_score %>% select(gene, aapos, aaref) %>% 
    dplyr::slice(rep(1:n(), each = length(aa_list)))
  df.out$aaalt = rep(aa_list, nrow(df.pos_score))
  
  # assign functional class
  if (show_func_class){
    df.out$functional_class = NA
    ind = which(df.out$aaref != df.out$aaalt & df.out$aaalt != '*')
    df.out$functional_class[ind] = "MIS"
    ind = which(df.out$aaref != "*" & df.out$aaalt == "*")
    df.out$functional_class[ind] = "LOF"
    ind = which(df.out$aaref == df.out$aaalt)
    df.out$functional_class[ind] = "SYN"
  }
  
  # assign norm_score, pos_score, sub_score to df.out
  df.out$gene_aa_str = paste0(df.out$gene, "---", df.out$aaref, df.out$aapos, df.out$aaalt)
  df.out$gene_aapos = paste0(df.out$gene, "---", df.out$aapos)
  df.out$aa_pair = paste0(df.out$aaref, df.out$aaalt)
  
  # add dssp SS annotation
  require(ptm)
  df.out$ss = NA
  df.out$acc = NA
  df.dssp = parse.dssp(file = dss_path, keepfiles = T)
  ind = match(df.out$aapos, table = df.dssp$respdb)
  df.out$ss = df.dssp$ss[ind]
  df.out$acc = df.dssp$sasa[ind]
  
  ind = match(df.out$gene_aa_str, table = df$gene_aa_str)
  df.out$raw_score = df$raw_score[ind]
  df.out$norm_raw_score = df$norm_raw_score[ind]
  ind = match(df.out$gene_aapos, table = df.pos_score$gene_aapos)
  df.out$pos_score = df.pos_score$pos_mean[ind]
  ind = match(df.out$aa_pair, table = df.sub_tb$aa_pair)
  df.out$sub_score = df.sub_tb$score[ind]
  
  ss_list = c("G" = "Helices", "H" = "Helices", "I" = "Helices", 
              "E" = "Strands", "B" = "Strands", 
              "T" = "Loops", "S" = "Loops", "C" = "Loops")
  df.out$sub_score_ss = NA
  for (ss in names(ss_list)){
    ind = which(df.out$ss == ss)
    if (length(ind)){
      df.sub_tb = funsum_to_subTable(ls.funsum_ss[[ss]]) #convert FUNSUM to tabular format
      ind2 = match(df.out$aa_pair[ind], table = df.sub_tb$aa_pair)
      df.out$sub_score_ss[ind] = df.sub_tb$score[ind2]
    }
  }
  
  df.out$final_score = df.out$pos_score + df.out$sub_score
  df.out$final_score_ss = df.out$pos_score + df.out$sub_score_ss
  ind = which(is.na(df.out$final_score_ss))
  df.out$final_score_ss[ind] = df.out$final_score[ind]
  df.out$sub_score_ss[ind] = df.out$sub_score[ind]
  
  return(df.out)
}

de_noise_ss <- function(df, pos_mean_method, df.funsum, ls.funsum_ss, include_LOF = T, show_func_class=F){
  # filter out row without amino acid changes
  ind = which(!is.na(df$aaref) & !is.na(df$aaalt))
  df = df[ind,]
  
  if (include_LOF){
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  } else {
    aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW", split = ""))
    ind = which((df$aaref != "*") & (df$aaalt != "*"))
    df = df[ind,]
  }
  df.sub_tb = funsum_to_subTable(df.funsum) #convert FUNSUM to tabular format
  if (is.null(df[["gene"]])){
    df[["gene"]] = "gene"
  }
  df$gene_aa_str = paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)
  
  # collapse rows with the same amino acid substitutions
  df = df %>% group_by(gene_aa_str) %>%
    summarise(gene = unique(gene), aapos = unique(aapos), aaref = unique(aaref), aaalt = unique(aaalt), 
              raw_score = mean(raw_score), norm_raw_score = mean(norm_raw_score))
  
  ## calculate positional component
  df.pos_score = df %>% group_by(gene, aapos) %>% summarise(aaref = unique(aaref), pos_mean = NA)
  df.pos_score$gene_aapos = paste0(df.pos_score$gene, "---", df.pos_score$aapos)
  
  temp = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
  colnames(temp) <- aa_list
  if (pos_mean_method == "funsum"){
    temp2 = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
    colnames(temp2) <- aa_list
  }
  
  pb = txtProgressBar(min = 1, max = nrow(df.pos_score), initial = 1, style = 3) 
  for (i in 1:nrow(df.pos_score)){
    ind = which(df$gene == df.pos_score$gene[i] & df$aapos == df.pos_score$aapos[i])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$norm_raw_score[ind[ind2]]
    
    # calculate pos_mean by funsum method
    if (pos_mean_method == "funsum"){
      ind = which(!is.na(temp[i,]))
      temp2[i,ind] = temp[i,ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
    }
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # calculate pos_mean by other methods
  if (pos_mean_method == "mean"){
    df.pos_score$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.pos_score$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.pos_score$pos_mean = get_js(t(temp))
  } else if (pos_mean_method == "funsum"){
    df.pos_score$pos_mean = get_js(t(temp2))
  }
  
  ## construct a new df with all possible substitutions
  df.out = df.pos_score %>% select(gene, aapos, aaref) %>% 
    dplyr::slice(rep(1:n(), each = length(aa_list)))
  df.out$aaalt = rep(aa_list, nrow(df.pos_score))
  
  # assign functional class
  if (show_func_class){
    df.out$functional_class = NA
    ind = which(df.out$aaref != df.out$aaalt & df.out$aaalt != '*')
    df.out$functional_class[ind] = "MIS"
    ind = which(df.out$aaref != "*" & df.out$aaalt == "*")
    df.out$functional_class[ind] = "LOF"
    ind = which(df.out$aaref == df.out$aaalt)
    df.out$functional_class[ind] = "SYN"
  }
  
  # assign norm_score, pos_score, sub_score to df.out
  df.out$gene_aa_str = paste0(df.out$gene, "---", df.out$aaref, df.out$aapos, df.out$aaalt)
  df.out$gene_aapos = paste0(df.out$gene, "---", df.out$aapos)
  df.out$aa_pair = paste0(df.out$aaref, df.out$aaalt)
  
  # add dssp SS annotation
  require(ptm)
  df.out$ss = NA
  df.out$acc = NA
  for (gene in unique(df.out$gene)){
    df.dssp = parse.dssp(file = paste0("./alphaFold/dssp_out/", gene, ".dss"), keepfiles = T)
    df.dssp$gene_aapos = paste0(gene, "---", df.dssp$respdb)
    ind = which(df.out$gene == gene)
    ind2 = match(df.out$gene_aapos[ind], table = df.dssp$gene_aapos)
    df.out$ss[ind] = df.dssp$ss[ind2]
    df.out$acc[ind] = df.dssp$sasa[ind2]
  }
  
  ind = match(df.out$gene_aa_str, table = df$gene_aa_str)
  df.out$raw_score = df$raw_score[ind]
  df.out$norm_raw_score = df$norm_raw_score[ind]
  ind = match(df.out$gene_aapos, table = df.pos_score$gene_aapos)
  df.out$pos_score = df.pos_score$pos_mean[ind]
  ind = match(df.out$aa_pair, table = df.sub_tb$aa_pair)
  df.out$sub_score = df.sub_tb$score[ind]
  
  ss_list = c("G" = "Helices", "H" = "Helices", "I" = "Helices", 
              "E" = "Strands", "B" = "Strands", 
              "T" = "Loops", "S" = "Loops", "C" = "Loops")
  df.out$sub_score_ss = NA
  for (ss in names(ss_list)){
    ind = which(df.out$ss == ss)
    if (length(ind)){
      df.sub_tb = funsum_to_subTable(ls.funsum_ss[[ss]]) #convert FUNSUM to tabular format
      ind2 = match(df.out$aa_pair[ind], table = df.sub_tb$aa_pair)
      df.out$sub_score_ss[ind] = df.sub_tb$score[ind2]
    }
  }
  
  df.out$final_score = df.out$pos_score + df.out$sub_score
  df.out$final_score_ss = df.out$pos_score + df.out$sub_score_ss
  ind = which(is.na(df.out$final_score_ss))
  df.out$final_score_ss[ind] = df.out$final_score[ind]
  df.out$sub_score_ss[ind] = df.out$sub_score[ind]
  
  return(df.out)
}

# normalization by sd
check_direction_2 <- function(df, normalize=T){
  msg = c()
  gene_id = unique(df$gene)
  
  ## normalization with synonymous and LOF substitutions
  ind = which(df$functional_class == "SYN")
  ind2 = which(df$functional_class == "LOF")
  if (length(ind)>0 & length(ind2)>0){
    m1 = median(df$raw_score[ind])
    m2 = median(df$raw_score[ind2])
    
    if (m2 < m1){
      df$raw_score = -df$raw_score
      msg = c(msg, paste0("Info: ", gene_id, " DMS data in opposite direction. Score inverted."))
    }
    m1 = median(df$raw_score[ind])
    m2 = median(df$raw_score[ind2])
    
    if (normalize){
      df$norm_raw_score = (df$raw_score - median(df$raw_score, na.rm = T))/sd(df$raw_score, na.rm = T)
    }
  } else {
    if (!length(ind)){
      msg = c(msg, paste0("Error: ", gene_id, " DMS data contains no synonymous variants!"))
    }
    if (!length(ind2)){
      msg = c(msg, paste0("Error: ", gene_id, " DMS data contains no LOF variants!"))
    }
    msg = c(msg, paste0("Warning: ", gene_id, " DMS data direction cannot be determined! Gene skipped."))
  }
  
  return(list(df, msg))
}

check_direction <- function(df, normalize=T){
  msg = c()
  gene_id = unique(df$gene)
  
  ## normalization with synonymous and LOF substitutions
  ind = which(df$functional_class == "SYN")
  ind2 = which(df$functional_class == "LOF")
  if (length(ind)>0 & length(ind2)>0){
    m1 = median(df$raw_score[ind])
    m2 = median(df$raw_score[ind2])
    
    if (m2 < m1){
      df$raw_score = -df$raw_score
      msg = c(msg, paste0("Info: ", gene_id, " DMS data in opposite direction. Score inverted."))
    }
    m1 = median(df$raw_score[ind])
    m2 = median(df$raw_score[ind2])
    
    if (normalize){
      df$norm_raw_score = (df$raw_score - m1)/(m2-m1)
    }
  } else {
    if (!length(ind)){
      msg = c(msg, paste0("Error: ", gene_id, " DMS data contains no synonymous variants!"))
    }
    if (!length(ind2)){
      msg = c(msg, paste0("Error: ", gene_id, " DMS data contains no LOF variants!"))
    }
    msg = c(msg, paste0("Warning: ", gene_id, " DMS data direction cannot be determined! Gene skipped."))
  }
  
  return(list(df, msg))
}

draw_variant_plot <- function(df.NSFP, score_types, title_str="", ks_text_xpos=0.8, gold_star_thres = 0, mw_test = F){
  df.NSFP$gold_star[is.na(df.NSFP$gold_star)] = 0
  df.NSFP = df.NSFP[df.NSFP$gold_star >= gold_star_thres,]
  
  for (i in 1:length(score_types)){
    score_type = score_types[i]
    has_g1 = T
    has_g2 = T
    
    ind = which(df.NSFP$cs_simple == "b")
    score_g1 = df.NSFP[[score_type]][ind]
    score_g1 = score_g1[which(!is.na(score_g1))]
    if (length(score_g1) > 0){
      score_g1.density = density(score_g1, bw = 0.5)
    } else {
      has_g1 = F
    }
    
    ind = which(df.NSFP$cs_simple == "p")
    score_g2 = df.NSFP[[score_type]][ind]
    score_g2 = score_g2[which(!is.na(score_g2))]
    if (length(score_g2) > 0){
      score_g2.density = density(score_g2, bw = 0.5)
    } else {
      has_g2 = F
    }
    
    if (has_g1 & has_g2){
      p1 = NA
      res = ks.test(x = score_g1, y = score_g2, alternative = "greater") # KS-test
      p1 = format(res$p.value, digits = 3)
      
      if (mw_test){
        p2 = NA
        res = wilcox.test(x = score_g1, y = score_g2, alternative = "less") # KS-test
        p2 = format(res$p.value, digits = 3)
      }
      
      ymax = max(c(score_g1.density$y, score_g2.density$y))
      xmax = max(abs(c(score_g1.density$x, score_g2.density$x)))
      
      plot(score_g1.density, col = "blue", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g1 is blue
      text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*ks_text_xpos, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.9, paste0("KS p=", p1))
      if (mw_test){
        text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*ks_text_xpos, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.8, paste0("MW p-val = ", p2))
      }
      
      lines(score_g2.density, col = "red") # g2 is red
      if (i == 1){
        mtext(title_str, side = 3, outer = T, line = -1.5)
      }
    } else if (has_g1) {
      plot(score_g1.density, col = "blue", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g1 is blue
    } else if (has_g2) {
      plot(score_g2.density, col = "red", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g2 is red
    }
    
    legend("topleft", legend = c(paste0("B/LB (n=", length(score_g1), ")"), paste0("P/LP (n=", length(score_g2), ")")), col = c("blue", "red"), lty = 1)
  }
}

draw_variant_ggplot <- function(df.NSFP, score_type, show_legend=T, show_stat=T){
  require(hrbrthemes)
  if (score_type %in% colnames(df.NSFP)){
    df.NSFP2 = df.NSFP %>% filter(!is.na(cs_simple) & !is.na(df.NSFP[[score_type]]))
    
    score_g1 = df.NSFP2[[score_type]][df.NSFP2$cs_simple == "b"]
    score_g2 = df.NSFP2[[score_type]][df.NSFP2$cs_simple == "p"]
    
    if (length(score_g1)>0 & length(score_g2)>0){
      res = ks.test(x = score_g1, y = score_g2, alternative = "greater") # KS-test
      p1 = format(res$p.value, digits = 3)
      res = wilcox.test(x = score_g1, y = score_g2, alternative = "less") # KS-test
      p2 = format(res$p.value, digits = 3)
    }
    
    plt <- ggplot(data = df.NSFP2, aes(x=.data[[score_type]], group=cs_simple, fill=cs_simple)) + 
      geom_density(alpha=0.5, show.legend = show_legend) + 
      scale_x_continuous(name="Functional score", limits=c(-3.5,3.5)) +
      scale_y_continuous(name="Density", limits=c(0,1)) + 
      scale_fill_manual(values = c("blue", "red"), name = "Clinvar", labels=c('Benign/Likely Benign', 'Pathogenic/Likely Pathogenic')) + 
      theme_ipsum()
    
    if (show_stat){
      plt <- plt + annotate("text", x=0, y=1, label= paste0("b=", length(score_g1), ", p=", length(score_g2), ", KS p=", p1, ", MW p=", p2))
    }
    
    return(plt)
  } else {
    simpleError("score_type not found in df.NSFP!")
  }
}

draw_variant_ggplot_v2 <- function(df.NSFP, score_type, score_type_name=score_type, ymax = 1, bw = 0.2, x_range = c(-5, 5), show_legend = F, legend_pos = c(0.3, 0.9), test="KS"){
  require(hrbrthemes)
  ind = which(df.NSFP$cs_simple == "b")
  ind2 = which(df.NSFP$cs_simple == "p")
  
  pval_str = "No test"
  if (test == "KS"){
    pval = get_KS_pval(x = df.NSFP[[score_type]][ind], y = df.NSFP[[score_type]][ind2])
    pval_str = paste0("KS p=", format(as.numeric(as.character(pval)), digits=2))
    print(pval_str)
  } else if (test == "MW"){
    pval = get_MW_pval(x = df.NSFP[[score_type]][ind], y = df.NSFP[[score_type]][ind2])
    pval_str = paste0("MW p=", format(as.numeric(as.character(pval)), digits=2))
    print(pval_str)
  }
  
  n_g1 = length(which(!is.na(df.NSFP[[score_type]][ind])))
  n_g2 = length(which(!is.na(df.NSFP[[score_type]][ind2])))
  plt <- df.NSFP %>% filter(!is.na(cs_simple) & !is.na(get(score_type))) %>% 
    ggplot(aes(x=get(score_type), group=cs_simple, fill=cs_simple)) + geom_density(alpha=0.5, show.legend = show_legend, bw = bw) + 
    scale_x_continuous(name="Functional score", limits=x_range) + scale_y_continuous(name="Density", limits=c(0,ymax)) + 
    scale_fill_manual(values = c("blue", "red"), name = NULL, labels=c('B/LB', 'P/LP')) + 
    geom_label(label=pval_str, x=max(x_range)-abs(x_range[2]-x_range[1])*0.02, y=ymax*0.98, color = "black", fill="white", hjust=1, vjust=1, label.size=NA) + 
    geom_label(label=paste0("N=", n_g1), x=min(x_range+abs(x_range[2]-x_range[1])*0.02), y=ymax*0.02, color = "blue", fill="white", hjust=0, vjust=0, label.size=NA) + 
    geom_label(label=paste0("N=", n_g2), x=max(x_range-abs(x_range[2]-x_range[1])*0.02), y=ymax*0.02, color = "red", fill="white", hjust=1, vjust=0, label.size=NA) +
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c", plot_margin = margin(10, 10, 10, 10)) + 
    ggtitle(score_type_name) + theme(plot.title = element_text(hjust = 0.5), legend.position = legend_pos, legend.background = element_rect(fill = "white", color = "white"), legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  return(plt)
}

draw_patient_ggplot <- function(df.p2v, score_type, score_type_name=score_type, ymax = 1, bw = 0.2, x_range = c(-5, 5), show_legend = F, legend_pos = c(0.3, 0.9), test="KS"){
  require(hrbrthemes)
  ind = which(df.p2v$cancer == 0)
  ind2 = which(df.p2v$cancer == 1)
  
  pval_str = "No test"
  if (test == "KS"){
    pval = get_KS_pval(x = df.p2v[[score_type]][ind], y = df.p2v[[score_type]][ind2])
    pval_str = paste0("KS p=", format(as.numeric(as.character(pval)), digits=2))
    print(pval_str)
  } else if (test == "MW"){
    pval = get_MW_pval(x = df.p2v[[score_type]][ind], y = df.p2v[[score_type]][ind2])
    pval_str = paste0("MW p=", format(as.numeric(as.character(pval)), digits=2))
    print(pval_str)
  }
  
  n_g1 = length(which(!is.na(df.p2v[[score_type]][ind])))
  n_g2 = length(which(!is.na(df.p2v[[score_type]][ind2])))
  plt <- df.p2v %>% filter(!is.na(cancer) & !is.na(get(score_type))) %>% 
    ggplot(aes(x=get(score_type), group=cancer, fill=cancer)) + geom_density(alpha=0.5, show.legend = show_legend, bw = bw) + 
    scale_x_continuous(name="Functional score", limits=x_range) + scale_y_continuous(name="Density", limits=c(0,ymax)) + 
    scale_fill_manual(values = c("blue", "red"), name = NULL, labels=c('No cancer', 'Cancer')) + 
    geom_label(label=pval_str, x=max(x_range)-abs(x_range[2]-x_range[1])*0.02, y=ymax*0.98, color = "black", fill="white", hjust=1, vjust=1, label.size=NA) + 
    geom_label(label=paste0("N=", n_g1), x=min(x_range+abs(x_range[2]-x_range[1])*0.02), y=ymax*0.02, color = "blue", fill="white", hjust=0, vjust=0, label.size=NA) + 
    geom_label(label=paste0("N=", n_g2), x=max(x_range-abs(x_range[2]-x_range[1])*0.02), y=ymax*0.02, color = "red", fill="white", hjust=1, vjust=0, label.size=NA) +
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c", plot_margin = margin(10, 10, 10, 10)) + 
    ggtitle(score_type_name) + theme(plot.title = element_text(hjust = 0.5), legend.position = legend_pos, legend.background = element_rect(fill = "white", color = "white"), legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  return(plt)
}


draw_patient_ggplot_v2 <- function(df.p2v, score_type, score_type_name=score_type, phenotype_col = "cancer", legend_labels = c('No cancer', 'Cancer'), ymax = 1, bw = 0.2, x_range = c(-5, 5), show_legend = F, legend_pos = c(0.3, 0.9), test="KS"){
  require(hrbrthemes)
  ind = which(df.p2v[[phenotype_col]] == 0)
  ind2 = which(df.p2v[[phenotype_col]] == 1)
  
  pval_str = "No test"
  if (test == "KS"){
    pval = get_KS_pval(x = df.p2v[[score_type]][ind], y = df.p2v[[score_type]][ind2])
    pval_str = paste0("KS p=", format(as.numeric(as.character(pval)), digits=2))
    print(pval_str)
  } else if (test == "MW"){
    pval = get_MW_pval(x = df.p2v[[score_type]][ind], y = df.p2v[[score_type]][ind2])
    pval_str = paste0("MW p=", format(as.numeric(as.character(pval)), digits=2))
    print(pval_str)
  }
  
  n_g1 = length(which(!is.na(df.p2v[[score_type]][ind])))
  n_g2 = length(which(!is.na(df.p2v[[score_type]][ind2])))
  plt <- df.p2v %>% filter(!is.na(get(phenotype_col)) & !is.na(get(score_type))) %>% 
    ggplot(aes(x=get(score_type), group=get(phenotype_col), fill=get(phenotype_col))) + geom_density(alpha=0.5, show.legend = show_legend, bw = bw) + 
    scale_x_continuous(name="Functional score", limits=x_range) + scale_y_continuous(name="Density", limits=c(0,ymax)) + 
    scale_fill_manual(values = c("blue", "red"), name = NULL, labels=legend_labels) + 
    geom_label(label=pval_str, x=max(x_range)-abs(x_range[2]-x_range[1])*0.02, y=ymax*0.99, color = "black", fill="white", hjust=1, vjust=1, label.size=NA) + 
    geom_label(label=paste0("N=", n_g1), x=min(x_range+abs(x_range[2]-x_range[1])*0.02), y=ymax*0.02, color = "blue", fill="white", hjust=0, vjust=0, label.size=NA) + 
    geom_label(label=paste0("N=", n_g2), x=max(x_range-abs(x_range[2]-x_range[1])*0.02), y=ymax*0.02, color = "red", fill="white", hjust=1, vjust=0, label.size=NA) +
    theme_ipsum(base_family = "Arial", plot_title_size = 14, axis_title_size = 12, axis_title_face = "bold", axis_title_just = "c", plot_margin = margin(10, 10, 10, 10)) + 
    ggtitle(score_type_name) + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = legend_pos, 
          legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
          legend.background = element_rect(fill = "white", color = "white"), 
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  return(plt)
}


draw_variant_plot_cancerGeneBE <- function(df, score_types, g1_keyword, g2_keyword, title_str="", ks_text_xpos=0.8){
  for (i in 1:length(score_types)){
    score_type = score_types[i]
    has_g1 = T
    has_g2 = T
    
    ind = which(df$class_simple == g1_keyword)
    score_g1 = df[[score_type]][ind]
    score_g1 = score_g1[which(!is.na(score_g1))]
    if (length(score_g1) > 0){
      score_g1.density = density(score_g1, bw = 0.5)
    } else {
      has_g1 = F
    }
    
    ind = which(df$class_simple == g2_keyword)
    score_g2 = df[[score_type]][ind]
    score_g2 = score_g2[which(!is.na(score_g2))]
    if (length(score_g2) > 0){
      score_g2.density = density(score_g2, bw = 0.5)
    } else {
      has_g2 = F
    }
    
    if (has_g1 & has_g2){
      p1 = NA
      res = ks.test(x = score_g1, y = score_g2, alternative = "greater") # KS-test
      p1 = format(res$p.value, digits = 3)
      
      ymax = max(c(score_g1.density$y, score_g2.density$y))
      xmax = max(abs(c(score_g1.density$x, score_g2.density$x)))
      
      plot(score_g1.density, col = "blue", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g1 is blue
      text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*ks_text_xpos, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.9, paste0("KS-test p-val = ", p1))
      lines(score_g2.density, col = "red") # g2 is red
      if (i == 1){
        mtext(title_str, side = 3, outer = T, line = -1.5)
      }
    } else if (has_g1) {
      plot(score_g1.density, col = "blue", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g1 is blue
    } else if (has_g2) {
      plot(score_g2.density, col = "red", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g2 is red
    }
    
    legend("topleft", legend = c(paste0(g1_keyword, " (n=", length(score_g1), ")"), paste0(g2_keyword, " (n=", length(score_g2), ")")), col = c("blue", "red"), lty = 1)
  }
}

draw_patient_plot_BRCA1 <- function(df.NSFP, score_types, title_str="", age_thres_g1=0, age_thres_g2=100, prs_filter=10, fh_breastcancer=c(0,1), ks_text_xpos=0.8, y_labels = NA, draw_plot=T){
  ymax = 1.5
  x_range = c(-4,4)
  plot_method = "density"
  bw = 0.25
  kernel = "gaussian"
  # load patient phenotypes, group patients by breastcancer=0 vs breastcancer=1
  df.patient = read.csv("./UKBB/updated_patient_to_phenotype_041721.csv") # load patient info
  # df.patient = df.patient[which(df.patient$breastcancer_age < age_thres & df.patient$normalized_breastcancer_PRS < prs_filter & (df.patient$fh_breastcancer %in% fh_breastcancer)),]
  df.patient = df.patient[which(df.patient$normalized_breastcancer_PRS < prs_filter & (df.patient$fh_breastcancer %in% fh_breastcancer)),]
  ind = which(df.patient$breastcancer == 0 & df.patient$sex == "F" & df.patient$breastcancer_age > age_thres_g1)
  patient_g1 = df.patient$new_id[ind] # g1 is non-cancer group
  ind = which(df.patient$breastcancer == 1 & df.patient$sex == "F" & df.patient$breastcancer_age < age_thres_g2)
  patient_g2 = df.patient$new_id[ind] # g2 is cancer group
  
  df.p2v = read.csv("./UKBB/BRCA1_missenseLOF_p2v.csv") # load UKBB BRCA1 variants
  df.p2v$patient_label = NA
  df.p2v$patient_label[df.p2v$patient %in% patient_g1] = "g1"
  df.p2v$patient_label[df.p2v$patient %in% patient_g2] = "g2"
  
  df.p2v[,score_types] = NA
  df.p2v$cs_simple = NA
  ind = grep(pattern = "\\|", x = df.p2v$variant)
  for (i in ind){
    variant_str = unlist(strsplit(x = df.p2v$variant[i], split = "\\|"))
    ind2 = match(variant_str, table = df.NSFP$variant_str)
    df.p2v[i,score_types] = colMeans(df.NSFP[ind2,score_types], na.rm = T)
    # df.p2v[i,score_types] = apply(df.NSFP[ind2,score_types], 2, median, na.rm = T)
    temp = unique(df.NSFP$cs_simple[ind2])
    temp = temp[!is.na(temp)]
    if (length(temp) == 1){
      df.p2v$cs_simple[i] = paste0(temp, collapse = "")
    }
  }
  
  ind2 = match(df.p2v$variant[-ind], table = df.NSFP$variant_str)
  df.p2v[-ind,score_types] = df.NSFP[ind2,score_types]
  df.p2v$cs_simple[-ind] = df.NSFP$cs_simple[ind2]
  
  # plot score distribution for breatcancer=0 vs breatcancer=1
  ind = which(df.p2v$patient_label == "g1")
  ind2 = which(df.p2v$patient_label == "g2")
  
  # # create data frame to return each score types p-val
  # df.res = c(age_thres_g1, age_thres_g2)
  # names(df.res) = c("age_thres_g1", "age_thres_g2")
  
  df.obs = c()
  # draw plot for each score types
  for (i in 1:length(score_types)){
    score_type = score_types[i]
    score_g1 = df.p2v[[score_type]][ind]
    score_g1 = score_g1[which(!is.na(score_g1))]
    score_g1.density = density(score_g1, bw = bw, kernel = kernel)
    score_g2 = df.p2v[[score_type]][ind2]
    score_g2 = score_g2[which(!is.na(score_g2))]
    score_g2.density = density(score_g2, bw = bw, kernel = kernel)
    
    df.obs = c(df.obs, length(score_g1)+length(score_g2))
    
    p1 = NA
    res = ks.test(x = score_g1, y = score_g2, alternative = "greater") # KS-test
    p1 = format(res$p.value, digits = 3)
    
    # # data frame to return each score types p-val
    # score_name = names(score_types)[i]
    # temp = c(length(score_g1), length(score_g2), as.numeric(p1))
    # names(temp) = c(paste0(score_name, "_g1"), paste0(score_name, "_g2"), paste0(score_name, "_pval"))
    # df.res = c(df.res, temp)
    # # ymax = max(c(score_g1.density$y, score_g2.density$y))
    # # xmax = max(abs(c(score_g1.density$x, score_g2.density$x)))
    
    if (is.na(y_labels)){
      y_label = names(score_types)[i]
    } else if (y_labels == F){
      y_label = ""
    } else {
      y_label = y_labels[i]
    }
    
    if (draw_plot){
      if (plot_method == "density"){
        plot(score_g1.density, col = "blue", ylim=c(0,ymax), xlim=x_range, las=1, main = "", xaxs="i", yaxs="i", ylab=y_label, xlab="") # g1 is blue
        lines(score_g2.density, col = "red") # g2 is black
        legend("topleft", legend = c(paste0("Non-bc (n=", length(score_g1), ")"), paste0("Bc (n=", length(score_g2), ")")), col = c("blue", "red"), lty = 1)
        if (i == 1){
          mtext(title_str, side = 3, outer = T, line = -1.5)
        }
      }
      
      if (plot_method == "histogram"){
        hist(score_g1, breaks = breaks, freq = F, col = add.alpha("blue", alpha = 0.5), ylim=c(0,ymax), xlim=x_range, las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="")
        hist(score_g2, breaks = breaks, freq = F, col = add.alpha("red", alpha = 0.5), ylim=c(0,ymax), xlim=x_range, add=T, las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="")
        legend("topleft", legend = c(paste0("Non-bc (n=", length(score_g1), ")"), paste0("Bc (n=", length(score_g2), ")")), fill = add.alpha(c("blue", "red"), 0.5))
        box()
      }
      
      text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*ks_text_xpos, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.9, paste0("KS-test p-val = ", p1))
    }
  }
  
  names(df.obs) = score_types
  return(df.obs)
}

draw_patient_plot_LDLR <- function(df.NSFP, score_types, title_str="", ldl_col="ldl", ldl_thres_g1=80, ldl_thres_g2=80, age_thres=0, ks_text_xpos=0.8, y_labels = NA, draw_plot=T){
  ymax = 1.5
  x_range = c(-4,4)
  plot_method = "density"
  bw = 0.25
  kernel = "gaussian"
  # load patient phenotypes, group patients with >75% or <75% LDL levels
  df.patient = read.csv("./UKBB/updated_patient_to_phenotype_042822.csv") # load patient info
  # decide ehcih set of ldl values to use
  if (ldl_col == "ldl"){
    df.patient$ldl_used = df.patient$ldl
  } else if (ldl_col == "ldl_James") {
    df.patient$ldl_used = df.patient$ldl_James
  }
  ind = which(!is.na(df.patient$ldl_used) & df.patient$ldl_used > 0)
  df.patient = df.patient[ind,]
  ind = which(df.patient$ldl < quantile(df.patient$ldl, ldl_thres_g1/100) & df.patient$death_age > age_thres)
  patient_g1 = df.patient$new_id[ind] # g1 is ldl less than 75%
  ind = which(df.patient$ldl > quantile(df.patient$ldl,ldl_thres_g2/100) & df.patient$death_age > age_thres)
  patient_g2 = df.patient$new_id[ind] # g2 is ldl more than 75%
  
  df.p2v = read.csv("./UKBB/LDLR_missenseLOF_p2v.csv") # load UKBB LDLR variants
  df.p2v$patient_label = NA
  df.p2v$patient_label[df.p2v$patient %in% patient_g1] = "g1"
  df.p2v$patient_label[df.p2v$patient %in% patient_g2] = "g2"
  
  df.p2v[,score_types] = NA
  df.p2v$cs_simple = NA
  ind = grep(pattern = "\\|", x = df.p2v$variant)
  for (i in ind){
    variant_str = unlist(strsplit(x = df.p2v$variant[i], split = "\\|"))
    ind2 = match(variant_str, table = df.NSFP$variant_str)
    df.p2v[i,score_types] = colMeans(df.NSFP[ind2,score_types], na.rm = T)
    # df.p2v[i,score_types] = apply(df.NSFP[ind2,score_types], 2, median, na.rm = T)
    temp = unique(df.NSFP$cs_simple[ind2])
    temp = temp[!is.na(temp)]
    if (length(temp) == 1){
      df.p2v$cs_simple[i] = paste0(temp, collapse = "")
    }
  }
  
  ind2 = match(df.p2v$variant[-ind], table = df.NSFP$variant_str)
  df.p2v[-ind,score_types] = df.NSFP[ind2,score_types]
  df.p2v$cs_simple[-ind] = df.NSFP$cs_simple[ind2]
  
  # plot score distribution for less than 75% vs more than 75% LDL levels
  ind = which(df.p2v$patient_label == "g1")
  ind2 = which(df.p2v$patient_label == "g2")
  
  # # create data frame to return each score types p-val
  # df.res = c(ldl_thres_g1, ldl_thres_g2)
  # names(df.res) = c("ldl_thres_g1", "ldl_thres_g2")
  
  df.obs = c()
  # draw plot for each score types
  for (i in 1:length(score_types)){
    score_type = score_types[i]
    score_g1 = df.p2v[[score_type]][ind]
    score_g1 = score_g1[which(!is.na(score_g1))]
    score_g1.density = density(score_g1, bw = bw, kernel = kernel)
    score_g2 = df.p2v[[score_type]][ind2]
    score_g2 = score_g2[which(!is.na(score_g2))]
    score_g2.density = density(score_g2, bw = bw, kernel = kernel)
    
    df.obs = c(df.obs, length(score_g1)+length(score_g2))
    
    p1 = NA
    res = ks.test(x = score_g1, y = score_g2, alternative = "greater") # KS-test
    p1 = format(res$p.value, digits = 3)
    
    # # data frame to return each score types p-val
    # score_name = names(score_types)[i]
    # temp = c(length(score_g1), length(score_g2), as.numeric(p1))
    # names(temp) = c(paste0(score_name, "_g1"), paste0(score_name, "_g2"), paste0(score_name, "_pval"))
    # df.res = c(df.res, temp)
    # # ymax = max(c(score_g1.density$y, score_g2.density$y))
    # # xmax = max(abs(c(score_g1.density$x, score_g2.density$x)))
    
    if (is.na(y_labels)){
      y_label = names(score_types)[i]
    } else if (y_labels == F){
      y_label = ""
    } else {
      y_label = y_labels[i]
    }
    
    if (draw_plot){
      if (plot_method == "density"){
        plot(score_g1.density, col = "blue", ylim=c(0,ymax), xlim=x_range, las=1, main = "", xaxs="i", yaxs="i", ylab=y_label, xlab="") # g1 is blue
        lines(score_g2.density, col = "red") # g2 is black
        legend("topleft", legend = c(paste0("Low LDL (n=", length(score_g1), ")"), paste0("High LDL (n=", length(score_g2), ")")), col = c("blue", "red"), lty = 1)
        if (i == 1){
          mtext(title_str, side = 3, outer = T, line = -1.5)
        }
      }
      
      if (plot_method == "histogram"){
        hist(score_g1, breaks = breaks, freq = F, col = add.alpha("blue", alpha = 0.5), ylim=c(0,ymax), xlim=x_range, las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="")
        hist(score_g2, breaks = breaks, freq = F, col = add.alpha("red", alpha = 0.5), ylim=c(0,ymax), xlim=x_range, add=T, las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="")
        legend("topleft", legend = c(paste0("Low LDL (n=", length(score_g1), ")"), paste0("High LDL (n=", length(score_g2), ")")), fill = add.alpha(c("blue", "red"), 0.5))
        box()
      }
      
      text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*ks_text_xpos, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.9, paste0("KS p=", p1))
    }
  }
  
  names(df.obs) = score_types
  return(df.obs)
}

draw_patient_plot_TP53 <- function(df.NSFP, score_types, title_str="", ks_text_xpos=0.8, y_labels = NA, draw_plot=T){
  ymax = 2
  x_range = c(-4,4)
  plot_method = "density"
  bw = 0.25
  kernel = "gaussian"
  
  # load UKBB patient with LFS cancer types
  df.cancer = c()
  fnames = list.files(path = "./UKBB/cancer_types/")
  for (fname in fnames){
    temp = read.csv(paste0("./UKBB/cancer_types/", fname))
    df.cancer = cbind(df.cancer, temp$cancer)
  }
  df.cancer = data.frame(df.cancer)
  colnames(df.cancer) = gsub(pattern = ".csv", replacement = "", x = fnames)
  df.cancer$any_cancer = rowSums(df.cancer)
  df.cancer = cbind(eid = temp$eid, gender = temp$gender, age = temp$age, df.cancer)
  
  ind = which(df.cancer$any_cancer == 0)
  patient_g1 = df.cancer$eid[ind] # g1 is non-cancer group
  ind = which(df.cancer$any_cancer >= 1)
  patient_g2 = df.cancer$eid[ind] # g2 is cancer group
  
  df.p2v = read.csv("./UKBB/TP53_missenseLOF_p2v.csv")
  df.p2v$patient_label = NA
  df.p2v$patient_label[df.p2v$patient %in% patient_g1] = "g1"
  df.p2v$patient_label[df.p2v$patient %in% patient_g2] = "g2"
  
  df.p2v[,score_types] = NA
  df.p2v$cs_simple = NA
  ind = grep(pattern = "\\|", x = df.p2v$variant)
  for (i in ind){
    variant_str = unlist(strsplit(x = df.p2v$variant[i], split = "\\|"))
    ind2 = match(variant_str, table = df.NSFP$variant_str)
    df.p2v[i,score_types] = colMeans(df.NSFP[ind2,score_types], na.rm = T)
    # df.p2v[i,score_types] = apply(df.NSFP[ind2,score_types], 2, median, na.rm = T)
    temp = unique(df.NSFP$cs_simple[ind2])
    temp = temp[!is.na(temp)]
    if (length(temp) == 1){
      df.p2v$cs_simple[i] = paste0(temp, collapse = "")
    }
  }
  
  ind2 = match(df.p2v$variant[-ind], table = df.NSFP$variant_str)
  df.p2v[-ind,score_types] = df.NSFP[ind2,score_types]
  df.p2v$cs_simple[-ind] = df.NSFP$cs_simple[ind2]
  
  # plot score distribution 
  ind = which(df.p2v$patient_label == "g1")
  ind2 = which(df.p2v$patient_label == "g2")
  
  df.obs = c()
  # draw plot for each score types
  for (i in 1:length(score_types)){
    score_type = score_types[i]
    score_g1 = df.p2v[[score_type]][ind]
    score_g1 = score_g1[which(!is.na(score_g1))]
    score_g1.density = density(score_g1, bw = bw, kernel = kernel)
    score_g2 = df.p2v[[score_type]][ind2]
    score_g2 = score_g2[which(!is.na(score_g2))]
    score_g2.density = density(score_g2, bw = bw, kernel = kernel)
    
    df.obs = c(df.obs, length(score_g1)+length(score_g2))
    
    p1 = NA
    res = ks.test(x = score_g1, y = score_g2, alternative = "greater") # KS-test
    p1 = format(res$p.value, digits = 3)
    
    # ymax = max(c(score_g1.density$y, score_g2.density$y))
    # xmax = max(abs(c(score_g1.density$x, score_g2.density$x)))
    
    if (is.na(y_labels)){
      y_label = names(score_types)[i]
    } else if (y_labels == F){
      y_label = ""
    } else {
      y_label = y_labels[i]
    }
    
    if (draw_plot){
      if (plot_method == "density"){
        plot(score_g1.density, col = "blue", ylim=c(0,ymax), xlim=x_range, las=1, main = "", xaxs="i", yaxs="i", ylab=y_label, xlab="") # g1 is blue
        lines(score_g2.density, col = "red") # g2 is black
        legend("topleft", legend = c(paste0("non-LFS (n=", length(score_g1), ")"), paste0("LFS (n=", length(score_g2), ")")), col = c("blue", "red"), lty = 1)
        if (i == 1){
          mtext(title_str, side = 3, outer = T, line = -1.5)
        }
      }
      
      if (plot_method == "histogram"){
        hist(score_g1, breaks = breaks, freq = F, col = add.alpha("blue", alpha = 0.5), ylim=c(0,ymax), xlim=x_range, las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="")
        hist(score_g2, breaks = breaks, freq = F, col = add.alpha("red", alpha = 0.5), ylim=c(0,ymax), xlim=x_range, add=T, las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="")
        legend("topleft", legend = c(paste0("non-LFS (n=", length(score_g1), ")"), paste0("LFS (n=", length(score_g2), ")")), fill = add.alpha(c("blue", "red"), 0.5))
        box()
      }
      
      text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*ks_text_xpos, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.9, paste0("KS p=", p1))
    }
  }
  
  names(df.obs) = score_types
  return(df.obs)
}

draw_score_distribution <- function(df, probability=F){
  temp = hist(df$norm_raw_score, breaks = 100, probability = probability)
  ind = which(df$functional_class == "SYN")
  hist(df$norm_raw_score[ind], breaks = temp$breaks, probability = probability, col = add.alpha("blue", 0.5), add=T)
  ind = which(df$functional_class == "LOF")
  hist(df$norm_raw_score[ind], breaks = temp$breaks, probability = probability, col = add.alpha("red", 0.5), add=T)
  ind = which(df$aaalt == "P" | df$aaref == "P")
  hist(df$norm_raw_score[ind], breaks = temp$breaks, probability = probability, col = add.alpha("purple", 0.5), add=T)
}

normalize_by_median_sd <- function(score, invert_score=F){
  if (invert_score){
    score_out = -score
  } else {
    score_out = score
  }
  
  score_out = (score_out - median(score_out, na.rm=T))/sd(score_out, na.rm=T)
  return(score_out)
}

prepare_dfNSFP_BRCA1 <- function(data = "BEHive", df.funsum, df.NSFP_all, AF_thres = 0.005, edit_thres = 0, draw_plot=F, ls.funsum_ss = NULL, df.clinvar = NULL, include_LOF=T){
  # load BRCA1 functional data
  if (grepl(pattern = "BEHive", x = data, ignore.case = T)){
    if (grepl(pattern = "CBE", x = data, ignore.case = T)){
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE_BEHive.tsv", header = T, sep = "\t")
      ind = which(df.example$edit_frequency > edit_thres)
      df = df.example[ind,]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else if (grepl(pattern = "ABE", x = data, ignore.case = T)){
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE_BEHive.tsv", header = T, sep = "\t")
      ind = which(df.example$edit_frequency > edit_thres)
      df = df.example[ind,]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else {
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE_BEHive.tsv", header = T, sep = "\t")
      ind = which(df.example$edit_frequency > edit_thres)
      df = df.example[ind,]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
      
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE_BEHive.tsv", header = T, sep = "\t")
      ind = which(df.example$edit_frequency > edit_thres)
      df2 = df.example[ind,]
      df2$norm_raw_score = normalize_by_median_sd(df2$raw_score, invert_score = T)
      
      df = rbind(df, df2)
    }
  } else if (grepl(pattern = "original", x = data, ignore.case = T)) {
    if (grepl(pattern = "CBE", x = data, ignore.case = T)){
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else if (grepl(pattern = "ABE", x = data, ignore.case = T)){
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else {
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-CBE.tsv", header = T, sep = "\t")
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    
      df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1-NG-ABE.tsv", header = T, sep = "\t")
      df2 = df.example
      df2$norm_raw_score = normalize_by_median_sd(df2$raw_score, invert_score = T)
      
      df = rbind(df, df2)
    }
  } else if (data == "DMS"){
    df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/BRCA1_DMS_v3.tsv", header = T, sep = "\t")
    # temp = check_direction(df = df.example, normalize = T)
    # df = temp[[1]]
    df = df.example
    df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
  } else {
    errorCondition("data seelction is invalid!")
  }
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  ## de-noise DMS data for each gene
  if (is.null(ls.funsum_ss)){
    df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = include_LOF)
  } else {
    df.out = de_noise_ss(df = df, pos_mean_method = "js", df.funsum = df.funsum, ls.funsum_ss = ls.funsum_ss, include_LOF = include_LOF)
  }
  
  ## load df.NSFP for each gene
  gene_id = "BRCA1"
  ind = which(df.NSFP_all$genename == gene_id)
  df.NSFP = df.NSFP_all[ind,]
  df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
  df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
  df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
  ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
  df.NSFP = df.NSFP[ind,]
  
  ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
  df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
  df.NSFP$pos_score = df.out$pos_score[ind]
  df.NSFP$sub_score = df.out$sub_score[ind]
  df.NSFP$final_score = df.out$final_score[ind]
  df.NSFP$final_score_lite = df.NSFP$final_score
  ind2 = which(is.na(df.NSFP$norm_raw_score))
  if (length(ind2)){
    df.NSFP$final_score_lite[ind2] = NA
  }
  
  if (!is.null(ls.funsum_ss)){
    df.NSFP$sub_score_ss = df.out$sub_score_ss[ind]
    df.NSFP$final_score_ss = df.out$final_score_ss[ind]
    df.NSFP$final_score_ss_lite = df.NSFP$final_score_ss
    ind2 = which(is.na(df.NSFP$norm_raw_score))
    if (length(ind2)){
      df.NSFP$final_score_ss_lite[ind2] = NA
    }
  }
  
  if (!is.null(df.clinvar)){
    df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
    ind = match(x = df.NSFP$gene_aa_str, table = df.clinvar$gene_aa_str)
    cs_temp = df.clinvar$cs_simple[ind]
    df.NSFP$cs_simple[cs_temp == "p"] = "p"
    df.NSFP$cs_simple[cs_temp == "b"] = "b"
    
    # add gold start ratings
    df.NSFP$gold_star = df.clinvar$gold_star[ind]
  }
  
  if (!include_LOF){
    ind = which((df.NSFP$aaref != "*") & (df.NSFP$aaalt != "*"))
    df.NSFP = df.NSFP[ind,]
  }
  
  return(df.NSFP)
}

prepare_dfNSFP_LDLR <- function(data = "BEHive", df.funsum, df.NSFP_all, AF_thres = 0.005, edit_thres = 0, draw_plot=F, df.clinvar = NULL, include_LOF=T){
  # load LDLR functional data
  # ABE does not have stop gain so cannot check direction
  if (data == "BEHive"){
    df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/LDLR-ABE-HEK293T_BEHive.tsv", header = T, sep = "\t")
    ind = which(df.example$edit_frequency > edit_thres)
    df = df.example[ind,]
    df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
  } else if (data == "DMS") {
    df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/LDLR_DMS_v3.tsv", header = T, sep = "\t")
    df = df.example
    df$norm_raw_score = normalize_by_median_sd(df$raw_score)
  } else {
    errorCondition("data seelction is invalid!")
  }
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  ## de-noise DMS data for each gene
  df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = include_LOF)
  
  ## load df.NSFP for each gene
  gene_id = "LDLR"
  ind = which(df.NSFP_all$genename == gene_id)
  df.NSFP = df.NSFP_all[ind,]
  df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
  df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
  df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
  ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
  df.NSFP = df.NSFP[ind,]
  
  ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
  df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
  df.NSFP$pos_score = df.out$pos_score[ind]
  df.NSFP$sub_score = df.out$sub_score[ind]
  df.NSFP$final_score = df.out$final_score[ind]
  df.NSFP$final_score_lite = df.NSFP$final_score
  ind = which(is.na(df.NSFP$norm_raw_score))
  if (length(ind)){
    df.NSFP$final_score_lite[ind] = NA
  }
  
  if (!is.null(df.clinvar)){
    df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
    ind = match(x = df.NSFP$gene_aa_str, table = df.clinvar$gene_aa_str)
    cs_temp = df.clinvar$cs_simple[ind]
    df.NSFP$cs_simple[cs_temp == "p"] = "p"
    df.NSFP$cs_simple[cs_temp == "b"] = "b"
    
    # add gold start ratings
    df.NSFP$gold_star = df.clinvar$gold_star[ind]
  }
  
  if (!include_LOF){
    ind = which((df.NSFP$aaref != "*") & (df.NSFP$aaalt != "*"))
    df.NSFP = df.NSFP[ind,]
  }
  
  return(df.NSFP)
}

prepare_dfNSFP_TP53 <- function(data = "DMS", df.funsum, df.NSFP_all, AF_thres = 0.005, edit_thres = 0, draw_plot=F, ls.funsum_ss = NULL, df.clinvar = NULL, include_LOF=T){
  # load TP53 functional data
  if (data == "DMS") {
    df.example = read.table("./shiny_app/DMS_denoiser/DMS_data/TP53_DMS_v3.tsv", header = T, sep = "\t")
    df = df.example
    # temp = check_direction(df)
    # df = temp[[1]]
    df$norm_raw_score = normalize_by_median_sd(df$raw_score)
  } else {
    errorCondition("data seelction is invalid!")
  }
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  ## de-noise DMS data for each gene
  if (is.null(ls.funsum_ss)){
    df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = include_LOF)
  } else {
    df.out = de_noise_ss(df = df, pos_mean_method = "js", df.funsum = df.funsum, ls.funsum_ss = ls.funsum_ss, include_LOF = include_LOF)
  }
  
  ## load df.NSFP for each gene
  gene_id = "TP53"
  ind = which(df.NSFP_all$genename == gene_id)
  df.NSFP = df.NSFP_all[ind,]
  df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
  df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
  df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
  ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
  df.NSFP = df.NSFP[ind,]
  
  ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
  df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
  df.NSFP$pos_score = df.out$pos_score[ind]
  df.NSFP$sub_score = df.out$sub_score[ind]
  df.NSFP$final_score = df.out$final_score[ind]
  df.NSFP$final_score_lite = df.NSFP$final_score
  ind2 = which(is.na(df.NSFP$norm_raw_score))
  if (length(ind2)){
    df.NSFP$final_score_lite[ind2] = NA
  }
  
  if (!is.null(ls.funsum_ss)){
    df.NSFP$sub_score_ss = df.out$sub_score_ss[ind]
    df.NSFP$final_score_ss = df.out$final_score_ss[ind]
    df.NSFP$final_score_ss_lite = df.NSFP$final_score_ss
    ind2 = which(is.na(df.NSFP$norm_raw_score))
    if (length(ind2)){
      df.NSFP$final_score_ss_lite[ind2] = NA
    }
  }
  
  if (!is.null(df.clinvar)){
    df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
    ind = match(x = df.NSFP$gene_aa_str, table = df.clinvar$gene_aa_str)
    cs_temp = df.clinvar$cs_simple[ind]
    df.NSFP$cs_simple[cs_temp == "p"] = "p"
    df.NSFP$cs_simple[cs_temp == "b"] = "b"
    
    # add gold start ratings
    df.NSFP$gold_star = df.clinvar$gold_star[ind]
  }
  
  if (!include_LOF){
    ind = which((df.NSFP$aaref != "*") & (df.NSFP$aaalt != "*"))
    df.NSFP = df.NSFP[ind,]
  }
  
  return(df.NSFP)
}

prepare_dbNSFP_cancer_genes <- function(data = "original", df.funsum, df.NSFP_all, AF_thres = 0.005, draw_plot=F, include_LOF=T){
  if (grepl(pattern = "cell_prolif", x = data, ignore.case = T)) {
    df.example = read.csv("./shiny_app/DMS_denoiser/DMS_data/Cancer_genes_cell_proliferation_ABE_CBE.csv", header = T)
    
    if (grepl(pattern = "absLFC", x = data, ignore.case = T)){
      df.example$raw_score = abs(df.example$raw_score)
    }
    
    if (grepl(pattern = "CBE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "CBE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score)
    } else if (grepl(pattern = "ABE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "ABE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score)
    } else {
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score)
    }
  } else if (grepl(pattern = "egf_dependency", x = data, ignore.case = T)) {
    df.example = read.csv("./shiny_app/DMS_denoiser/DMS_data/Cancer_genes_EGF_dependency_ABE_CBE.tsv", header = T)
    
    if (grepl(pattern = "absLFC", x = data, ignore.case = T)){
      df.example$raw_score = abs(df.example$raw_score)
    }
    
    if (grepl(pattern = "CBE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "CBE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score)
    } else if (grepl(pattern = "ABE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "ABE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score)
    } else {
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score)
    }
  } else {
    errorCondition("data seelction is invalid!")
  }
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  ## de-noise DMS data for each gene
  df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = include_LOF)
  
  ## load df.NSFP for each gene
  gene_ids = unique(df.out$gene)
  df.NSFP = df.NSFP_all[df.NSFP_all$genename %in% gene_ids,]
  df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
  df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
  df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
  ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
  df.NSFP = df.NSFP[ind,]
  
  ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
  df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
  df.NSFP$pos_score = df.out$pos_score[ind]
  df.NSFP$sub_score = df.out$sub_score[ind]
  df.NSFP$final_score = df.out$final_score[ind]
  df.NSFP$final_score_lite = df.NSFP$final_score
  ind = which(is.na(df.NSFP$norm_raw_score))
  if (length(ind)){
    df.NSFP$final_score_lite[ind] = NA
  }
  
  if (!include_LOF){
    ind = which((df.NSFP$aaref != "*") & (df.NSFP$aaalt != "*"))
    df.NSFP = df.NSFP[ind,]
  }
  
  return(df.NSFP)
}

prepare_dbNSFP_DDR_genes <- function(df.funsum, df.NSFP_all, AF_thres = 0.005, draw_plot=F, df.clinvar = NULL, include_LOF=T){
  df.example = read.csv("./shiny_app/DMS_denoiser/DMS_data/DDR_genes_BE.csv", header = T)
  df = df.example
  df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  ## de-noise DMS data for each gene
  df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = include_LOF)
  
  ## load df.NSFP for each gene
  gene_ids = unique(df.out$gene)
  df.NSFP = df.NSFP_all[df.NSFP_all$genename %in% gene_ids,]
  df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
  df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
  df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
  ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
  df.NSFP = df.NSFP[ind,]
  
  ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
  df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
  df.NSFP$pos_score = df.out$pos_score[ind]
  df.NSFP$sub_score = df.out$sub_score[ind]
  df.NSFP$final_score = df.out$final_score[ind]
  df.NSFP$final_score_lite = df.NSFP$final_score
  ind = which(is.na(df.NSFP$norm_raw_score))
  if (length(ind)){
    df.NSFP$final_score_lite[ind] = NA
  }
  
  ind = match(x = df.NSFP$gene_aa_str, table = df$gene_aa_str)
  cs_temp = df$cs_simple[ind]
  df.NSFP$cs_simple[cs_temp == "p"] = "p"
  df.NSFP$cs_simple[cs_temp == "b"] = "b"
  
  if (!is.null(df.clinvar)){
    df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
    ind = match(x = df.NSFP$gene_aa_str, table = df.clinvar$gene_aa_str)
    cs_temp = df.clinvar$cs_simple[ind]
    df.NSFP$cs_simple[cs_temp == "p"] = "p"
    df.NSFP$cs_simple[cs_temp == "b"] = "b"
    
    # add gold start ratings
    df.NSFP$gold_star = df.clinvar$gold_star[ind]
  }
  
  if (!include_LOF){
    ind = which((df.NSFP$aaref != "*") & (df.NSFP$aaalt != "*"))
    df.NSFP = df.NSFP[ind,]
  }
  
  return(df.NSFP)
}

prepare_df_cancer_genes <- function(data = "original", df.funsum, draw_plot=F){
  if (grepl(pattern = "cell_prolif", x = data, ignore.case = T)) {
    df.example = read.csv("./shiny_app/DMS_denoiser/DMS_data/Cancer_genes_cell_proliferation_ABE_CBE.csv", header = T)
    
    if (grepl(pattern = "absLFC", x = data, ignore.case = T)){
      df.example$raw_score = abs(df.example$raw_score)
    }
    
    if (grepl(pattern = "CBE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "CBE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else if (grepl(pattern = "ABE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "ABE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else {
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    }
  } else if (grepl(pattern = "egf_dependency", x = data, ignore.case = T)) {
    df.example = read.csv("./shiny_app/DMS_denoiser/DMS_data/Cancer_genes_EGF_dependency_ABE_CBE.tsv", header = T)
    
    if (grepl(pattern = "absLFC", x = data, ignore.case = T)){
      df.example$raw_score = abs(df.example$raw_score)
    }
    
    if (grepl(pattern = "CBE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "CBE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else if (grepl(pattern = "ABE", x = data, ignore.case = T)){
      df = df.example[df.example$base_editor == "ABE",]
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    } else {
      df = df.example
      df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
    }
  } else {
    errorCondition("data seelction is invalid!")
  }
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  # simplify classification
  df$class_simple = NA
  ind = which(df$classification == "Depleting")
  df$class_simple[ind] = "Depleting"
  ind = which(df$classification == "Outgrowing")
  df$class_simple[ind] = "Outgrowing"
  ind = which(df$classification == "Likely depleting")
  df$class_simple[ind] = "Likely depleting"
  ind = which(df$classification == "Likely outgrowing")
  df$class_simple[ind] = "Likely outgrowing"
  ind = which(df$classification == "Likely neutral (possibly depleting)")
  df$class_simple[ind] = "Possibly depleting"
  ind = which(df$classification == "Likely neutral (possibly outgrowing)")
  df$class_simple[ind] = "Possibly outgrowing"
  
  ## de-noise DMS data for each gene
  df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = T)
  df$gene_aa_str = paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)
  ind = match(x = df$gene_aa_str, table = df.out$gene_aa_str)
  df$pos_score = df.out$pos_score[ind]
  df$sub_score = df.out$sub_score[ind]
  df$final_score = df.out$final_score[ind]
  df$final_score_lite = df$final_score
  ind = which(is.na(df$norm_raw_score))
  if (length(ind)){
    df$final_score_lite[ind] = NA
  }
  
  return(df)
}

prepare_df_DDR_genes <- function(df.funsum, AF_thres = 0.005, draw_plot=F){
  df.example = read.csv("./shiny_app/DMS_denoiser/DMS_data/DDR_genes_BE.csv", header = T)
  df = df.example
  df$norm_raw_score = normalize_by_median_sd(df$raw_score, invert_score = T)
  
  if (draw_plot){
    draw_score_distribution(df)
  }
  
  # creat cs_simple
  df$cs_simple = NA
  ind = which(df$ClinVar_from_data == "P/LP")
  df$cs_simple[ind] = "p"
  ind = which(df$ClinVar_from_data == "B/LB")
  df$cs_simple[ind] = "b"
  
  ## de-noise DMS data for each gene
  df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum, include_LOF = T)
  
  ## since there's no df.NSFP for all DDR genes, copy the cs_simple from the original DDR dataset to df.out
  df$gene_aa_str = paste0(df$gene, "---", df$aaref, df$aapos, df$aaalt)
  ind = match(x = df.out$gene_aa_str, table = df$gene_aa_str)
  df.out$cs_simple = df$cs_simple[ind]
  df.out$final_score_lite = df.out$final_score
  ind = which(is.na(df.out$norm_raw_score))
  if (length(ind)){
    df.out$final_score_lite[ind] = NA
  }
  
  # remove aa_strs not included in original dataset
  df.out = df.out[-ind,]
  return(df.out)
}

prepare_dfNSFP_aggregate <- function(df.clinvar = NULL, include_LOF=T){
  # load df.NSFP for 14 genes
  df.NSFP_aggregate = read.xlsx("./shiny_app/DMS_denoiser/dfNFSP_14genes_050622.xlsx")

  # load raw scores from 12 genes in standard format
  ls.DMS_score = readRDS("DMS_score_14genes.rds")

  # calculate FUMSUM using a set of included genes' functional data
  gene_ids = names(ls.DMS_score)
  gene_ids.exclude = c("LDLR", "SUMO1", "CALM1", "TPK1")
  gene_ids.no_invert = c("TP53", "MAPK1")
  AF_thres = 0.005

  df.NSFP_out = c()
  for (gene_id.test in gene_ids){
    gene_ids.pool = gene_ids[!(gene_ids %in% c(gene_ids.exclude, gene_id.test))]
    df.all_sub_score = c()
    for (gene_id in gene_ids.pool){
      df = ls.DMS_score[[gene_id]]
      if (gene_id %in% gene_ids.no_invert){
        df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF)
      } else {
        df.sub_score = get_norm_score(df, pos_mean_method = "js", invert_score = T, include_LOF = include_LOF)
      }

      df.all_sub_score = rbind(df.all_sub_score, cbind(gene_name = gene_id, df.sub_score))
    }
    df.funsum_all = get_FUNSUM(df.all_sub_score[,-1], avg_method = "js", include_LOF = T) # get combined funsum

    # normalize data by median and sd
    df = ls.DMS_score[[gene_id.test]]
    ifelse(gene_id.test %in% gene_ids.no_invert, {df$norm_raw_score = df$raw_score}, {df$norm_raw_score = -df$raw_score})
    df$norm_raw_score = (df$norm_raw_score - median(df$norm_raw_score, na.rm = T))/sd(df$norm_raw_score, na.rm = T)

    ## de-noise DMS data for each gene
    df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF)

    ## load df.NSFP for each gene
    ind = which(df.NSFP_aggregate$genename == gene_id.test)
    df.NSFP = df.NSFP_aggregate[ind,]
    df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
    df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
    df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
    ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
    df.NSFP = df.NSFP[ind,]

    ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
    df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
    df.NSFP$pos_score = df.out$pos_score[ind]
    df.NSFP$sub_score = df.out$sub_score[ind]
    df.NSFP$final_score = df.out$final_score[ind]
    df.NSFP$final_score_lite = df.NSFP$final_score
    ind = which(is.na(df.NSFP$norm_raw_score))
    if (length(ind)){
      df.NSFP$final_score_lite[ind] = NA
    }

    df.NSFP_out = rbind(df.NSFP_out, df.NSFP)
  }

  df.NSFP = df.NSFP_out
  if (!is.null(df.clinvar)){
    df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
    ind = match(x = df.NSFP$gene_aa_str, table = df.clinvar$gene_aa_str)
    cs_temp = df.clinvar$cs_simple[ind]
    df.NSFP$cs_simple[cs_temp == "p"] = "p"
    df.NSFP$cs_simple[cs_temp == "b"] = "b"

    # add gold start ratings
    df.NSFP$gold_star = df.clinvar$gold_star[ind]
  }
  
  if (!include_LOF){
    ind = which((df.NSFP$aaref != "*") & (df.NSFP$aaalt != "*"))
    df.NSFP = df.NSFP[ind,]
  }

  return(df.NSFP)
}



draw_FUNSUM_heatmap <- function(df.funsum, title_str="FUNSUM", include_LOF=F, show_diag=T, aa_colored=T){
  require(tidyverse)
  require(hrbrthemes)
  
  # color the amino acids according to their class
  if (include_LOF){
    aa_order = unlist(str_split("*CSTAGPDEQNHRKMILVWYF", pattern = ""))
    aa_col = c("black", "purple", rep("blue",5), rep("green",4), rep("orange",3), rep("red",4), rep("pink",3))
  } else {
    aa_order = unlist(str_split("CSTAGPDEQNHRKMILVWYF", pattern = ""))
    aa_col = c("purple", rep("blue",5), rep("green",4), rep("orange",3), rep("red",4), rep("pink",3))
  }
  
  # convert funsum to long format
  df.funsum_long = df.funsum %>% as_tibble(rownames = "aaref") %>% 
    pivot_longer(cols = !aaref, names_to = "aaalt", values_to = "score") %>%
    filter((aaref %in% aa_order) & (aaalt %in% aa_order))
  
  # normalize
  score_abs_max = max(abs(df.funsum_long$score), na.rm=T)
  df.funsum_long$score = df.funsum_long$score/score_abs_max # scale score to [-1,1] but not changing median&mean
  # df.funsum_long$score = (df.funsum_long$score - median(df.funsum_long$score, na.rm=T))/sd(df.funsum_long$score, na.rm=T) # scale score by median=0 and sd=1
  # df.funsum_long$score = (df.funsum_long$score)/sd(df.funsum_long$score, na.rm=T) # scale score by sd=1
  
  if (!show_diag){
    ind = which(df.funsum_long$aaref == df.funsum_long$aaalt)
    df.funsum_long$score[ind] = NaN
  }
  
  # draw heatmap
  plt.funsum = ggplot(data = df.funsum_long, mapping = aes(aaref, aaalt, fill= score)) + geom_tile(color="black") + 
    scale_fill_gradient2(low="blue", mid = "white", high="red", limits=c(-1,1), na.value = "grey75") + coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10)) + 
    labs(title=title_str, x ="Reference amino acid", y = "Alternative amino acid") + 
    xlim(aa_order) + ylim(aa_order) + 
    theme_ipsum(axis_title_size=12, plot_title_size=15)
  
  if (aa_colored){
    plt.funsum = plt.funsum + theme(axis.text = element_text(colour = aa_col))
  }

  return(plt.funsum)
}


load_maveDB_csv <- function(file_path){
  require(Biostrings)
  require(tidyverse)
  AMINO_ACID_CODE_v2 = c(AMINO_ACID_CODE, "*"="Ter", "="="=")
  df.raw = read_csv(file_path, skip = 4)
  raw_score = df.raw$score
  df.raw$temp = gsub(pattern = "p.|\\[|\\]", replacement = "", df.raw$hgvs_pro)
  temp = str_split(df.raw$temp, pattern = "\\d+", simplify = T)
  if (ncol(temp) > 1){
    ind = match(temp[,1], table = AMINO_ACID_CODE_v2)
    aaref = names(AMINO_ACID_CODE_v2)[ind]
    ind = match(temp[,2], table = AMINO_ACID_CODE_v2)
    aaalt = names(AMINO_ACID_CODE_v2)[ind]
    ind = which(aaalt == "=")
    aaalt[ind] = aaref[ind]
    aapos = as.numeric(gsub(pattern = "\\D+", replacement = "", x = df.raw$temp))
    df.out = data.frame(aapos, aaref, aaalt, functional_class = NA, raw_score)
    ind = which(df.out$aaref != df.out$aaalt & df.out$aaalt != '*')
    df.out$functional_class[ind] = "MIS"
    ind = which(df.out$aaref != "*" & df.out$aaalt == "*")
    df.out$functional_class[ind] = "LOF"
    ind = which(df.out$aaref == df.out$aaalt)
    df.out$functional_class[ind] = "SYN"
    df.out = df.out %>% drop_na()
  } else {
    df.out = NULL
  }
  
  return(df.out)
}

draw_score_distribution_ggplot <- function(df, score_type = "raw_score", bin_ct = 30){
  require(ggridges)
  
  plt = df %>% ggplot(aes(x=.data[[score_type]], y=functional_class, fill=functional_class)) +
    geom_density_ridges(alpha=0.6, stat="binline", bins = bin_ct) + ylab("") + xlab(score_type) +  
    theme_ridges()
  return(plt)
}

normalize_scoreset <- function(df, lower_bound = 0.05, upper_bound = 0.95, force_quantile=F){
  if (!"functional_class" %in% colnames(df)){
    df$functional_class = NA
    ind = which(df$aaref != df$aaalt & df$aaalt != '*')
    df$functional_class[ind] = "MIS"
    ind = which(df$aaref != "*" & df$aaalt == "*")
    df$functional_class[ind] = "LOF"
    ind = which(df$aaref == df$aaalt)
    df$functional_class[ind] = "SYN"
  }
  
  ind_lof = which(df$functional_class == "LOF")
  if (length(ind_lof)>10){
    median_lof = median(df$raw_score[ind_lof], na.rm = T)
  } else {
    median_lof = NA
  }
  ind_mis = which(df$functional_class == "MIS")
  median_mis = median(df$raw_score[ind_mis], na.rm = T)
  ind_syn = which(df$functional_class == "SYN")
  if (length(ind_syn)>10){
    median_syn = median(df$raw_score[ind_syn], na.rm = T)
  } else {
    median_syn = NA
  }
  
  temp = df$raw_score
  if (!is.na(median_syn) & !is.na(median_lof)){
    if (median_lof < median_syn){
      temp = -temp
      median_lof = -median_lof
      median_syn = -median_syn
    }
    df$norm_raw_score = (temp-median_syn)/(median_lof-median_syn)
  } else {
    mis_scores = df$raw_score[ind_mis]
    if (!is.na(median_syn) & is.na(median_lof)){
      if (median_mis < median_syn){
        temp = -temp
        median_syn = -median_syn
        mis_scores = -mis_scores
      }
      median_lof = quantile(mis_scores, upper_bound, na.rm=T)
      df$norm_raw_score = (temp-median_syn)/(median_lof-median_syn)
    } else if (is.na(median_syn) & !is.na(median_lof)){
      if (median_mis > median_lof){
        temp = -temp
        median_lof = -median_lof
        mis_scores = -mis_scores
      }
      median_syn = quantile(mis_scores, lower_bound, na.rm=T)
      df$norm_raw_score = (temp-median_syn)/(median_lof-median_syn)
    } else {
      ind_p = which(df$aaalt == "P")
      median_p = median(df$raw_score[ind_p], na.rm=T)
      if (median_p < median_mis){
        temp = -temp
        median_syn = -median_syn
        mis_scores = -mis_scores
      }
      median_syn = quantile(mis_scores, lower_bound, na.rm = T)
      median_lof = quantile(mis_scores, upper_bound, na.rm = T)
      df$norm_raw_score = (temp-median_syn)/(median_lof-median_syn)
    }
  }
  
  if (force_quantile){
    ind_mis = which(df$functional_class == "MIS")
    temp = df$norm_raw_score
    mis_scores = df$norm_raw_score[ind_mis]
    median_syn = quantile(mis_scores, lower_bound, na.rm = T)
    median_lof = quantile(mis_scores, upper_bound, na.rm = T)
    df$norm_raw_score = (temp-median_syn)/(median_lof-median_syn)
  }
  return(df)
}



clinvar_align_aapos <- function(df_gene, df.clinvar_gene, test_range = c(-5, 5), diagnose=F){
  require(tidyverse)
  df.align_res = c()
  for (i in seq(test_range[1], test_range[2])){
    df_gene_2 = df_gene %>% 
      mutate(test_aapos = aapos+i, aapos_aaref = paste0(test_aapos, aaref)) %>%
      group_by(aapos_aaref) %>%
      summarise(test_aapos = unique(test_aapos)) %>% 
      arrange(test_aapos)
    
    df.clinvar_gene_2 = df.clinvar_gene %>%
      mutate(aapos_aaref = paste0(aapos, aaref)) %>%
      group_by(aapos_aaref) %>%
      summarise(aapos = unique(aapos))
    
    ind = match(x = df_gene_2$aapos_aaref, table = df.clinvar_gene_2$aapos_aaref)
    df_gene_2$clinvar_aapos_aaref = df.clinvar_gene_2$aapos_aaref[ind]
    temp = data.frame(i, matched_ct = sum(!is.na(ind)))
    df.align_res = rbind(df.align_res, temp)
  }
  
  if (diagnose){
    print(df.align_res)
  }
  
  ind = which(df.align_res$matched_ct == max(df.align_res$matched_ct))
  pos_diff = df.align_res$i[ind]
  if (0 %in% pos_diff){
    return(0)
  } else if (length(pos_diff)>0){
    ind = which(abs(pos_diff) == min(abs(pos_diff)))
    return(pos_diff[ind][1])
  }
}





# prepare_dfNSFP_aggregate <- function(df.clinvar = NULL){
#   # load df.NSFP for 14 genes
#   df.NSFP_aggregate = read.xlsx("./shiny_app/DMS_denoiser/dfNFSP_14genes_050622.xlsx")
#   
#   # load raw scores from 12 genes in standard format
#   ls.DMS_score = readRDS("DMS_score_14genes.rds")
#   
#   # calculate FUMSUM using a set of included genes' functional data
#   gene_ids = names(ls.DMS_score)
#   gene_ids.exclude = c("LDLR") # LDLR excluded, but all other genes okay
#   gene_ids.no_invert = c("TP53", "MAPK1")
#   AF_thres = 0.005
#   
#   df.NSFP_out = c()
#   gene_ids = gene_ids[!(gene_ids %in% gene_ids.exclude)]
#   for (gene_id.test in gene_ids){
#     gene_ids.pool = gene_ids[!(gene_ids %in% gene_id.test)]
#     df.all_sub_score = c()
#     for (gene_id in gene_ids.pool){
#       df = ls.DMS_score[[gene_id]]
#       if (gene_id %in% gene_ids.no_invert){
#         df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = T)
#       } else {
#         df.sub_score = get_norm_score(df, pos_mean_method = "js", invert_score = T, include_LOF = T)
#       }
#       
#       df.all_sub_score = rbind(df.all_sub_score, cbind(gene_name = gene_id, df.sub_score))
#     }
#     df.funsum_all = get_FUNSUM(df.all_sub_score[,-1], avg_method = "js", include_LOF = T) # get combined funsum
#     
#     # normalize data by median and sd
#     df = ls.DMS_score[[gene_id.test]]
#     ifelse(gene_id.test %in% gene_ids.no_invert, {df$norm_raw_score = df$raw_score}, {df$norm_raw_score = -df$raw_score})
#     df$norm_raw_score = (df$norm_raw_score - median(df$norm_raw_score, na.rm = T))/sd(df$norm_raw_score, na.rm = T)
#     
#     ## de-noise DMS data for each gene
#     df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = T)
#     
#     ## load df.NSFP for each gene
#     ind = which(df.NSFP_aggregate$genename == gene_id.test)
#     df.NSFP = df.NSFP_aggregate[ind,]
#     df.NSFP$gene_aa_str = paste0(df.NSFP$genename, "---", df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
#     df.NSFP$aa_pair = paste0(df.NSFP$aaref, df.NSFP$aaalt)
#     df.NSFP$variant_str = paste(df.NSFP$chr, df.NSFP$pos, df.NSFP$ref, df.NSFP$alt, sep = "-")
#     ind = which(df.NSFP$gnomAD_exomes_AF < AF_thres | is.na(df.NSFP$gnomAD_exomes_AF)) # filter AF<0.005
#     df.NSFP = df.NSFP[ind,]
#     
#     ind = match(df.NSFP$gene_aa_str, df.out$gene_aa_str)
#     df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
#     df.NSFP$pos_score = df.out$pos_score[ind]
#     df.NSFP$sub_score = df.out$sub_score[ind]
#     df.NSFP$final_score = df.out$final_score[ind]
#     df.NSFP$final_score_lite = df.NSFP$final_score
#     ind = which(is.na(df.NSFP$norm_raw_score))
#     if (length(ind)){
#       df.NSFP$final_score_lite[ind] = NA
#     }
#     
#     df.NSFP_out = rbind(df.NSFP_out, df.NSFP)
#   }
#   
#   df.NSFP = df.NSFP_out
#   if (!is.null(df.clinvar)){
#     df.clinvar$gene_aa_str = paste0(df.clinvar$gene, "---", df.clinvar$aaref, df.clinvar$aapos, df.clinvar$aaalt)
#     ind = match(x = df.NSFP$gene_aa_str, table = df.clinvar$gene_aa_str)
#     cs_temp = df.clinvar$cs_simple[ind]
#     df.NSFP$cs_simple[cs_temp == "p"] = "p"
#     df.NSFP$cs_simple[cs_temp == "b"] = "b"
#     
#     # add gold start ratings
#     df.NSFP$gold_star = df.clinvar$gold_star[ind]
#   }
#   
#   return(df.NSFP)
# }





