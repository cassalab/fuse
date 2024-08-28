library(shiny)
library(shinyjs)
library(shinydashboard)
library(tidyverse)
library(patchwork)
library(hrbrthemes)
library(vroom)
library(DT)

# df.af = af_predictions(uniprot_ids = "Q9BYF1")



## global variables
version_number = 'v1.0.0'
df.example = read.table("./example_DMS.txt", header = T, sep = "\t")

include_LOF = F
ss_list = c("G" = "Helices", "H" = "Helices", "I" = "Helices", 
            "E" = "Strands", "B" = "Strands", 
            "T" = "Loops", "S" = "Loops", "C" = "Loops")
df.funsum_all = readRDS("./funsum/funsum_combined_noLOF_051624.rds")
ls.funsum_ss = readRDS("./funsum/funsum_SS_combined_noLOF_051624.rds")
example_txt = paste(apply(df.example[1:5,], MARGIN = 1, FUN = paste, collapse="\t"), collapse = "\n")



## functions
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

get_js <- function(m){
  mbar = mean(colMeans(m, na.rm = T), na.rm = T) # global mean
  mu0 = colMeans(m, na.rm = T) # column means
  s2 = var(as.vector(m), na.rm = T)/(nrow(m)*ncol(m)-length(which(is.na(as.vector(m)) == T))) # global variance
  cval = 1 - (ncol(m)-2)*s2/sum((mu0 - mbar)^2, na.rm = T) # adjusted value
  js_est = mbar + cval*(mu0 - mbar)
  return(js_est)
}

funsum_to_subTable <- function(df.funsum){
  df.sub_tb = c()
  for (i in 1:nrow(df.funsum)){
    temp = data.frame(aa_pair = paste0(rownames(df.funsum)[i], colnames(df.funsum)), score = as.vector(df.funsum[i,]))
    df.sub_tb = rbind(df.sub_tb, temp)
  }
  
  return(df.sub_tb)
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
  
  for (i in 1:nrow(df.pos_score)){
    ind = which(df$gene == df.pos_score$gene[i] & df$aapos == df.pos_score$aapos[i])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$norm_raw_score[ind[ind2]]
    
    # calculate pos_mean by funsum method
    if (pos_mean_method == "funsum"){
      ind = which(!is.na(temp[i,]))
      temp2[i,ind] = temp[i,ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
    }
  }
  
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

parse_text_input <- function(text_input){
  temp = matrix(unlist(strsplit(text_input, split = ",|\t|\n")), ncol = 5, byrow = T)
  colnames(temp) = c("gene", "aapos", "aaref", "aaalt", "raw_score")
  df_all = data.frame(temp)
  df_all$aapos = as.numeric(df_all$aapos)
  df_all$raw_score = as.numeric(df_all$raw_score)
  return(df_all)
}

check_DMS_input <- function(df){
  aa_list = unlist(strsplit("RHKDESTNQCGPAVILMFYW*", split = ""))
  msg = c()
  if (ncol(df) < 5){
    msg = c(msg, "Error: DMS data has less than 5 columns!")
    return(list(df, msg))
  }
  
  col_names = c("gene", "aapos", "aaref", "aaalt", "raw_score")
  if (!all(col_names %in% colnames(df))){
    missing_col = paste(col_names[!(col_names %in% colnames(df))], collapse = ",")
    msg = c(msg, paste0("Error: DMS data missing ", missing_col, " columns!"))
    return(list(df, msg))
  }
  
  if ("Z" %in% df$aaref | "Z" %in% df$aaalt){
    msg = c(msg, "Warning: amino acid Z detected! Converted to stop codon '*'.")
    df$aaref[df$aaref == "Z"] = "*"
    df$aaalt[df$aaalt == "Z"] = "*"
  }
  
  if (!all(df$aaref %in% aa_list) | !all(df$aaalt %in% aa_list)){
    msg = c(msg, "Warning: unknown amino acid detected! Rows removed.")
    df = df[df$aaref %in% aa_list,]
    df = df[df$aaalt %in% aa_list,]
  }
  
  if (class(df$aapos) != "numeric"){
    msg = c(msg, "Warning: non-numeric animo acid positions detected! Rows removed.")
    ind = grep(pattern = "\\D", x = df$aapos)
    df = df[-ind,]
    df$aapos = as.numeric(df$aapos)
  }
  
  if (class(df$raw_score) != "numeric"){
    msg = c(msg, "Warning: non-numeric DMS score detected! Rows removed.")
    ind = grep(pattern = "\\D", x = df$aapos)
    df = df[-ind,]
    df$aapos = as.numeric(df$aapos)
  }
  
  return(list(df, msg))
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

correct_direction <- function(df, df.funsum, include_LOF=F){
  gene = df$gene[1]
  df.sub_score = get_norm_score(df, pos_mean_method = "js", include_LOF = include_LOF, score_type = "raw_score", sd_norm = F)
  df.funsum_gene = get_FUNSUM(df.sub_score, avg_method = "js", include_LOF = include_LOF)
  res = cor.test(as.numeric(df.funsum_gene), as.numeric(df.funsum_all), method = "pearson")
  if (res$estimate < 0){
    msg = paste0("Info: Functional data for ", gene, " in opposite direction, score inverted.")
    return(list("inverted", msg))
  } else {
    return(list("same", ""))
  }
}

normalize_scoreset_lite <- function(df, lower_bound = 0.05, upper_bound = 0.95, direction = "same"){
  if (!"functional_class" %in% colnames(df)){
    df$functional_class = NA
    ind = which(df$aaref != df$aaalt & df$aaalt != '*')
    df$functional_class[ind] = "MIS"
    ind = which(df$aaref != "*" & df$aaalt == "*")
    df$functional_class[ind] = "LOF"
    ind = which(df$aaref == df$aaalt)
    df$functional_class[ind] = "SYN"
  }
  
  ind_mis = which(df$functional_class == "MIS")
  temp = df$raw_score
  
  if (direction == "inverted"){
    temp = -temp
  }
  mis_scores = temp[ind_mis]
  median_syn = quantile(mis_scores, lower_bound, na.rm = T)
  median_lof = quantile(mis_scores, upper_bound, na.rm = T)
  df$norm_raw_score = (temp-median_syn)/(median_lof-median_syn)
  return(df)
}



# Define UI for application
ui <- dashboardPage(
  dashboardHeader(title = "FUSE"), # Functional Substitution Estimation
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName = "intro", icon = icon("circle-info")),
      menuItem("Online tool", tabName = "tool", icon = icon("magnifying-glass-chart"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "intro",
              h2(HTML('<b>Fu</b>nctional <b>S</b>ubstitution <b>E</b>stimation (<b>FUSE</b>)')),
              span(version_number),
              br(),
              h4(HTML('<b>About</b>')),
              p(HTML('Read the full manuscript <a ref="https://www.medrxiv.org/content/10.1101/2023.01.06.23284280v1">here</a>.')),
              br(),
              h4(HTML('<b>Abstract</b>')),
              imageOutput("image1"),
              p(HTML('Deep mutational scanning enables the functional assessment of variants in high throughput. 
                     Although phenotypic measurements from screening assays are broadly concordant with clinical outcomes, systematic sources of statistical and experimental noise may affect the accuracy of individual variant estimates. 
                     We have developed a framework that leverages measurements collectively within a functional screening assay to improve the estimation of individual variant impacts. 
                     We estimate the mean functional effect per amino acid residue position, and make estimates for individual allelic variants with a pipeline called <b>FUSE</b> (<b>FU</b>nctional <b>S</b>ubstitution <b>E</b>stimation), which draws on data from 115 prior functional assays. 
                     FUSE enhances the correlation of functional effects from different assay platforms covering the same variants and increases classification accuracy of missense variants in ClinVar across 29 genes (AUC 0.83 to 0.90). 
                     In UK Biobank patients with a rare missense variant in <i>BRCA1</i>, <i>LDLR</i>, or <i>TP53</i>, FUSE improves the classification accuracy of associated phenotypes. 
                     Additionally, FUSE imputes predicted variant effects for substitutions not experimentally screened. 
                     This approach improves accuracy and broadens the utility of data generated from functional screening assays.')),
              br(),
              h4(HTML('<b>Disclaimer</b>')),
              p(HTML('The information presented is strictly for research use only, and should not be construed or used as a substitute for professional medical advice or diagnosis. We encourage you to share any information you find relevant with a trained medical professional, including results from this site, manuscript, or accompanying software tools and downloads. This application, software tools, and accompanying downloads are released under the Creative Commons non-commercial use license.')),
              p(HTML('<a href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="License: CC BY-NC-ND 4.0" src="https://licensebuttons.net/l/by-nc-nd/4.0/80x15.png"></a>')),
              br(),
              h4(HTML('<b>Contact</b>')),
              p(HTML('Please direct comments or questions to tyu7@bwh[dot]harvard[dot]edu.'))
              
      ),
      
      # Second tab content
      tabItem(tabName = "tool",
              fluidRow(
                column(width = 5, 
                       textAreaInput(inputId = "txt_in", label = "Input functional data", placeholder = example_txt, value = "", height='200px', width = '100%'),
                       actionButton(inputId = "btn_example", label = "Load example", width = '150px'),
                       actionButton(inputId = "btn_clear", label = "Clear", width = '150px'),
                       br(), br(),
                       fileInput(inputId = "file_in", label = "Upload DMS data (.csv/.tsv)", multiple = F, accept = c(".csv", ".tsv")),
                       hr(),
                       radioButtons(inputId = "radio_direction", label = "Direction of the funtional scores", 
                                    choices = list("Same" = "same", "Inverted" = "inverted", "Auto-detect" = "auto"), selected = "auto"),
                       hr(),
                       fileInput(inputId = "btn_dss", label = "Upload DSSP data (.dss)", multiple = F, accept = c(".dss")),
                       checkboxInput(inputId = "check_funsum2d", label = "Use 2nd structure specific FUNSUM", value = F),
                       hr(),
                       actionButton(inputId = "btn_submit", label = "Estimate functional scores", width = '250px'), br(),br(),
                       verbatimTextOutput(outputId = "msg_1"),
                ),
                column(width = 7,
                       tabsetPanel(
                         tabPanel(title = "Data output", br(), br(),
                                  DTOutput("tb_1")),
                         tabPanel(title = "Help", br(),
                                  h4(HTML('<b>How to use this tool</b>')),
                                  p(HTML('<ol>
                                            <li>There are two ways to read fucntional data - input through the textarea or upload a .tsv or .csv file. The input data should contains at least the following five columns:</li>
                                              <ul>
                                                <li><b>gene</b>  - Gene names or identifiers for independent assays</li>
                                                <li><b>aapos</b> - The amino acid position of the substitution</li>
                                                <li><b>aaref</b> - The reference amino acid in the substitution</li>
                                                <li><b>aaalt</b> - The alternative (mutated) amino acid in the substitution</li>
                                                <li><b>raw_score</b> - The assay measured functional impact of the substitution</li>
                                              </ul>
                                            <br>
                                            <li>Choose the direction of the functional scores, or allow FUSE to auto-detect the direction via correlation with the overall FUNSUM scores. By default, more detrimental variants have higher FUSE estiamted scores. </li>
                                            <br>
                                            <li>Choose whether or not to calculate FUSE estimated scores using 2nd structure-specific FUNSUM. To use this option, the input data must only contain dataset from one gene (or assay), and the corresponding DSSP output file must be uploaded. DSSP output files can be obtained by running DSSP using a .pdb structure file. By default, FUSE uses the overall FUNSUM. </li>
                                            <br>
                                            <li>Finally, click <b>Estimate functional score</b> and FUSE will output the estimated functional scores! </li>
                                          </ol>'))
                         )
                       )
                )
              )
      )
    )
  )
)




# Define server logic 
server <- function(input, output, session) {
  ## empty initial log
  textLog <- reactiveVal("")
  data_source <- reactiveVal(NULL)
  dssFile <- reactiveVal(NULL)
  
  observeEvent(input$btn_example, {
    example_str = paste(apply(df.example, MARGIN = 1, FUN = paste, collapse="\t"), collapse = "\n")
    example_str = gsub(pattern = " ", replacement = "", x = example_str)
    updateTextAreaInput(inputId = "txt_in", value = example_str)
    
  })
  
  observeEvent(input$txt_in, {
    data_source("txt_in")
    textLog(c(textLog(), "Info: data source changed to textarea!"))
  })
  
  observeEvent(input$file_in, {
    data_source("file_in")
    textLog(c(textLog(), "Info: Data source changed to file upload!"))
  })
  
  observeEvent(input$btn_clear, {
    updateTextAreaInput(inputId = "txt_in", value = "")
  })
  
  observe({
    if (input$check_funsum2d) {
      shinyjs::enable("btn_dss")  # Enable file input when checkbox is checked
    } else {
      shinyjs::disable("btn_dss")  # Disable (grey out) file input when checkbox is not checked
    }
  })
  
  observeEvent(input$btn_dss, {
    req(input$btn_dss)  # Ensure the file is uploaded
    textLog(c(textLog(), "Info: DSSP file is uploaded!"))
    dssFile(input$btn_dss$datapath)
  })
  
  observeEvent(input$btn_submit, {
    ## load table
    check_input = F
    if (is.null(data_source())){
      textLog(c(textLog(), "Error: DMS data not found!"))
    }
    else {
      if (data_source() == "file_in"){
        ext <- tools::file_ext(input$file_in$name)
        if (ext == "csv"){
          df_all = vroom(input$file_in$datapath, delim = ",")
          check_input = T
        } else if (ext == "tsv"){
          df_all = vroom(input$file_in$datapath, delim = "\t")
          check_input = T
        } else {
          textLog(c(textLog(), "Error: invalid file! Please upload a .csv or .tsv file."))
        }
      } else if (data_source() == "txt_in"){
        df_all = parse_text_input(text_input = input$txt_in)
        check_input = T
      }
      
      ## check DMS input requirements
      temp = check_DMS_input(df = df_all)
      df_all = temp[[1]]
      msg = temp[[2]]
      textLog(c(textLog(), msg))
      
      if (length(msg) > 0){
        if (grepl(pattern = "error", x = msg, ignore.case = T)){
          check_input = F
        }
      }
      
      if (check_input){
        gene_ids = unique(df_all$gene)
        textLog(c(textLog(), paste0("Info: Total genes/datasets from input: ", length(gene_ids))))
        
        ls.df = list()
        df_all.out = c()
        if (length(gene_ids) > 1){
          if (input$check_funsum2d){
            textLog(c(textLog(), paste0("Warning: 2nd-structure specific FUSE estimations are only avalible for single gene/assay dataset. Fall back to overall FUNSUM.")))
          } else {
            textLog(c(textLog(), paste0("Info: Estimating functional impact using overall FUNSUM.")))
          }
          
          for (gene_id in gene_ids){
            df = df_all %>% filter(gene == gene_id)
            
            ## check and correct DMS data direction
            if (input$radio_direction == "auto"){
              temp = correct_direction(df = df, df.funsum = df.funsum_all, include_LOF = include_LOF) 
              direction = temp[[1]]
              msg = temp[[2]]
              textLog(c(textLog(), msg))
            } else {
              direction = input$radio_direction
            }
            
            ## normalize DMS data
            df = normalize_scoreset_lite(df = df, lower_bound = 0.1, upper_bound = 0.9, direction = direction)
            
            ## de-noise using overall FUNSUM
            df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF, show_func_class = T)
            ls.df[[gene_id]] = df.out
            df_all.out = rbind(df_all.out, df.out)
          }
          
          df_all.out = df_all.out %>% select(gene, aapos, aaref, aaalt, functional_class, raw_score, norm_raw_score, pos_score, sub_score, FUSE_score = final_score)
          
        } else {
          df = df_all
          gene_id = df$gene[1]
          
          ## check and correct DMS data direction
          if (input$radio_direction == "auto"){
            print("before correct")
            temp = correct_direction(df = df, df.funsum = df.funsum_all, include_LOF = include_LOF) 
            print("after correct")
            direction = temp[[1]]
            msg = temp[[2]]
            textLog(c(textLog(), msg))
          } else {
            direction = input$radio_direction
          }
          
          ## normalize DMS data
          df = normalize_scoreset_lite(df = df, lower_bound = 0.1, upper_bound = 0.9, direction = direction)
          
          if (input$check_funsum2d){
            textLog(c(textLog(), paste0("Info: Estimating functional impact using 2nd-structure specific FUNSUM.")))
            
            req(input$btn_dss)  # Ensure the file is uploaded
            dir.create("dssp", showWarnings = FALSE)
            tempFile <- input$btn_dss$datapath
            df.out = de_noise_ss_1gene(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, ls.funsum_ss = ls.funsum_ss, 
                                       dss_path = dssFile(), include_LOF = include_LOF, show_func_class = T)
            
            ls.df[[gene_id]] = df.out
            df_all.out = rbind(df_all.out, df.out)
            df_all.out = df_all.out %>% select(gene, aapos, aaref, aaalt, functional_class, DSSP_ss = ss, DSSP_acc = acc, raw_score, norm_raw_score, pos_score, sub_score = sub_score_ss, FUSE_score = final_score_ss)
          } else {
            textLog(c(textLog(), paste0("Info: Estimating functional impact using overall FUNSUM.")))
            
            df.out = de_noise(df = df, pos_mean_method = "js", df.funsum = df.funsum_all, include_LOF = include_LOF, show_func_class = T)
            ls.df[[gene_id]] = df.out
            df_all.out = rbind(df_all.out, df.out)
            df_all.out = df_all.out %>% select(gene, aapos, aaref, aaalt, functional_class, raw_score, norm_raw_score, pos_score, sub_score, FUSE_score = final_score)
          }
        }
        
        ## render de-noised DMS data table
        output$tb_1 <- {renderDT(server = F,
                                 datatable(df_all.out, selection = "single", style = 'default', rownames=F, 
                                           extensions = c('Buttons'),
                                           options = list(pageLength=20, dom = 'Bfrtip', searchHighlight = TRUE, scrollX = TRUE))
        )}
      }
    }
  })
  
  observe({
    msg = textLog()
    print(msg)
    if(length(msg)<10){1
      msg = c(msg, rep("", 10-length(msg)))
    } else if (length(msg)>10) {
      msg = tail(msg, n = 10)
    }
    output$msg_1 <- renderText({paste(msg, collapse = "\n")})
  })
  
  # Send a pre-rendered image, and don't delete the image after sending it
  output$image1 <- renderImage({
    filename <- normalizePath(file.path('www/images/graphical_abstract_v2.png'))
    list(src = filename, contentType = "png", alt = 'Graphical abstract', height = '400px')
  }, deleteFile = FALSE)
}

# Run the application 
shinyApp(ui = ui, server = server)
