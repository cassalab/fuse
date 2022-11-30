library(shiny)
library(vroom)
library(DT)

## global variables
df.example = read.table("./example_DMS.txt", header = T, sep = "\t")
ls.funsum_all = readRDS("./funsum.rds")
df.NSFP_all = read.table("./dfNFSP_031022.txt", header = T, sep = "\t")
df.NSFP_all$aaref[which(df.NSFP_all$aaref == "X")] = "*"
df.NSFP_all$aaalt[which(df.NSFP_all$aaalt == "X")] = "*"
example_txt = paste(apply(df.example[1:5,], MARGIN = 1, FUN = paste, collapse="\t"), collapse = "\n")
score_types = c("norm_raw_score", "pos_score", "sub_score", "final_score")
score_names = c("Normalized DMS score", "FUNSUM positional score", "FUNSUM substitution score", "FUNSUM final score")
names(score_names) = score_types


## functions
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
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

funsum_to_subTable <- function(df.funsum){
  df.sub_tb = c()
  for (i in 1:nrow(df.funsum)){
    temp = data.frame(aa_pair = paste0(rownames(df.funsum)[i], colnames(df.funsum)), score = as.vector(df.funsum[i,]))
    df.sub_tb = rbind(df.sub_tb, temp)
  }
  
  return(df.sub_tb)
}

parse_text_input <- function(text_input){
  temp = matrix(unlist(strsplit(text_input, split = ",|\t|\n")), ncol = 5, byrow = T)
  colnames(temp) = c("gene", "aapos", "aaref", "aaalt", "raw_score")
  df_all = data.frame(temp)
  df_all$aapos = as.numeric(df_all$aapos)
  df_all$raw_score = as.numeric(df_all$raw_score)
  return(df_all)
}

de_noise <- function(df, pos_mean_method, df.funsum, sd_norm = T){
  ## need to add test for direction
  
  gene_id = unique(df$gene)
  aa_list = colnames(df.funsum)
  df.sub_tb = funsum_to_subTable(df.funsum) #convert FUNSUM to tabular format
  df$aa_str = paste0(df$aaref, df$aapos, df$aaalt)
  df$aa_pair = paste0(df$aaref, df$aaalt)
  
  ## normalize the data by sd
  if (sd_norm){
    df$norm_raw_score = (df$raw_score-median(df$raw_score, na.rm = T))/sd(df$raw_score, na.rm = T)
  }
  
  ## calculate positional component
  df.pos_score = data.frame(aapos = unique(df$aapos), aaref = NA, pos_mean = NA)
  temp = matrix(NA, nrow = nrow(df.pos_score), ncol = length(aa_list))
  colnames(temp) <- aa_list
  
  df.out = c()
  for (i in 1:nrow(df.pos_score)){
    ind = which(df$aapos == df.pos_score$aapos[i])
    df.pos_score$aaref[i] = unique(df$aaref[ind])
    ind2 = df$aaalt[ind] %in% aa_list
    temp[i, df$aaalt[ind[ind2]]] = df$norm_raw_score[ind[ind2]]
    
    # calculate pos_mean by funsum method
    if (pos_mean_method == "funsum"){
      ind = which(!is.na(temp[i,]))
      temp2 = temp[i,ind] - df.funsum[df.pos_score$aaref[i], aa_list[ind]]
      df.pos_score$pos_mean[i] = mean(as.numeric(temp2))
    }
    
    ## construct a new df with all possible substitutions
    df.temp = data.frame(aapos = df.pos_score$aapos[i], aaref = df.pos_score$aaref[i], aaalt = aa_list)
    df.out = rbind(df.out, df.temp)
  }
  
  # calculate pos_mean by other methods
  if (pos_mean_method == "mean"){
    df.pos_score$pos_mean = rowMeans(temp, na.rm = T)
  } else if (pos_mean_method == "median"){
    df.pos_score$pos_mean = apply(temp, 1, FUN = median, na.rm = T)
  } else if (pos_mean_method == "js"){
    df.pos_score$pos_mean = get_js(t(temp))
  } 
  
  # assign norm_score, pos_score, sub_score to df.out
  df.out$aa_str = paste0(df.out$aaref, df.out$aapos, df.out$aaalt)
  df.out$aa_pair = paste0(df.out$aaref, df.out$aaalt)
  ind = match(df.out$aa_str, table = df$aa_str)
  df.out$raw_score = df$raw_score[ind]
  df.out$norm_raw_score = df$norm_raw_score[ind]
  ind = match(df.out$aapos, table = df.pos_score$aapos)
  df.out$pos_score = df.pos_score$pos_mean[ind]
  ind = match(df.out$aa_pair, table = df.sub_tb$aa_pair)
  df.out$sub_score = df.sub_tb$score[ind]
  df.out$final_score = df.out$pos_score + df.out$sub_score
  df.out = cbind(gene = gene_id, df.out)
  
  return(df.out)
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

check_direction <- function(df){
  msg = c()
  gene_id = unique(df$gene)
  ind = which((df$aaref != "P" & df$aaalt == "P") | (df$aaalt == "*"))
  if (length(ind) > 0){
    if (mean(df$raw_score) > mean(df$raw_score[ind])){
      df$raw_score = - df$raw_score
      msg = c(msg, paste0("Info: ", gene_id, " DMS data in opposite direction. Score inverted."))
    }
  } else {
    msg = c(msg, paste0("Warning: ", gene_id, " DMS data direction cannot be determined!"))
  }
  
  return(list(df, msg))
}

draw_validation_plot <- function(df.NSFP, score_types){
  par(mar = c(1.1, 3.1, 2.1, 1.1), oma = c(3.1, 1.1, 1.1, 1.1), mfcol=c(4,1))
  for (score_type in score_types){
    p1 = NA
    
    ind = which(df.NSFP$cs_simple == "b")
    score_g1 = df.NSFP[[score_type]][ind]
    score_g1 = score_g1[which(!is.na(score_g1))]
    score_g1.density = density(score_g1, bw = 0.5)
    
    ind = which(df.NSFP$cs_simple == "p")
    score_g2 = df.NSFP[[score_type]][ind]
    score_g2 = score_g2[which(!is.na(score_g2))]
    score_g2.density = density(score_g2, bw = 0.5)
    
    res = ks.test(x = score_g1, y = score_g2) # KS-test
    p1 = format(res$p.value, digits = 3)
    
    ymax = max(c(score_g1.density$y, score_g2.density$y))
    xmax = max(abs(c(score_g1.density$x, score_g2.density$x)))
    plot(score_g1.density, col = "blue", ylim=c(0,1), xlim=c(-5,5), las=1, main = "", xaxs="i", yaxs="i", ylab="", xlab="") # g1 is blue
    
    text(x = par("usr")[1]+(par("usr")[2]-par("usr")[1])*0.8, y=par("usr")[3]+(par("usr")[4]-par("usr")[3])*0.9, paste0("KS-test p-val = ", p1))
    lines(score_g2.density, col = "red") # g2 is red
    legend("topleft", legend = c(paste0("B/LB (n=", length(score_g1), ")"), paste0("P/LP (n=", length(score_g2), ")")), col = c("blue", "red"), lty = 1)
  }
}



# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("DMS data de-noiser"),
  
  fluidRow(
    column(width = 5, 
           textAreaInput(inputId = "txt_1", label = "Input DMS data", placeholder = example_txt, value = "", height='200px', width = '100%'),
           actionButton(inputId = "btn_1", label = "Load example", width = '150px'),
           actionButton(inputId = "btn_2", label = "Clear", width = '150px'),
           br(), br(),
           fileInput(inputId = "file_1", label = "Upload DMS data (.csv/.tsv)", multiple = F, accept = c(".csv", ".tsv")),
           hr(),
           radioButtons(inputId = "radio_1", label = "Method to calculate positional component", choices = list("FUNSUM" = "funsum", "James-Stein" = "js", "Median" = "median", "Mean" = "mean")),
           radioButtons(inputId = "radio_2", label = "Include stop-codon substitutions?", choices = list("Yes"=T, "No"=F)),
           radioButtons(inputId = "radio_3", label = "Validation data source", choices = c("ClinVar", "Other")),
           hr(),
           actionButton(inputId = "btn_3", label = "De-noise DMS data", width = '150px'), br(),br(),
           verbatimTextOutput(outputId = "msg_1"),
    ),
    column(width = 7,
           tabsetPanel(
             tabPanel(title = "Data output", br(), br(),
                      DTOutput("tb_1")),
             tabPanel(title = "Visualization", br(), br(),
                      selectInput(inputId = "select_1", label = "Gene to visualize", choices = c("None")),
                      checkboxInput(inputId = "checkbox_1", label = "Use only the variants avaliable in original DMS data", value = F),
                      plotOutput(outputId = "plot_1"))
           )
    )
  ),
)

# Define server logic 
server <- function(input, output, session) {
  ## empty initial log
  textLog <- reactiveVal("")
  data_source <- reactiveVal(NULL)
  ls.NSFP_reac <- reactiveValues()
  
  observeEvent(input$btn_1, {
    example_str = paste(apply(df.example, MARGIN = 1, FUN = paste, collapse="\t"), collapse = "\n")
    example_str = gsub(pattern = " ", replacement = "", x = example_str)
    updateTextAreaInput(inputId = "txt_1", value = example_str)
    
  })
  
  observeEvent(input$txt_1, {
    data_source("txt_1")
    textLog(c(textLog(), "Info: data source changed to textarea!"))
  })
  
  observeEvent(input$file_1, {
    data_source("file_1")
    textLog(c(textLog(), "Info: Data source changed to file upload!"))
  })
  
  observeEvent(input$btn_2, {
    updateTextAreaInput(inputId = "txt_1", value = "")
  })
  
  # observeEvent(input$radio_2, {
  #   textLog(paste(textLog(), paste0("\nRadio_2: ", input$radio_2)))
  # })
  
  observeEvent(input$btn_3, {
    ## choose FUNSUM matrix
    if(input$radio_2){
      df.funsum_all = ls.funsum_all$with_stop
    } else {
      df.funsum_all = ls.funsum_all$no_stop
    }
    
    ## load table
    check_input = F
    if (is.null(data_source())){
      textLog(c(textLog(), "Error: DMS data not found!"))
    }
    else {
      if (data_source() == "file_1"){
        ext <- tools::file_ext(input$file_1$name)
        if (ext == "csv"){
          df_all = vroom(input$file_1$datapath, delim = ",")
          check_input = T
        } else if (ext == "tsv"){
          df_all = vroom(input$file_1$datapath, delim = "\t")
          check_input = T
        } else {
          textLog(c(textLog(), "Error: invalid file! Please upload a .csv or .tsv file."))
        }
      } else if (data_source() == "txt_1"){
        df_all = parse_text_input(text_input = input$txt_1)
        check_input = T
      }
      
      temp = check_DMS_input(df = df_all)
      df_all = temp[[1]]
      msg = temp[[2]]
      textLog(c(textLog(), msg))
      
      print(paste0("df_all nrow: ", nrow(df_all)))
      
      if (length(msg) > 0){
        if (grepl(pattern = "error", x = msg, ignore.case = T)){
          check_input = F
        }
      }
      
      ## remove stop-codon substitutions if choose to not include stop-codon
      if (input$radio_2 == F){
        ind = which(df_all$aaref != "*" & df_all$aaalt != "*")
        df_all = df_all[ind,]
      }
      
      if (check_input){
        gene_ids = unique(df_all$gene)
        textLog(c(textLog(), paste0("Info: number of genes: ", length(gene_ids))))
        
        ls.df = list()
        df_all.out = c()
        df.NSFP_all_out = c()
        for (gene_id in gene_ids){
          ## de-noise DMS data for each gene
          ind = which(df_all$gene == gene_id)
          df = df_all[ind,]
          
          temp = check_direction(df) # check and correct DMS data direction
          df = temp[[1]]
          msg = temp[[2]]
          textLog(c(textLog(), msg))
          
          df.out = de_noise(df, pos_mean_method = input$radio_1, df.funsum = df.funsum_all)
          ls.df[[gene_id]] = df.out
          df_all.out = rbind(df_all.out, df.out)
          
          ## load df.NSFP for each gene
          ind = which(df.NSFP_all$genename == gene_id)
          if (length(ind)){
            df.NSFP = df.NSFP_all[ind,]
            df.NSFP$aa_str = paste0(df.NSFP$aaref, df.NSFP$aapos, df.NSFP$aaalt)
            ind = match(df.NSFP$aa_str, df.out$aa_str)
            df.NSFP$norm_raw_score = df.out$norm_raw_score[ind]
            df.NSFP$pos_score = df.out$pos_score[ind]
            df.NSFP$sub_score = df.out$sub_score[ind]
            df.NSFP$final_score = df.out$final_score[ind]
            ls.NSFP_reac[[gene_id]] = df.NSFP
            df.NSFP_all_out = rbind(df.NSFP_all_out, df.NSFP)
          } else {
            textLog(c(textLog(), paste0("Warning: no ClinVar record found for ", gene_id, "!")))
          }
        }
        ls.NSFP_reac[["All"]] = df.NSFP_all_out
        
        # update selectInput values
        updateSelectInput(session = session, inputId = "select_1", choices = names(ls.NSFP_reac), selected = "All")
        
        ## render de-noised DMS data table
        ind = which(df_all.out$aaref == df_all.out$aaalt)
        df_all.out = df_all.out[-ind, c("gene", "aapos", "aaref", "aaalt", "raw_score", "norm_raw_score", "pos_score", "sub_score", "final_score")]
        output$tb_1 <- {renderDT(
          datatable(df_all.out, selection = "single", style = 'default', rownames=F,
                    extensions = c('Buttons'),
                    options = list(pageLength=10, dom = 'Bfrtip', searchHighlight = TRUE))
        )}
      }
    }
  })
  
  # observeEvent(input$select_1, {
  #   if (input$select_1 != "None"){
  #     df.NSFP = ls.NSFP_reac[[input$select_1]]
  #     
  #     ## render validation plot
  #     output$plot_1 <- renderPlot({
  #       draw_validation_plot(df.NSFP, score_types)
  #     }, height = 800, width = 600, res = 120)
  #   }
  # })
  
  observe({
    if (input$select_1 != "None"){
      df.NSFP = ls.NSFP_reac[[input$select_1]]
      if (input$checkbox_1){
        ind = which(!is.na(df.NSFP$norm_raw_score))
        df.NSFP = df.NSFP[ind,]
        print("checkbox on!")
      }
      
      ## render validation plot
      output$plot_1 <- renderPlot({
        draw_validation_plot(df.NSFP, score_types)
        print(paste0("df.NSFP nrow: ", nrow(df.NSFP)))
        print("plot updated!")
      }, height = 800, width = 600, res = 120)
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
}

# Run the application 
shinyApp(ui = ui, server = server)
