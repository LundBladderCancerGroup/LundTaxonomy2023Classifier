#' @title Calculate Scores.
#'
#' @description Calculate proliferation and progression scores.
#'
#' @details Internal function called by [LundTax2023Classifier::lundtax_calc_sigscore()]. 
#' Not meant for out of package use. Takes a data frame of matrix with expression values and 
#' calculates scores based on gene expression.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param variable Required parameter. Input should be one of the following; 
#' proliferation, or progression.
#' @param gene_id Specify the type of gene identifier used in `this_data`. 
#' Accepted values are; hgnc_symbol (default) or ensembl_gene_id.
#' @param verbose A logical value indicating whether processing messages will be 
#' printed or not. Default is TRUE.
#' 
#' @return A list with two objects. 1, A data frame with scores for the selected variable. 
#' 2, A data frame indicating what genes from the incoming data are missing, based on the expected 
#' genes for signature calculations.
#' 
#' @import dplyr
#'
int_ratio_score = function(this_data = NULL,
                           variable = NULL,
                           gene_id = "hgnc_symbol",
                           verbose = TRUE){
  
  #check the incoming data
  if(!class(this_data)[1] %in% c("data.frame","matrix")){
    stop("Data must be in dataframe or matrix format...")
  }
  
  #check valid variables
  if(!variable %in% c("immune", "score141up", "proliferation", "progression")){
    stop("Variable must be one of the following: immune, score141up...")
  }
  
  #check gene format
  if(!gene_id %in% c("hgnc_symbol", "ensembl_gene_id")){
    stop("gene_id must be one of the following: 'hgnc_symbol' or 'ensembl_gene_id'")
  }

  if(variable == "proliferation"){
    if(verbose){
      message("Calculating proliferation scores...") 
    }
    
    #get all genes with signature LateCellCycle
    up_genes = filter(signatures$proliferation, signature == "LateCellCycle") %>% 
      select(signature, !!as.symbol(gene_id))
    
    #get all genes with signature EarlyCellCycle
    down_genes = filter(signatures$proliferation, signature == "EarlyCellCycle") %>% 
      select(signature, !!as.symbol(gene_id))
  }else if(variable == "progression"){
    if(verbose){
      message("Calculating progression scores...") 
    }
    
    #get all genes with signature LateCellCycle
    up_genes = filter(signatures$progression, direction_in_prog == "Up") %>% 
      select(direction_in_prog, !!as.symbol(gene_id))
    
    #get all genes with signature EarlyCellCycle
    down_genes = filter(signatures$progression, direction_in_prog == "Down") %>% 
      select(direction_in_prog, !!as.symbol(gene_id))
  }
    
  #intersect with incoming data
  int_up_genes = this_data %>% 
    dplyr::filter(rownames(this_data) %in% unique(up_genes[,gene_id]))
  
  int_down_genes = this_data %>% 
    dplyr::filter(rownames(this_data) %in% unique(down_genes[,gene_id]))
  
  #get genes that are not in the incoming data
  diff_genes_up = setdiff(up_genes[,gene_id], rownames(int_up_genes))
  diff_genes_down = setdiff(down_genes[,gene_id], rownames(int_down_genes))
  
  #create data frame with missing gene information
  missing_genes = data.frame(genes = as.character(), 
                             signature = as.character(), 
                             process = as.character())
  
  #append missing genes information to the data frame
  missing_genes = add_row(missing_genes, 
                          genes = diff_genes_up, 
                          signature = variable, 
                          process = "up")
  
  missing_genes = add_row(missing_genes, 
                          genes = diff_genes_down, 
                          signature = variable, 
                          process = "down")
  
  #notify the user what genes are missing
  if(length(diff_genes_up > 0)){
    if(verbose){
      message(paste0(length(diff_genes_up), " out of ", length(unique(up_genes[,gene_id]))," genes are missing from the data..."))
      print(diff_genes_up) 
    }
  }
  
  if(length(diff_genes_down > 0)){
    if(verbose){
      message(paste0(length(diff_genes_down), " out of ", length(unique(down_genes[,gene_id]))," genes are missing from the data..."))
      print(diff_genes_down)
    }
  }
    
  #ratio of ranks
  rank_data = apply(this_data, 2, rank)
  rank_data = rank_data/nrow(rank_data)
  
  #get vector with relevant genes
  up_vector = rownames(int_up_genes)
  down_vector = rownames(int_down_genes)
    
  median_up_genes = apply(rank_data[up_vector,,drop=F],2,median)
  median_down_genes = apply(rank_data[down_vector,,drop=F],2,median)
    
  median_up_down = median_up_genes/median_down_genes
    
  r_score = data.frame(Score = median_up_down, 
                       row.names = colnames(this_data))
  
  return(list(sig_score = r_score, na_genes = missing_genes))
}
