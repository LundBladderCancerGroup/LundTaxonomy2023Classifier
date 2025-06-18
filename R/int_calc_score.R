#' @title Calculate Scores.
#'
#' @description Calculate immune and infiltration scores.
#'
#' @details Internal function called by [LundTax2023Classifier::lundtax_calc_sigscore()]. 
#' Not meant for out of package use. Takes a data frame of matrix with expression values and 
#' calculates scores based on gene expression.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param variable Required parameter. Input should be one of the following; immune, 
#' score141up, proliferation, or progression.
#' @param gene_id Specify the type of gene identifier used in `this_data`. 
#' Accepted values are; hgnc_symbol (default) or ensembl_gene_id.
#' @param adjust Boolean parameter. If TRUE, the function will proceed with 
#' adjusting the scores based on stable genes. If FALSE (default), no adjustment 
#' will be made and the original score values will be retained. 
#' @param adj_factor Only applicable if adjust is set to TRUE. Allows users to 
#' apply a proportional adjustment to the normalized scores, enabling finer 
#' control over the final output values. After dividing each score by the mean 
#' expression of stable genes, the result is multiplied by this factor. 
#' Default is 5.1431
#' @param verbose A logical value indicating whether processing messages will be 
#' printed or not. Default is TRUE.
#' 
#' @return A list with two objects. 1, A data frame with scores for the selected variable. 
#' 2, A data frame indicating what genes from the incoming data are missing, based on the expected 
#' genes for signature calculations.
#' 
#' @import dplyr
#'
int_calc_score = function(this_data = NULL,
                          variable = NULL,
                          gene_id = "hgnc_symbol",
                          adjust = TRUE,
                          adj_factor = 5.1431,
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
  
  #get scores for selected varaible
  if(variable == "immune"){
    if(verbose){
      message("Calculating immune scores...")
    }
    
    these_signatures = dplyr::select(signatures$immune, signature, !!as.symbol(gene_id))
  }else if(variable == "score141up"){
    if(verbose){
      message("Calculating 141 scores...") 
    }
    
    these_signatures = filter(signatures$signatures_plot, signature %in% c("Stromal141_UP", "Immune141_UP")) %>% 
      select(signature, !!as.symbol(gene_id))
  }

  #subset expression data to relevant prostate genes
  genes_int = this_data %>% 
    dplyr::filter(rownames(this_data) %in% unique(these_signatures[,gene_id]))
  
  #get genes that are not in the incoming data
  diff_genes = setdiff(these_signatures[,gene_id], rownames(genes_int))
  
  #create data frame with missing gene information
  missing_genes = data.frame(genes = as.character(), 
                             signature = as.character(), 
                             process = as.character())
  
  #append missing genes information to the data frame
  missing_genes = add_row(missing_genes, 
                          genes = diff_genes, 
                          signature = variable, 
                          process = "NA")
  
  #notify the user what genes are missing
  if(length(diff_genes > 0)){
    if(verbose){
      message(paste0(length(diff_genes), " out of ", length(unique(these_signatures[,gene_id]))," genes are missing from the data..."))
      print(diff_genes) 
    }
  }
  
  #create results object
  score_results <- as.data.frame(matrix(nrow = ncol(this_data),
                                        ncol = length(unique(these_signatures$signature)),
                                        dimnames = list(colnames(this_data), unique(these_signatures$signature))))
  
  for(i in unique(these_signatures$signature)){
    genes_signature_int = intersect(rownames(this_data),these_signatures[these_signatures$signature == i, gene_id])
    res = unlist(lapply(1:ncol(this_data),function(x) {mean(this_data[genes_signature_int,x])}))
    score_results[,i] = res
  }
  
  if(adjust){
    StableGenes <- signatures$stable_genes
    stable_genes_int <- intersect(rownames(this_data),StableGenes[,gene_id])
    score_results <- do.call("rbind",lapply(1:nrow(score_results),function(x){
      (score_results[x,]/mean(this_data[stable_genes_int,x]))*adj_factor
    }))
  }
  
  #calculate proportions
  if(variable == "immune"){
    immune__prop = int_calc_immune_proportions(immune_results = score_results)
    score_results = cbind(score_results, immune__prop)
  }
  
  return(list(sig_score = score_results, na_genes = missing_genes))
}

