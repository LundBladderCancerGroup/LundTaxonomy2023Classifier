#' @title Predict Grade
#'
#' @description Function for predicting grade using random forest and the bundled classifiers.
#'
#' @details Internal function called by [LundTax2023Classifier::lundtax_calc_sigscore()]. 
#' Not meant for out of package use. This function is internally calling 
#' [multiclassPairs::predict_RF()] to run a random forest prediction using one of the two bundled 
#' classifiers for grade prediction.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param grade_predictor Required parameter, the predictor needed for grading.
#' Should be one of the following bundled classifiers (classifier_GRADE3 or classifier_HG).
#' @param gene_id Specify the type of gene identifier used in `this_data`. 
#' Accepted values are; hgnc_symbol (default) or ensembl_gene_id.
#' @param impute From [multiclassPairs::predict_RF()]. Boolean. To determine if 
#' missed genes and NA values should be imputed or not. The non missed rules will 
#' be used to detemine the closest samples in the training binary matrix (i.e. 
#' which is stored in the classifier object). For each sample, the mode value for 
#' nearest samples in the training data will be assigned to the missed rules. Default is FALSE.
#' @param impute_reject From [multiclassPairs::predict_RF()]. A number between 0 
#' and 1 indicating the threshold of the missed rules in the sample. Based on this 
#' threshold the sample will be rejected (i.e. skipped if higher than the 
#' impute_reject threshold) and the missed rules will not be imputed in this sample. 
#' Default is 0.67. NOTE, The results object will not have any results for this sample.
#' @param impute_kNN From [multiclassPairs::predict_RF()]. Integer determines 
#' the number of the nearest samples in the training data to be used in the 
#' imputation. Default is 5. It is not recommended to use large number (i.e. >10).
#' @param verbose A logical value indicating whether processing messages will be 
#' printed or not. Default is TRUE.
#'
#' @return A data frame with grade results.
#'
#' @importFrom multiclassPairs predict_RF
#' 
int_predict_grade = function(this_data = NULL, 
                             grade_predictor = NULL,
                             gene_id = "hgnc_symbol",
                             impute = FALSE, 
                             impute_reject = 0.67, 
                             impute_kNN = 5, 
                             verbose = TRUE){
  
  #check the incoming data
  if(!class(this_data)[1] %in% c("data.frame","matrix")){
    stop("Data must be in dataframe or matrix format...")
  }
  
  #check the predictor class
  if(class(grade_predictor)[1] != "rule_based_RandomForest"){
    stop("Classifier must be a multiclassPairs::rule_based_RandomForest object...")
  }
  
  #update gene names in incoming data from hgnc_symbols to ensembl_gene_id, which is required by the rpedictor.
  if(gene_id == "hgnc_symbol"){
    
    #convert the rownames to first column
    mutated_data = dplyr::as_tibble(this_data, 
                                    rownames = "hgnc_symbol")
    
    #left join with gene list to get ensembl IDs
    mutated_data = dplyr::left_join(gene_list, 
                                    mutated_data, 
                                    by = "hgnc_symbol")
    
    #remove duplicated rows
    mutated_data = dplyr::filter(mutated_data, 
                                 duplicated(ensembl_gene_id) == FALSE)
    
    #convert the first column back to rownames
    row.names(mutated_data) = unique(mutated_data$ensembl_gene_id)
    mutated_data[1:2] <- NULL
    
    #convert back to expected name
    this_data = mutated_data
  }
  
  #run predict_RF from multiclassPairs
  grade_results = multiclassPairs::predict_RF(classifier = grade_predictor, 
                                              Data = this_data, 
                                              impute = impute, 
                                              impute_reject = impute_reject, 
                                              impute_kNN = impute_kNN, 
                                              verbose = verbose)
  
  #return results
  return(grade_results)
}
