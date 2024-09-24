#' @title Lund Taxonomy Predict Subtype
#'
#' @description Predict Lund Taxonomy subtypes based on rule-based Random Forest classifiers.
#'
#' @details This function uses 2 classifiers to classify the samples: 5-class 
#' classifier first  classifies samples into Uro, GU, BaSq, Mes or ScNE.
#' Samples classified as Uro receive a second classification as UroA, B or C by 
#' the second classifier.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param gene_id Specify the type of gene identifier used in `this_data`. Accepted values are; 
#' hgnc_symbol (default) or ensembl_gene_id.
#' @param threshold_progression Threshold to flag a sample as high risk of progression, default is 0.58.
#' @param log_transform Boolean parameter. If TRUE, the function log transforms the incoming expression
#' values. Default is FALSE.
#' @param adjust Boolean parameter. If TRUE, the function will proceed with adjusting the scores based
#' on stable genes. If FALSE (default), no adjustment will be made and the original score values will be retained. 
#' @param adj_factor Only applicable if adjust is set to TRUE. Allows users to apply a proportional 
#' adjustment to the normalized scores, enabling finer control over the final output values. 
#' After dividing each score by the mean expression of stable genes, the result is multiplied by this factor. 
#' @param impute From [multiclassPairs::predict_RF()]. Boolean. To determine if missed genes and NA 
#' values should be imputed or not. The non missed rules will be used to determine the closest samples
#' in the training binary matrix (i.e. which is stored in the classifier object). For each sample, 
#' the mode value for nearest samples in the training data will be assigned to the missed rules. Default is FALSE.
#' @param impute_kNN From [multiclassPairs::predict_RF()]. Integer determines the number of the nearest
#' samples in the training data to be used in the imputation. Default is 5. It is not recommended to
#'  use large number (i.e. >10).
#' @param impute_reject From [multiclassPairs::predict_RF()]. A number between 0 and 1 indicating the 
#' threshold of the missed rules in the sample. Based on this threshold the sample will be rejected 
#' (i.e. skipped if higher than the impute_reject threshold) and the missed rules will not be imputed 
#' in this sample. Default is 0.67. NOTE, The results object will not have any results for this sample.
#' @param verbose A logical value indicating whether processing messages will be 
#' printed or not. Default is TRUE.
#' @param subtype_only Boolean parameter. Set to TRUE to return subtypes and nothing else. Default is FALSE.
#' @param include_data Boolean parameter. Set to TRUE to include data in output, FALSE is default.
#' @param include_pred_scores Boolean parameter. Set to TRUE (default) to include 
#' prediction scores for each sample and class in output.
#'
#' @return Returns a list object including: Data (optional, not included by default), 
#' Prediction scores for all classes (optional, included by default), Predicted LundTax 
#' class for 7-class system, Predicted LundTax class for 5-class system
#'
#' @import dplyr
#' @importFrom stats setNames median na.omit quantile
#' @importFrom multiclassPairs predict_RF
#'
#' @export
#' 
#' @examples
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                         impute = TRUE)
#'
lundtax_predict_sub = function(this_data = NULL,
                               gene_id = "hgnc_symbol",
                               threshold_progression = 0.58,
                               log_transform = TRUE,
                               adjust = FALSE,
                               adj_factor = 5.1431,
                               impute = FALSE, 
                               impute_reject = 0.67, 
                               impute_kNN = 5,
                               subtype_only = FALSE,
                               include_data = FALSE,
                               include_pred_scores = TRUE,
                               verbose = TRUE){
  #checks
  #check the incoming data
  if(!class(this_data)[1] %in% c("data.frame","matrix")){
    stop("Data must be in dataframe or matrix format...")
  }
  
  #check for duplicated samples
  if(ncol(this_data) != length(unique(colnames(this_data)))){
    stop("Sample names (column names) should not be duplicated")
  }
  
  #check gene format
  if(!gene_id %in% c("hgnc_symbol", "ensembl_gene_id")){
    stop("gene_id must be one of the following: 'hgnc_symbol' or 'ensembl_gene_id'")
  }
  
  #log transform
  if(log_transform){
    this_data = log2(this_data + 1)
  }
  
  #removed gene ID conversions tep, might be redundant since it's taken place in other internal functions

  #initiate results object
  results_suburo <- list(data = this_data,
                         subtype_scores = NULL,
                         predictions_7classes = NULL,
                         predictions_5classes = NULL,
                         scores = NULL)
  
  #predict 5 class
  prediction = predict_RF(classifier = classifier_lundtax_5c,
                          Data = this_data,
                          impute = impute, 
                          impute_reject = impute_reject, 
                          impute_kNN = impute_kNN, 
                          verbose = verbose)
  
  #reorder scores
  #set correct order
  pred_order = prediction$predictions[,c("Uro","GU","BaSq","Mes","ScNE"), drop = FALSE]
  
  #apply order
  prediction$predictions = pred_order
  prediction$predictions_classes = colnames(pred_order)[max.col(replace(pred_order, is.na(pred_order), -Inf), ties.method = "first")]
  
  #get uro samples
  if("Uro" %in% prediction$predictions_classes){
    
    #subset poredictions based on Uro classification
    pred_uro = this_data[,which(prediction$predictions_classes == "Uro"), drop = FALSE]
    pred_nouro = this_data[,which(prediction$predictions_classes != "Uro"), drop = FALSE]
    
    #predict uro samples
    prediction_suburo = predict_RF(classifier = classifier_lundtax_7c,
                                   Data = pred_uro,
                                   impute = impute, 
                                   impute_reject = impute_reject, 
                                   impute_kNN = impute_kNN, 
                                   verbose = verbose)
    
    #get column names for both data subsets
    names_uro = colnames(pred_uro)
    names_all = colnames(this_data)
    
    #merge score matrix
    score_matrix = int_merge_suburo_matrix(score_matrix1 = prediction_suburo$predictions,
                                           score_matrix2 = prediction$predictions,
                                           row.names = list(names_uro,names_all))
  }else{
    score_matrix <- cbind("Uro" = prediction$predictions[,"Uro"],
                          "UroA" = NA,
                          "UroB" = NA,
                          "UroC" = NA,
                          "GU" = prediction$predictions[,"GU"],
                          "BaSq" = prediction$predictions[,"BaSq"],
                          "Mes" = prediction$predictions[,"Mes"],
                          "ScNE" = prediction$predictions[,"ScNE"])
  }
  
  #calculate additional scores
  if(!subtype_only){
    if(verbose){
      message("Calculating signature scores...")
    }
  all_scores = lundtax_calc_sigscore(this_data = this_data,
                                     gene_id = gene_id,
                                     threshold_progression = threshold_progression,
                                     log_transform = log_transform,
                                     adjust = adjust,
                                     adj_factor = adj_factor,
                                     impute = impute, 
                                     impute_reject = impute_reject, 
                                     impute_kNN = impute_kNN,
                                     verbose = verbose)
  }else{
    if(verbose){
      message("Signature scores will be skipped, set subtype_only = FALSE to return scores...")
    }
    all_scores = NULL
  }
  
  #gather return
  #additional scores
  results_suburo$scores = all_scores
  
  #score matrices
  results_suburo$subtype_scores = score_matrix
  
  score_matrix_suburo = score_matrix[,2:4, drop = FALSE]
  score_matrix_5c = score_matrix[,c(1,5:8), drop = FALSE]
  
  #5 class level
  results_suburo$predictions_5classes = colnames(score_matrix_5c)[max.col(replace(score_matrix_5c,is.na(score_matrix_5c),-Inf), ties.method = "first")]
  
  #7 class level
  results_suburo$predictions_7classes = results_suburo$predictions_5classes
  max_suburo = colnames(score_matrix_suburo)[max.col(replace(score_matrix_suburo, is.na(score_matrix_suburo), -Inf), ties.method = "first")]
  
  for(i in 1:length(results_suburo$predictions_7classes)){
    p = results_suburo$predictions_7classes[i]
    if(p == "Uro"){
      suburo = max_suburo[i]
      results_suburo$predictions_7classes[i] = suburo
    }
  }
  
  #score ties
  #5 class level
  first5 = setNames(colnames(score_matrix_5c)[max.col(replace(score_matrix_5c, is.na(score_matrix_5c), -Inf), ties.method = "first")], rownames(score_matrix_5c))
  last5  = setNames(colnames(score_matrix_5c)[max.col(replace(score_matrix_5c, is.na(score_matrix_5c), -Inf), ties.method = "last")], rownames(score_matrix_5c))
  
  #run check.ties
  if(sum(first5 != last5) > 0){
    int_check_ties(first5,last5)
  }
  
  #7 class level
  score_matrix_suburo_ties = score_matrix_suburo[!is.na(score_matrix_suburo[,1]),]
  first7 = setNames(colnames(score_matrix_suburo_ties)[max.col(score_matrix_suburo_ties, ties.method = "first")], rownames(score_matrix_suburo_ties))
  last7 = setNames(colnames(score_matrix_suburo_ties)[max.col(score_matrix_suburo_ties, ties.method = "last")], rownames(score_matrix_suburo_ties))
  
  #run check.ties
  if(sum(first7 != last7) > 0){
    int_check_ties(first7,last7)
  }
  
  #final results
  names(results_suburo$predictions_7classes) = colnames(this_data)
  names(results_suburo$predictions_5classes) = colnames(this_data)
  
  if(subtype_only){
    predictions_suburo = list(predictions_7classes = results_suburo$predictions_7classes,
                              predictions_5classes = results_suburo$predictions_5classes)
    
    results_suburo_nodata = list(subtype_scores = results_suburo$subtype_scores,
                                 predictions_7classes = results_suburo$predictions_7classes,
                                  predictions_5classes = results_suburo$predictions_5classes)
  }else{
    predictions_suburo = list(predictions_7classes = results_suburo$predictions_7classes,
                              predictions_5classes = results_suburo$predictions_5classes,
                              scores = results_suburo$scores)
    
    results_suburo_nodata = list(subtype_scores = results_suburo$subtype_scores,
                                 predictions_7classes = results_suburo$predictions_7classes,
                                 predictions_5classes = results_suburo$predictions_5classes,
                                 scores = results_suburo$scores)
  }
  
  if(include_data & include_pred_scores){
    result = results_suburo
  }else if(include_data == FALSE & include_pred_scores){
    result = results_suburo_nodata
  }else if(include_data == FALSE & include_pred_scores == FALSE){
    result = predictions_suburo
  }
  
  if(verbose){
    message("Success!")
  }
  return(result)
}
