#' @title Lund Taxonomy Scores.
#'
#' @description Wrapper function for calculating Lund Taxonomy scores. 
#'
#' @details This function internally calls a set of non-exported scoring functions starting with the
#' int_ prefix. See individual function documentation for more information on the individual function calls.
#'
#' @param this_data Required parameter. Data frame or matrix with expression values.
#' @param gene_id Specify the type of gene identifier used in `this_data`. Accepted values are; 
#' hgnc_symbol (default) or ensembl_gene_id.
#' @param threshold_progression Threshold to flag a sample as high risk of progression, default is 0.58.
#' @param log_transform Boolean parameter. If TRUE (default), the function log transforms the incoming expression
#' values.
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
#' 
#' @return A data frame with scores for the selected variable.
#' 
#' @import dplyr multiclassPairs
#
#' @examples
#' #load packages
#' library(dplyr, multiclassPairs)
#' 
#' #get scores for bundled data
#' my_scores = score_lundtax(this_data = sjodahl_2017, impute = TRUE)
#'                             
score_lundtax = function(this_data = NULL,
                         gene_id = "hgnc_symbol",
                         threshold_progression = 0.58,
                         logTransform = TRUE,
                         adjust = FALSE,
                         adj_factor = 5.1431,
                         impute = FALSE, 
                         impute_reject = 0.67, 
                         impute_kNN = 5, 
                         verbose = TRUE){
  
  #proliferation
  results_proliferation = int_ratio_score(this_data = this_data,
                                          variable = "proliferation",
                                          gene_id = gene_id, 
                                          logTransform = logTransform)
  
  #progression
  score_progression = int_ratio_score(this_data = this_data,
                                      variable = "progression",
                                      logTransform = logTransform,
                                      gene_id = gene_id)
  
  results_progression = ifelse(score_progression$Score >= threshold_progression, "HR", "LR")
  
  #grade
  #WHO 1999 (G3 vs G1/2)
  results_g3 = int_predict_grade(this_data = this_data, 
                                 grade_predictor = classifier_grade3, 
                                 gene_id = gene_id, 
                                 impute = impute, 
                                 impute_reject = impute_reject, 
                                 impute_kNN = impute_kNN, 
                                 verbose = verbose)

  #WHO 2004/2016 (HG vs LG)
  results_hg = int_predict_grade(this_data = this_data, 
                                 grade_predictor = classifier_hg, 
                                 gene_id = gene_id, 
                                 impute = impute, 
                                 impute_reject = impute_reject, 
                                 impute_kNN = impute_kNN, 
                                 verbose = verbose)

  #immune
  results_immune = int_calc_score(this_data = this_data, 
                                  variable = "immune", 
                                  logTransform = logTransform, 
                                  gene_id = gene_id, 
                                  adjust = adjust, 
                                  adj_factor = adj_factor)

  #141UP
  scores141up = int_calc_score(this_data = this_data, 
                               variable = "score141up", 
                               logTransform = logTransform, 
                               gene_id = gene_id, 
                               adjust = adjust, 
                               adj_factor = adj_factor)

  #add 141UP to immune results
  results_immune <- cbind(Immune141_UP=scores141up$Immune141_UP,
                          results_immune[,1:10,drop=FALSE],
                          Stromal141_UP=scores141up$Stromal141_UP,
                          results_immune[,11:13,drop=FALSE])
  
  #merge scores
  message("Merging scores...")
  merge_scores <- cbind(ProliferationScore=results_proliferation$Score,
                        MolecularGradeWHO1999=results_g3$predictions_classes,
                        MolecularGradeWHO1999_score=results_g3$predictions[,"G3"],
                        MolecularGradeWHO2016=results_hg$predictions_classes,
                        MolecularGradeWHO2016_score=results_hg$predictions[,"HG"],
                        ProgressionScore = score_progression$Score,
                        ProgressionRisk = results_progression,
                        results_immune)
  
  return(merge_scores)
}
