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
#' my_scores = score_lundtax(this_data = sjodahl_2017, 
#'                           impute = TRUE)
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

  #merge scores
  message("Merging scores...")
  merge_scores <- cbind(proliferation_score = results_proliferation$Score,
                        molecular_grade_who_1999 = results_g3$predictions_classes,
                        molecular_grade_who_1999_score = results_g3$predictions[,"G3"],
                        molecular_grade_who_2016 = results_hg$predictions_classes,
                        molecular_grade_who_2016_score = results_hg$predictions[,"HG"],
                        progression_score = score_progression$Score,
                        progression_risk = results_progression,
                        results_immune, 
                        scores141up)
  
  #set new column names for easy access
  oldnames = c("B-cells", "T-cells","T-cells CD8+", "NK-cells", "Cytotoxicity Score", "Neutrophils",
               "Monocytic lineage", "Macrophages","M2 macrophage", "Myeloid Dendritic Cells", 
               "Endothelial cells", "Fibroblasts", "Smooth muscle", "B-cells Proportion", 
               "T-cells Proportion", "T-cells CD8+ Proportion","NK-cells Proportion", 
               "Cytotoxicity Score Proportion", "Neutrophils Proportion",
               "Monocytic lineage Proportion", "Macrophages Proportion","M2 macrophage Proportion", 
               "Myeloid Dendritic Cells Proportion", "Endothelial cells Proportion", 
               "Fibroblasts Proportion", "Smooth muscle Proportion", "Immune141_UP", "Stromal141_UP")
  
  newnames = c("b_cells", "t_cells","t_cells_cd8", "nk_cells", "cytotoxicity_score", "neutrophils",
               "monocytic_lineage", "macrophages","m2_macrophage", "myeloid_dendritic_cells", 
               "endothelial_cells", "fibroblasts", "smooth_muscle", "b_cells_proportion", 
               "t_cells_proportion", "t_cells_cd8_proportion","nk_cells_proportion", 
               "cytotoxicity_score_proportion", "neutrophils_proportion",
               "monocytic_lineage_proportion", "macrophages_proportion","m2_macrophage_proportion", 
               "myeloid_dendritic_cells_proportion", "endothelial_cells_proportion", 
               "fibroblasts_proportion", "smooth_muscle_proportion", "immune141_up", "stromal141_up")
  
  merge_scores = merge_scores %>% 
    rename_at(vars(oldnames), ~ newnames)
  
  #set column types
  factor_cols = c("molecular_grade_who_1999", "molecular_grade_who_2016", "progression_risk")
  merge_scores[factor_cols] = lapply(merge_scores[factor_cols], factor)
  
  #set custom order of columns
  merge_scores = merge_scores %>% 
    dplyr::select(proliferation_score, progression_score, progression_risk, 
                  molecular_grade_who_2016_score, molecular_grade_who_2016, 
                  molecular_grade_who_1999_score,molecular_grade_who_1999, stromal141_up, immune141_up, 
                  b_cells, b_cells_proportion, t_cells, t_cells_proportion,
                  t_cells_cd8, t_cells_cd8_proportion, nk_cells, nk_cells_proportion,
                  cytotoxicity_score, cytotoxicity_score_proportion, neutrophils, neutrophils_proportion,
                  monocytic_lineage, monocytic_lineage_proportion, macrophages, macrophages_proportion,
                  m2_macrophage, m2_macrophage_proportion, myeloid_dendritic_cells, 
                  myeloid_dendritic_cells_proportion, endothelial_cells, endothelial_cells_proportion,
                  fibroblasts, fibroblasts_proportion, smooth_muscle, smooth_muscle_proportion)
  
  return(merge_scores)
}
