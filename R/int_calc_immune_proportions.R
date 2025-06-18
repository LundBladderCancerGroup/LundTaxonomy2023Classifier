#' @title Calculate Immune Proportions.
#'
#' @description Calculate immune score proportions.
#'
#' @details Internal function called by [LundTax2023Classifier::int_calc_score()], if variable is 
#' set to immune. Not meant for out-of-package use.
#'
#' @param immune_results Required parameter. Data frame with immune scores from 
#' [LundTax2023Classifier::int_calc_score()].
#' 
#' @return A data frame with scores as percentages.
#'
int_calc_immune_proportions = function(immune_results){
  immune_proportions <- t(apply(immune_results,1,function(x){x/sum(x)}))
  colnames(immune_proportions) <- paste0(colnames(immune_proportions)," Proportion")
  as.data.frame(immune_proportions)
  return(immune_proportions)
}
