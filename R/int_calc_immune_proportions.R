#' @title Calculate Immune Proportions.
#'
#' @description Calculate immune score proportions.
#'
#' @details Internal function called by `in_calc_score` if varaible is set to immune.
#'
#' @param immune_results Required parameter. Data frame with immune scores from int_calc_score.
#' 
#' @return A data frame with scores as percentages.
#
#' @examples
#' #calculate immune scores from hgnc symbols 
#' immune_scores = int_calc_score(this_data = sjodahl_2017, 
#'                                 variable = "immune", 
#'                                 gene_id = "hgnc_symbol")
#' 
#' immune_prop = int_calc_immune_proportions(immune_results = immune_scores)
#
int_calc_immune_proportions = function(immune_results){
  immune_proportions <- t(apply(immune_results,1,function(x){x/sum(x)}))
  colnames(immune_proportions) <- paste0(colnames(immune_proportions)," Proportion")
  as.data.frame(immune_proportions)
  return(immune_proportions)
}
