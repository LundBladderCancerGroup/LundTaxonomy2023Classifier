#' @title Check Ties.
#'
#' @description Check score ties between predictions scores.
#'
#' @details Internal function called by [LundTax2023Classifier::lundtax_calc_sigscore()]. 
#' Not meant for out of package use. Checks score ties in the prediction scores and prints a message 
#' indicating the sample where the tie occurred, the two subtypes with the tied scores and the 
#' subtype that is reported in the output object.
#'
#' @param first Predictions when setting ties.method = "first"
#' @param last Predictions when setting ties.method = "last"
#'
#' @return Nothing.
#'
int_check_ties = function(first, last){
  for(tie in which(first != last)){
    sample = names(first)[tie]
    subtype1 = first[tie]
    subtype2 = last[tie]
    message(paste0("Sample ",sample, ": tie between ", subtype1, " and ", subtype2, "\nOutput prediction is ", subtype1))
  }
}
