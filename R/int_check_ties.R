#' @title Check Ties.
#'
#' @description Internal function called by `score_lundtax`. 
#' Not meant for out of package use.
#'
#' @details Checks score ties in the prediction scores and prints a message indicating
#' the sample where the tie occurred, the two subtypes with the tied scores and 
#' the subtype that is reported in the output object.
#'
#' @param first Predictions when setting ties.method = "first"
#' @param last Predictions when setting ties.method = "last"
#'
#' @return Nothing.
#'
check_ties = function(first, last){
  for(tie in which(first != last)){
    sample = names(first)[tie]
    subtype1 = first[tie]
    subtype2 = last[tie]
    message(paste0("Sample ",sample, ": tie between ", subtype1, " and ", subtype2, "\nOutput prediction is ", subtype1))
  }
}
