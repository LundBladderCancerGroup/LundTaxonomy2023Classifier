#' @title Plot Ranked Score
#'
#' @description Return a point plot with ranked scores for a set variable.
#'
#' @details This function takes predictions returned with 
#' [LundTax2023Classifier::lundtax_predict_sub()], together with a score variable,
#' subtype class and returns a point plot (grob) with ranked scores along the x axis and 
#' scores along the y axis. The returned plot will also color fill the points based on subtype 
#' classification.
#'
#' @param these_predictions Required parameter. Predictions returned with 
#' [LundTax2023Classifier::lundtax_predict_sub()].
#' @param this_score Required parameter, should be one of the numeric columns in the scores data 
#' frame from the list returned with [LundTax2023Classifier::lundtax_predict_sub()].
#' @param subtype_class Required, one of the following; 5_class or 7_class. Needed for coloring the 
#' points based on subtype classification.
#' @param plot_title Optional parameter, the name of the returned plot. If NULL (default), no title 
#' will be printed to the generated figure.
#'
#' @return Nothing.
#' 
#' @import ggplot2 dplyr
#' @importFrom tibble rownames_to_column
#'
#' @export
#'
#' @examples
#' sjodahl_2017_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                              impute = TRUE)
#'
#' plot_ranked_score(these_predictions = sjodahl_2017_predicted, 
#'                   this_score = "proliferation_score", 
#'                   subtype_class = "7_class", 
#'                   plot_title = "Sjodahl 2017")
#'
plot_ranked_score = function(these_predictions = NULL, 
                             this_score = NULL, 
                             subtype_class = NULL, 
                             plot_title = NULL){
  
  #checks
  if(is.null(these_predictions)){
    stop("these_predictions is missing, should be the return from lund_predict_sub...")
  }
  
  #impute plot title
  if(is.null(plot_title)){
    plot_title = this_score
  }
  
  if(is.null(this_score)){
    stop("this_score is missing, should be a numeric variable in the scores data frame within these_predictions...")
  }else if(!this_score %in% colnames(these_predictions$scores)){
    stop(paste0(this_score, " is non-existing in the scores data frame within these_predictions..."))
  }else if(!is.numeric(these_predictions$scores[,this_score])){
    stop(paste0(this_score, " is not of type numeric..."))
  }
  
  #subset on subtype
  if(subtype_class == "5_class"){
    this_subtype = data.frame(these_predictions$predictions_5classes) %>% 
      tibble::rownames_to_column("sample_id") 
  }else if(subtype_class == "7_class"){
    this_subtype = data.frame(these_predictions$predictions_7classes) %>% 
      tibble::rownames_to_column("sample_id") 
  }else{
    stop("A valid subtype class must be provided, one of the following; 5_class or 7_class...")
  }
  
  #subset score, left join with subtypes, add rank
  score_sub = these_predictions$scores %>% dplyr::select(all_of(this_score)) %>% 
    tibble::rownames_to_column("sample_id") %>% 
    dplyr::left_join(this_subtype, by = "sample_id") %>% 
    dplyr::rename(score = 2, subtype = 3) %>% 
    dplyr::mutate_at(vars(subtype), factor) %>% 
    mutate(rank = as.numeric(as.factor(score)))

  #draw plot
  my_plot = ggplot(data = score_sub, aes(x = rank, y = score, color = subtype)) +
    geom_point(size = 2, shape = 16) +
    scale_color_manual(values = lund_colors$lund_colors) + 
    labs(title = plot_title, subtitle = paste0("All samples (n = ", nrow(score_sub), ")")) +
    xlab("Index") + 
    ylab(this_score) +
    theme_classic()
  
  print(my_plot)
  
  return()
}
