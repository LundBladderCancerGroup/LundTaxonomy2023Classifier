#' @title Get General Linear Model Scores.
#'
#' @description Run a Mann-Whitney and general linear regression on a set of
#' prediction scores.
#'
#' @details This function takes a data frame with prediction data `these_predictions` and executes 
#' statistical tests to find significant prediction scores based on a categorical variable (i.e response,
#' tumour grade, etc.) from a provided metadata table. Currently, the function expects the incoming 
#' data to be the score output from [LundTax2023Classifier::lundtax_predict_sub()], together with metadata
#' information of interest (e.g two level categorical) and subtype classification information. The user
#' have the option to point the function to the categorical variable with `categorical_factor`. The return
#' can be further subset by subtype by using the `this_subtype` variable, should be one of the valid
#' subtypes within the specified class.
#'
#' @param these_predictions Required parameter. A data frame with predictions scores, subtype and 
#' metadata information.
#' @param these_samples_metadata Required, a data frame with metadata associated with the prediction 
#' calls, Note, the function will subset the return to samples are included in this object.
#' @param subtype_class Can be one of the following; 5_class or 7_class. Default is 5_class.
#' @param this_subtype Optional parameter. Allows the user to subset the return to a specific subtype
#' within the selected class. If not specified, the function will return a data frame with subtype 
#' information for all the subtypes within the specified class.
#' @param categorical_factor Required parameter. Specify the two level categorical variable you want to test for.
#' @param row_to_col Boolean parameter. Set to TRUE to transform row names of the metadata to a new 
#' column called sample_id. Default is FALSE.
#' @param sample_id_col Parameter dictating the column name with sample IDs, the function expects this
#' column to be sample_id but the user can override this if they know the name for this column.
#'
#' @return A data frame with statistical scores for the selected samples/subtypes.
#'
#' @import dplyr
#' @rawNamespace import(stats, except = c(filter, lag))
#'
#' @export
#'
#' @examples
#' #load packages
#' library(dplyr, stats)
#' 
#' #get prediction calls
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                         impute = TRUE)
#'
#' #run general linear models
#' sjodahl_glm = get_glm_scores(these_predictions = sjodahl_predicted,
#'                              these_samples_metadata = sjodahl_2017_meta,
#'                              subtype_class = "5_class",
#'                              this_subtype = "Uro",
#'                              categorical_factor = "adj_chemo")
#'
get_glm_scores = function(these_predictions = NULL,
                          these_samples_metadata = NULL,
                          subtype_class = "5_class",
                          this_subtype = NULL,
                          categorical_factor = NULL,
                          row_to_col = FALSE,
                          sample_id_col = NULL){
  
  if(length(this_subtype) > 1){
    stop("Currently, only one subtype at the time is supported. For now, consider running this function multiple times for each subtype...")
  }

  #run helper
  this_object = int_prediction_wrangler(these_predictions = these_predictions,
                                        these_samples_metadata = these_samples_metadata,
                                        subtype_class = subtype_class,
                                        this_subtype = this_subtype,
                                        categorical_factor = categorical_factor,
                                        row_to_col = row_to_col,
                                        surv_time = NULL, 
                                        surv_event = NULL,
                                        sample_id_col = sample_id_col)

  #get the columns to test
  these_columns = colnames(dplyr::select_if(this_object, is.numeric))
  
  #remove all proportion columns
  these_columns = these_columns[!grepl(x = these_columns, pattern = "_proportion")]
  
  #check if data frame is empty and return it
  if(nrow(this_object) == 0){
    empty_df = data.frame(score = character(),
                          p_value = numeric(),
                          odds_ratio = numeric(),
                          conf_2.5 = numeric(),
                          conf_97.5 = numeric(),
                          subtype = factor())
    return(empty_df)
  }
  
  #mann whitney u-test
  message("Running Mann-Whitney U-test...")
  wilcox = lapply(this_object[,these_columns], function(x) wilcox.test(x ~ this_object[,categorical_factor]))

  #generalized linear model
  message("Running Generalized Linear Model...")
  glm = lapply(this_object[,these_columns], function(x) glm(x ~ this_object[,categorical_factor]))
  
  #helper function for running the stats on all scores
  extract_stats = function(wilcox, 
                           glm){
    
    #initialize an empty data frame to store the results
    my_stats = data.frame()
    
    #loop over each element in the wilcox object
    for (name in names(wilcox)) {
      #create a new data frame for the current statistic
      stats = data.frame(score = name)
      
      #extract p-value, odds ratio, and confidence intervals
      stats$p_value = wilcox[[name]]$p.value
      stats$odds_ratio = unname(unlist(exp(suppressMessages(glm[[name]]$coefficients))))[2]
      stats$conf_2.5 = unname(unlist(exp(suppressMessages(confint(glm[[name]])))))[2]
      stats$conf_97.5 = unname(unlist(exp(suppressMessages(confint(glm[[name]])))))[4]
      
      #bind the current results to the overall statistics data frame
      my_stats = rbind(my_stats, stats)
    }
    return(my_stats)
  }

  #run helper
  stat_results = extract_stats(wilcox = wilcox,
                                glm = glm) 


  if(!is.null(this_subtype)){
    stat_results$subtype = as.factor(this_subtype)
  }else{
    stat_results$subtype = as.factor("all")
  }

  #return results
  return(stat_results)
}
