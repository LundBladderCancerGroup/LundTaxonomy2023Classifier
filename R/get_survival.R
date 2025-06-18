#' @title Get Survival.
#'
#' @description Run a Cox model to calculate hazard ratio with confidence intervalls.
#'
#' @details This function takes a data frame with prediction data `these_predictions` and executes 
#' a cox model to retrieve hazard ratio with confidence intervals based on a categorical variable (i.e response,
#' tumor grade, etc.) from a provided metadata table. The function expects the incoming 
#' data to be the output from [LundTax2023Classifier::lundtax_predict_sub()], together with metadata
#' information of interest (e.g two level categorical) and subtype classification information. The user
#' have the option to point the function to the categorical variables with `predictor_columns`. The return
#' can be further subset by subtype by using the `this_subtype` variable, should be one of the valid
#' subtypes within the specified class. The user is required to provide the funciton with correct 
#' column names in the metadata for survival time and survival event (see `surv_time` and `surv_event`).
#'
#' @param these_predictions Required parameter. A data frame with predictions scores, subtype and 
#' metadata information.
#' @param these_samples_metadata Required, a data frame with metadata associated with the prediction 
#' calls, Note, the function will subset the return to samples are included in this object.
#' @param subtype_class Can be one of the following; 5_class or 7_class. Default is 5_class.
#' @param this_subtype Optional parameter. Allows the user to subset the return to a specific subtype
#' within the selected class. If not specified, the function will return a data frame with subtype 
#' information for all the subtypes within the specified class.
#' @param scale Optional parameter. A numeric value to scale the numeric scores. If provided, all 
#' numeric scores will be multiplied by this value.
#' @param bin_scores Boolean parameter. Set to TRUE to bin the numeric scores into discrete bins. Default is FALSE.
#' @param n_bins Optional parameter. The number of bins to use when binning numeric scores. Default is 10.
#' @param predictor_columns Optional, should be a vector with column names, either from the provided 
#' metadata or signature score object, to be tested for. If not provided, the function will subset 
#' data to the signature scores returned with `lundtax_predict_sub`.
#' @param exclude_columns Optional argument, specify columns you wish to exclude from the standard 
#' predictor columns. Note, this parameter is only validated if predictor_columns is NULL (default).
#' @param surv_time Required parameter, should be the name of the column in the metadata with survival 
#' time. Should be of value numeric.
#' @param surv_event Required parameter, should be the name of the column in the metadata with survival 
#' event. Should be of value factor, with two levels.
#' @param row_to_col Boolean parameter. Set to TRUE to transform row names of the metadata to a new 
#' column called sample_id. Default is FALSE.
#' @param sample_id_col Parameter dictating the column name with sample IDs, the function expects this
#' column to be sample_id but the user can override this if they know the name for this column.
#'
#' @return A data frame with statistical scores for the selected samples/subtypes.
#'
#' @import dplyr survival
#' @rawNamespace import(stats, except = c(filter, lag))
#'
#' @export
#'
#' @examples
#' #load packages
#' library(dplyr, stats survival)
#' 
#' #get prediction calls
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                         impute = TRUE)
#'
#' #run general linear models
#' sjodahl_surv = get_survival(these_predictions = sjodahl_predicted,
#'                            these_samples_metadata = sjodahl_2017_meta,
#'                            subtype_class = "5_class",
#'                            this_subtype = "Uro",
#'                            surv_time = "surv_css_time",
#'                            surv_event = "surv_css_event")
#'
get_survival = function(these_predictions = NULL,
                        these_samples_metadata = NULL,
                        subtype_class = "5_class",
                        scale = NULL,
                        bin_scores = FALSE,
                        n_bins = 10,
                        this_subtype = NULL,
                        predictor_columns = NULL,
                        exclude_columns = NULL,
                        surv_time = NULL,
                        surv_event = NULL,
                        row_to_col = FALSE,
                        sample_id_col = NULL){
  
  #checks
  if(length(this_subtype) > 1){
    stop("Currently, only one subtype at the time is supported. For now, consider running this 
         function multiple times for each subtype. Or set this_subtype to NULL to return all subtypes 
         in the specified class...")
  }
  
  if(is.null(surv_time)){
    stop("The user must provide a column name corresponding to survival time...")
  }else{
    if(!surv_time %in% colnames(these_samples_metadata)){
      stop(paste0(surv_time, " is not a valid column in the incoming metadata..."))
      if(!is.numeric(these_samples_metadata[,surv_time])){
        stop("surv_event must be of value numeric...")
      }
    }
  }
  
  if(is.null(surv_event)){
    stop("The user must provide a column name corresponding to survival event...")
  }else{
    if(!surv_event %in% colnames(these_samples_metadata)){
      stop(paste0(surv_event, " is not a valid column in the incoming metadata..."))
      if(!is.factor(these_samples_metadata[,surv_event])){
        stop("surv_event must be a factor...")
        if(length(levels(these_samples_metadata[,surv_event])) != 2){
          stop("Levels of surv_event must be exactly 2...")
        }
      }
    }
  }
  
  #run helper
  this_object = int_prediction_wrangler(these_predictions = these_predictions,
                                        these_samples_metadata = these_samples_metadata,
                                        subtype_class = subtype_class,
                                        this_subtype = this_subtype,
                                        scale = scale, 
                                        bin_scores = bin_scores, 
                                        n_bins = n_bins,
                                        categorical_factor = categorical_factor,
                                        row_to_col = row_to_col,
                                        surv_time = surv_time, 
                                        surv_event = surv_event,
                                        sample_id_col = sample_id_col, 
                                        return_all = TRUE)
  
  if(is.null(predictor_columns)){
    #get the columns to test
    these_columns = c("proliferation_score", "progression_score", 
                      "stromal141_up", "immune141_up", "b_cells", "t_cells", "t_cells_cd8", 
                      "nk_cells", "cytotoxicity_score", "neutrophils", "monocytic_lineage", 
                      "macrophages", "m2_macrophage", "myeloid_dendritic_cells", "endothelial_cells", 
                      "fibroblasts", "smooth_muscle", "molecular_grade_who_1999_score", "molecular_grade_who_2022_score")
    
    #define all possible subtypes
    if(is.null(this_subtype)){
      if(subtype_class == "5_class"){
        all_subtypes = c("Uro", "GU", "BaSq", "Mes", "ScNE")
      }else if(subtype_class == "7_class"){
        all_subtypes = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")
      }
      
      these_columns = c(these_columns, all_subtypes)
      
      #exclude columns
      if(!is.null(exclude_columns)){
        if(any(!exclude_columns %in% colnames(this_object))){
          stop("One or more columns specified in exclude_columns is not in the incoming data...")
        }else{
          message(paste0("Excluding the following predictor columns: ", exclude_columns))
          these_columns = setdiff(these_columns, exclude_columns) 
        }
      }
    }

    
  }else{
    these_columns = predictor_columns
  }

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
  
  #convert column with survival event to numeric
  this_object[[surv_event]] = as.numeric(levels(this_object[[surv_event]])[this_object[[surv_event]]])
  
  #function to fit Cox models for multiple predictors
  fit_cox_models = function(my_data,
                            time_column,
                            event_column,
                            pred_columns){
    
    #check if the necessary columns exist in the data frame
    if (!all(c(time_column, event_column, pred_columns) %in% colnames(my_data))) {
      stop("One or more specified columns do not exist in the data frame.")
    }
    
    #fit Cox proportional hazards models using lapply
    cox_models = lapply(pred_columns, function(col){
      formula = as.formula(paste("Surv(", time_column, ", ", event_column, ") ~", col))
      coxph(formula, data = my_data)
    })
    
    #name the list elements with the predictor column names
    names(cox_models) = pred_columns
    
    return(cox_models)
  }
  
  cox_models = fit_cox_models(my_data = this_object, 
                              time_column = surv_time, 
                              event_column = surv_event, 
                              pred_columns = these_columns)
  
  #helper function for running the stats on all scores
  extract_surv_stats = function(cox){
    
    #initialize an empty data frame to store the results
    my_stats = data.frame()
    
    #loop over each element in the cox object
    for(name in names(cox)){
      #create a new data frame for the current statistic
      stats = data.frame(score = name)
      
      #extract hazard ratio and confidence intervals from Cox model
      cox_model = cox[[name]]
      stats$p_value = summary(cox_model)$coefficients[5]
      stats$hazard_ratio = exp(coef(cox_model))
      conf_int = exp(confint(cox_model))
      stats$hazard_conf_2.5 <- conf_int[1]
      stats$hazard_conf_97.5 <- conf_int[2]
      
      #bind the current results to the overall statistics data frame
      my_stats = rbind(my_stats, stats)
    }
    return(my_stats)
  }
  
  #run helper
  stat_results = extract_surv_stats(cox = cox_models)
  
  if(!is.null(this_subtype)){
    stat_results$subtype = as.factor(this_subtype)
  }else{
    stat_results$subtype = as.factor("all")
  }
  
  #return results
  return(stat_results)
}
