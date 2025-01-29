#' @title Prediction Wrangler
#'
#' @description Internal function called by `get_glm_scores`. 
#' Not meant for out-of-package usage.
#'
#' @details Takes an output from [LundTax2023Classifier::lundtax_predict_sub()] together with 
#' associated metadata and pre-processing it so that statistical tests can be executed on this data.
#' For more info, see the documentation of `get_glm_scores`.
#'
#' @param these_predictions Required. Output from [LundTax2023Classifier::lundtax_predict_sub()].
#' @param these_samples_metadata Required. Metadata associated with he prediction output. Also possible
#' for the user to provide a metadata subset with samples of interest, the return will be restricted
#'  to the samples within the specified group.
#' @param sample_id_col Optional parameter dictating the column name with sample IDs, the function 
#' expects this column to be sample_id but the user can override this if they know the name for this
#' column.
#' @param row_to_col Boolean parameter. Set to TRUE to transform row names of the metadata to a new 
#' column called sample_id. Default is FALSE.
#' @param categorical_factor Required parameter. This should be the categorical variable that is intended 
#' for testing. In addition, this should also be a variable of type factor, with exactly 2 levels.
#' @param scale Optional parameter. A numeric value to scale the numeric scores. If provided, all 
#' numeric scores will be multiplied by this value.
#' @param bin_scores Boolean parameter. Set to TRUE to bin the numeric scores into discrete bins. Default is FALSE.
#' @param n_bins Optional parameter. The number of bins to use when binning numeric scores. Default is 10.
#' @param surv_time Optional, should be the column name for the survival time (numeric) in the metadata.
#' @param surv_event Optional, should be the column name for survival event (factor) in the metadata.
#' @param subtype_class Can be one of the following; 5_class or 7_class. Default is 5_class.
#' @param this_subtype Optional parameter. Allows the user to subset the return to a specific subtype
#' within the selected class. If not specified, the function will return a data frame based on all 
#' subtypes in the spcified class.
#' @param return_all Boolean parameter, set to TRUE to return all metadata columns. Default is FALSE.
#'
#' @return A data frame ready for `get_glm`, or `get_survival`.
#'
#' @import dplyr
#'
#' @examples
#' #load packages
#' library(dplyr)
#' 
#' #get prediction calls
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                         impute = TRUE)
#'   
#' #run helper                                  
#' my_out = int_prediction_wrangler(these_predictions = sjodahl_predicted, 
#'                                  these_samples_metadata = sjodahl_2017_meta, 
#'                                  subtype_class = "5_class", 
#'                                  this_subtype = "Uro", 
#'                                  categorical_factor = "adj_chemo")
#'                           
int_prediction_wrangler = function(these_predictions = NULL,
                                   these_samples_metadata = NULL,
                                   sample_id_col = NULL,
                                   row_to_col = FALSE,
                                   subtype_class = "5_class",
                                   categorical_factor = NULL,
                                   scale = NULL,
                                   bin_scores = FALSE,
                                   n_bins = 10,
                                   surv_time = NULL,
                                   surv_event = NULL,
                                   this_subtype = NULL, 
                                   return_all = FALSE){
    
    #check predictions
    if(is.null(these_predictions)){
      stop("You must provide a valid dataset with predictions scores!")
    }

    #subset to the score columns of interest
    my_scores = these_predictions$scores %>% 
      tibble::rownames_to_column("sample_id")
    
    if(!is.null(scale)){
      my_scores = my_scores %>% 
        mutate(across(where(is.numeric), ~ . * scale))
    }
    
    if(bin_scores){
      my_scores = int_bin_numeric_variables(this_data = my_scores, 
                                            num_bins = n_bins)
    }

    #get the subtype information
    if(subtype_class == "5_class"){
      my_subtypes = as.data.frame(these_predictions$predictions_5classes) %>% 
        tibble::rownames_to_column("sample_id")
    }else if(subtype_class == "7_class"){
      my_subtypes = as.data.frame(these_predictions$predictions_7classes) %>% 
        tibble::rownames_to_column("sample_id")
    }else{
      stop("only 5_class and 7_class subtypes are available...")
    }
    
    #rename subtype column
    colnames(my_subtypes)[2] = "subtype"
    
    #filter on subtypes
    if(!is.null(this_subtype)){
      message("Filtering prediction data on the selected subtype(s)...")
      my_subtypes = dplyr::filter(my_subtypes, subtype %in% this_subtype)
    }
    
    #subset to other relevant metadata columns
    if(is.null(these_samples_metadata)){
      stop("User must provide metadata with these_samples_metadata...")
    }else{
      if(row_to_col){
        these_samples_metadata = these_samples_metadata %>% 
          tibble::rownames_to_column("sample_id")
      }
      
      if(is.null(surv_time) && is.null(surv_event)){
        #check the categorical variable
        if(!categorical_factor %in% colnames(these_samples_metadata)){
          stop(paste0(categorical_factor, " is not a valid column in the incoming metadata..."))
          if(!is.factor(these_samples_metadata[,categorical_factor])){
            stop("categorical_factor must be a factor...")
            if(length(levels(these_samples_metadata[,categorical_factor])) != 2){
              stop("Levels of categorical_factor must be exactly 2...")
            }
          }
        } 
      }
      
      #check sample ID in metadata
      if(!"sample_id" %in% names(these_samples_metadata)){
        if(is.null(sample_id_col)){
          stop("Metadata has no column namned sample_id, consider running the function with sample_col, specifying the name of the column with sample IDs...")
        }else{
          colnames(these_samples_metadata)[colnames(these_samples_metadata) == sample_id_col] = "sample_id"
        }
      }
      
      if(!return_all){
        if(is.null(surv_time) && is.null(surv_event)){
          my_metadata = these_samples_metadata %>%
            dplyr::select(sample_id, categorical_factor)
        }else{
          message("Adding survival data to metadata")
          my_metadata = these_samples_metadata %>%
            dplyr::select(sample_id, categorical_factor, surv_time, surv_event)
        }
      }else{
        my_metadata = these_samples_metadata
      }
    }

    #combine data
    my_object = left_join(my_subtypes, my_scores, by = "sample_id")
    my_object = inner_join(my_metadata, my_object, by = "sample_id")
    message(paste0(nrow(my_object), " samples kept after filtering on subtype and metadata..."))
    
    #warn user how many samples are in the metadata and how many are kept after filtering
    no_scores = setdiff(my_metadata$sample_id, my_scores$sample_id)
    
    if(length(no_scores > 1)){
      message(paste0("Warning, ", length(no_scores), " samples are removed from the metadata, since they have no scores reported..."))
    }
    
    if(is.null(surv_time) && is.null(surv_event)){
      #reinstate factor type in categorical column
      my_object[,categorical_factor] = as.factor(my_object[,categorical_factor])
      
      #check that levels of subset data and remove subtypes that does not have two factors in the categorical_factor
      if(length(levels(my_object[,categorical_factor])) != 2){
        message("The resulting subset filter does not have 2 levels in the selected categorical_factor...")
      }
      
      #check if there are actually two levels present of the selected variable
      if(length(unique(my_object[,categorical_factor])) != 2){
        message("There are not two levels present in the selected column after filtering steps..")
        message("Returning empty data frame...")
        my_object[,] = matrix(ncol = ncol(my_object), rep(NA, prod(dim(my_object))))
        empty_object = my_object[complete.cases(my_object), ]
        return(empty_object)
      }
    }
    
    #define all possible subtypes
    if(subtype_class == "5_class"){
      all_subtypes = c("Uro", "GU", "BaSq", "Mes", "ScNE")
    }else if(subtype_class == "7_class"){
      all_subtypes = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")
    }
    
    #convert the 'subtype' column to a factor with all possible levels
    my_object$subtype = factor(my_object$subtype, levels = all_subtypes)
    
    #create a binary matrix for subtype
    binary_matrix = model.matrix(~ subtype - 1, data = my_object)
    
    #convert the binary matrix to a data frame
    binary_df = as.data.frame(binary_matrix)
    colnames(binary_df) = gsub("subtype", "", colnames(binary_df))
    
    #combine the binary columns with the original data frame
    my_object = cbind(my_object, binary_df)
    
    #ensure subtype are factors
    my_object[all_subtypes] = lapply(my_object[all_subtypes], factor)
    
    #convert specified columns to factors with both levels "0" and "1"
    my_object[all_subtypes] = lapply(my_object[all_subtypes], function(x) {
      factor(x, levels = c(0, 1))
    })
    
    #re-level categorical columns
    my_object$progression_risk = relevel(my_object$progression_risk, ref = "LR")
    my_object$molecular_grade_who_2022 = relevel(my_object$molecular_grade_who_2022, ref = "LG")
    my_object$GU = relevel(my_object$GU, ref = "0")
    my_object$BaSq = relevel(my_object$BaSq, ref = "0")
    my_object$Mes = relevel(my_object$Mes, ref = "0")
    my_object$ScNE = relevel(my_object$ScNE, ref = "0")
    
    if(subtype_class == "5_class"){
      my_object$Uro = relevel(my_object$Uro, ref = "0")
    }else if(subtype_class == "7_class"){
      my_object$UroA = relevel(my_object$UroA, ref = "0")
      my_object$UroB = relevel(my_object$UroB, ref = "0")
      my_object$UroC = relevel(my_object$UroC, ref = "0")
    }
    
    return(my_object)
}
