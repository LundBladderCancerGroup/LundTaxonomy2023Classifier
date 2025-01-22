#' @title Subset Predictions Data
#'
#' @description This function subsets the data frames and named vectors in a list based on a set of sample IDs.
#' 
#' @details This function is designed to subset various components of a predictions list based on a 
#' set of sample IDs. If `these_sample_ids` is not provided, the function will attempt to extract 
#' sample IDs from `these_samples_metadata`. The function handles different structures within the 
#' list, including data frames and named vectors, ensuring that only the specified samples are retained.
#'
#' @param these_predictions A list containing data frames and named vectors to be subsetted.
#' @param these_sample_ids A vector of sample IDs to subset the data. If NULL, sample IDs will be 
#' extracted from `these_samples_metadata`.
#' @param these_samples_metadata A data frame containing sample metadata. If provided and 
#' `these_sample_ids` is NULL, sample IDs will be extracted from this data frame.
#' @param samples_rownames A logical value indicating whether the sample IDs are stored in the row 
#' names of `these_samples_metadata`. Default is TRUE.
#'
#' @return A list with the same structure as `these_predictions`, but subsetted to the specified sample IDs.
#' 
#' @import dplyr tibble
#' 
#' @export
#'
#' @examples
#  sample_ids <- c("X18KFU0001", "X18KFU0002", "X18KFU0003")
#' subsetted_data <- subset_predictions(these_predictions = pred_uroscanseq_high_low, 
#'                                        these_sample_ids = sample_ids)
#'   
#' # Using metadata to extract sample IDs
#' subsetted_data <- subset_predictions(these_predictions = pred_uroscanseq_high_low, 
#'                                      these_samples_metadata = uc_high)
#'
subset_predictions <- function(these_predictions = NULL,
                               these_sample_ids = NULL,
                               these_samples_metadata = NULL,
                               samples_rownames = TRUE){
  
  if(!is.null(these_samples_metadata) && is.null(these_sample_ids)){
    message("Metadata provided, the function will subset to the sample IDs in this object...")
    if(samples_rownames){
      these_sample_ids = these_samples_metadata %>% 
        rownames_to_column("sample_id") %>% 
        pull(sample_id)
    }else{
      these_sample_ids = these_samples_metadata %>% 
        pull(sample_id)
    }
  }else{
    message("Sample IDs detected...")
  }
  
  #subset the 'data' data frame
  these_predictions$data <- these_predictions$data %>%
    select(all_of(these_sample_ids))
  
  #subset the 'subtype_scores' data frame
  these_predictions$subtype_scores <- these_predictions$subtype_scores[these_sample_ids, , drop = FALSE]
  
  #subset the 'predictions_7classes' named vector
  these_predictions$predictions_7classes <- these_predictions$predictions_7classes[these_sample_ids]
  
  #subset the 'predictions_5classes' named vector
  these_predictions$predictions_5classes <- these_predictions$predictions_5classes[these_sample_ids]
  
  #subset the 'scores' data frame
  these_predictions$scores <- these_predictions$scores %>%
    filter(rownames(these_predictions$scores) %in% these_sample_ids)
  
  return(these_predictions)
}
