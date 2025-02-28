#' @title Get Sample Metrics
#'
#' @description This function tables a metadata column based on subtype classification.
#' 
#' @details The function performs the following steps:
#' \itemize{
#'   \item Ensures the sample IDs are present in both metadata and predictions.
#'   \item Checks if the desired metadata column is valid and is a factor.
#'   \item Filters metadata and predictions to include only common samples.
#'   \item Combines metadata and predictions into a single data frame.
#'   \item Counts the number of samples and progression events for each subtype.
#' }
#'
#' @param this_metadata A data frame containing metadata with a `sample_id` column.
#' @param this_metadata_variable A string specifying the column name in the metadata to be used for progression events.
#' @param these_predictions A named vector of predictions with sample IDs as names.
#' @param factor_level The level of the factor in `this_metadata_variable` to summarize (default is 1).
#' @param subtype_class The classification system, default is 5 class.
#'
#' @return A data frame with the number of samples and progression events for each subtype.
#' 
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' my_predictions = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                      gene_id = "hgnc_symbol", 
#'                                      impute = TRUE, 
#'                                      adjust = TRUE)
#'
#' get_sample_metrics(these_predictions = my_predictions, 
#'                    this_metadata = sjodahl_2017_meta, 
#'                    this_metadata_variable = "surv_os_event", 
#'                    factor_level = 1)
#'
get_sample_metrics <- function(this_metadata = NULL,
                               this_metadata_variable = NULL,
                               these_predictions = NULL,
                               factor_level = NULL,
                               subtype_class = "5_class") {
  
  # Check if sample_id is present in metadata
  if (!"sample_id" %in% colnames(this_metadata)) {
    stop("The metadata must contain a 'sample_id' column.")
  }
  
  # Check if sample_id is present in predictions
  if (is.null(names(these_predictions))) {
    stop("The predictions must have sample IDs as names.")
  }
  
  # Check if the desired metadata column is valid
  if (!this_metadata_variable %in% colnames(this_metadata)) {
    stop("The specified metadata variable is not a valid column in the metadata.")
  }
  
  # Check if the metadata column is a factor
  if (!is.factor(this_metadata[[this_metadata_variable]])) {
    stop("The specified metadata variable must be a factor.")
  }
  
  if(subtype_class == "5_class"){
    #ensure the sample IDs are the same in both metadata and predictions
    common_samples <- intersect(names(these_predictions$predictions_5classes), this_metadata$sample_id)
    
    # Filter metadata and predictions to include only common samples
    filtered_metadata <- this_metadata %>% filter(sample_id %in% common_samples)
    filtered_predictions <- these_predictions$predictions_5classes[common_samples]
  }else if(subtype_class == "7_class"){
    #ensure the sample IDs are the same in both metadata and predictions
    common_samples <- intersect(names(these_predictions$predictions_7classes), this_metadata$sample_id)
    
    # Filter metadata and predictions to include only common samples
    filtered_metadata <- this_metadata %>% filter(sample_id %in% common_samples)
    filtered_predictions <- these_predictions$predictions_7classes[common_samples]
  }

  # Combine metadata and predictions into a single data frame
  combined_data <- filtered_metadata %>%
    mutate(prediction = filtered_predictions[match(sample_id, names(filtered_predictions))])
  
  # Count the number of samples and progression events for each subtype
  result <- combined_data %>%
    group_by(prediction) %>%
    summarise(
      num_samples = n(),
      num_progression_events = sum(.data[[this_metadata_variable]] == factor_level)
    )
  
  return(result)
}
