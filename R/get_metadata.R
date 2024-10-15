#' @title Get Metadata
#'
#' @description Return filtered metadata, or sample IDs from a provided metadata data frame 
#' based on a set of criteria.
#'
#' @details Takes a metadata data frame as the main input and allows easy subset for various variables.
#' For more information, see parameter descriptions and examples.
#'
#' @param this_metadata Required parameter, should be a data frame with metadata information
#' @param these_sample_ids Optional parameter. If provided, the incoming metadata will be restricted 
#' to the sample IDs provided here. Should be a vector of sample IDs of interest.
#' @param return_type Should be on of the following; metadata (default) or sample_id. If sample IDs 
#' are selected, a vector with selected sample IDs will be returned, nothing else.
#' @param rows_to_column Boolean parameter. Set to TRUE to convert row names to a new column called 
#' sample_id. Default is FALSE.
#' @param sample_column Optional, lets the user define a specific column in the data frame with sample IDs.
#' @param first_variable Filtering options, should be the name of the column of interest.
#' @param first_value Filtering options, the filter value for the first variable.
#' @param second_variable Filtering options, should be the name of the column of interest.
#' @param second_value Filtering options, the filter value for the second variable.
#' @param third_variable Filtering options, should be the name of the column of interest.
#' @param thid_value Filtering options, the filter value for the third variable.
#' 
#' @return A data frame with subset metadata, or a vector of characters with sample IDs.
#' 
#' @import dplyr
#'
#' @export
#'
#' @examples
#' #example 1
#' #return all sample IDs from bundled metadata
#' get_metadata(this_metadata = sjodahl_2017_meta, 
#'              return_type = "sample_id")
#' 
#' #use more filter options and return metadata data frame
#' get_metadata(this_metadata = sjodahl_2017_meta, 
#'              first_variable = "gender", 
#'              first_value = "Female", 
#'              second_variable = "turb_stage", 
#'              second_value = 1, 
#'              third_variable = "adj_chemo", 
#'              third_value = 0)
#'
get_metadata = function(this_metadata = NULL,
                        these_sample_ids = NULL,
                        return_type = "metadata",
                        rows_to_column = FALSE,
                        sample_column = NULL,
                        first_variable = NULL,
                        first_value = NULL,
                        second_variable = NULL,
                        second_value = NULL,
                        third_variable = NULL,
                        third_value = NULL){
  
  #check parameters
  if(is.null(this_metadata)){
    stop("No metadata is provided with `this_metadata`, this is a required parameter for this function...")
  }
  
  #rename the column with sample IDs
  if(!is.null(sample_column)){
    this_metadata = this_metadata %>%
      rename(sample_column = "sample_id")
  }
  
  #convert rownames to sample ID column
  if(rows_to_column){
    this_metadata = tibble::rownames_to_column(this_metadata, "sample_id")
  }
  
  #check if there is a valid column with sample IDs
  if(!"sample_id" %in% colnames(this_metadata)){
    stop("There is no column named 'sample_id' in the provided metadata, please adjust this. 
         See `rows_to_col` or `sample_col` for more info...")
  }
  
  if(is.null(these_sample_ids)){
    message("No sample IDs provided, the function will subset the metadata based on the provided filtering criterias...")
  }else{
    this_metadata = dplyr::filter(this_metadata, sample_id %in% these_sample_ids)
    
    #check if metadata is empty
    if(nrow(this_metadata) == 0){
      stop("Provided sample ID(s) were not found in the provided metadata, try a different sample ID...")
    }
    
    #check the existence of provided sample IDs in the metadata
    not_in_meta = setdiff(these_sample_ids, this_metadata$sample_id)
    
    if(length(not_in_meta) > 0){
      message("WARNING! The following sample IDs were not found in the metadata: ")
      print(not_in_meta)
    }
  }
  
  #filter the metadata
  if(!is.null(first_variable) && !is.null(first_value)){
    if(!first_variable %in% colnames(this_metadata)){
      message(paste0(first_variable, " is not a valid column in the incoming metadata, no filtering will take place..."))
    }
    this_metadata = this_metadata %>% 
      dplyr::filter(!!as.symbol(first_variable) %in% first_value)
  }
  
  if(!is.null(second_variable) && !is.null(second_value)){
    if(!second_variable %in% colnames(this_metadata)){
      message(paste0(second_variable, " is not a valid column in the incoming metadata, no filtering will take place..."))
    }
    this_metadata = this_metadata %>% 
      dplyr::filter(!!as.symbol(second_variable) %in% second_value)
  }
  
  if(!is.null(third_variable) && !is.null(third_value)){
    if(!third_variable %in% colnames(this_metadata)){
      message(paste0(third_variable, " is not a valid column in the incoming metadata, no filtering will take place..."))
    }
    this_metadata = this_metadata %>% 
      dplyr::filter(!!as.symbol(third_variable) %in% third_value)
  }
  
  if(return_type == "metadata"){
    return(this_metadata)
  }else if(return_type == "sample_id"){
    return(this_metadata$sample_id)
  }else{
    stop("Valid vlaues are 'metadata' and 'sample_id'...")
  }
}
 