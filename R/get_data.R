#' @title Get Data
#'
#' @description subset and filter incoming data based on a set of criteria. 
#'
#' @details Takes a metadata data frame as the main input and allows easy subset for various variables.
#' For more information, see parameter descriptions and examples.
#'
#' @param this_metadata Required parameter, should be a data frame with metadata information
#' @param these_sample_ids Optional parameter. If provided, the incoming metadata will be restricted 
#' to the sample IDs provided here. Should be a vector of sample IDs of interest.
#' @param return_type Should be on of the following; metadata (default), sample_id or expression_data. If sample IDs 
#' are selected, a vector with selected sample IDs will be returned, nothing else. If expression_data 
#' is selected, the user needs to specify the incoming expression data with `this_data`.
#' @param this_data Required parameter if return_type is set to `expression_data`. Should be a data 
#' frame with sample IDs as column names.
#' @param rows_to_column Boolean parameter. Set to TRUE to convert row names to a new column called 
#' sample_id. Default is FALSE.
#' @param sample_column Optional, lets the user define a specific column in the data frame with sample IDs.
#' @param exclude_id If metadata is provided and no sample IDs are specified, the user can control 
#' if any sample IDs should be omitted from the return. This parameter expects a character of vectors.
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
#' get_data(this_metadata = sjodahl_2017_meta, 
#'          return_type = "sample_id")
#' 
#' #use more filter options and return metadata data frame
#' get_data(this_metadata = sjodahl_2017_meta, 
#'          first_variable = "gender", 
#'          first_value = "Female", 
#'          second_variable = "turb_stage", 
#'          second_value = 1, 
#'          third_variable = "adj_chemo", 
#'          third_value = 0)
#'
get_data = function(this_metadata = NULL,
                    these_sample_ids = NULL,
                    this_data = NULL,
                    return_type = "metadata",
                    rows_to_column = FALSE,
                    sample_column = NULL,
                    exclude_id = NULL,
                    first_variable = NULL,
                    first_value = NULL,
                    second_variable = NULL,
                    second_value = NULL,
                    third_variable = NULL,
                    third_value = NULL){
  
  #check parameter combinations
  if(!is.null(these_sample_ids) && !is.null(exclude_id)){
    these_sample_ids = setdiff(these_sample_ids, exclude_id)
  }
  
  if(return_type == "metadata" && is.null(this_metadata)){
    stop("Metadata was requested as the return, but no incoming metadata is provided...")
  }
  
  #SCENARIO 1: user does not provide metadata or sample IDs
  if(is.null(this_metadata) && is.null(these_sample_ids)){
    stop("No metadata or sample IDs are provided. One or both are required for the funciton...")
  }
  
  #SCENARIO 2: user provides metadata (regardless of if sample IDs are provided or not)
  if(!is.null(this_metadata)){
    message("Metadata detected...")
    
    #rename the column with sample IDs
    if(!is.null(sample_column)){
      this_metadata = this_metadata %>%
        rename(sample_column = "sample_id")
    }
    
    #convert row names to sample ID column
    if(rows_to_column){
      this_metadata = tibble::rownames_to_column(this_metadata, "sample_id")
    }
    
    #check if there is a valid column with sample IDs
    if(!"sample_id" %in% colnames(this_metadata)){
      stop("There is no column named 'sample_id' in the provided metadata, please adjust this. 
         See `rows_to_col` or `sample_col` for more info...")
    }
    

    #remove specific samples
    if(!is.null(exclude_id)){
      this_metadata = dplyr::filter(this_metadata, !sample_id %in% exclude_id)
    }
    
    #filter the metadata, filter 1
    if(!is.null(first_variable) && !is.null(first_value)){
      if(!first_variable %in% colnames(this_metadata)){
        message(paste0(first_variable, " is not a valid column in the incoming metadata, 
                       no filtering will take place..."))
      }
      this_metadata = this_metadata %>% 
        dplyr::filter(!!as.symbol(first_variable) %in% first_value)
    }
    
    #filter the metadata, filter 2
    if(!is.null(second_variable) && !is.null(second_value)){
      if(!second_variable %in% colnames(this_metadata)){
        message(paste0(second_variable, " is not a valid column in the incoming metadata, 
                       no filtering will take place..."))
      }
      this_metadata = this_metadata %>% 
        dplyr::filter(!!as.symbol(second_variable) %in% second_value)
    }
    
    #filter the metadata, filter 3
    if(!is.null(third_variable) && !is.null(third_value)){
      if(!third_variable %in% colnames(this_metadata)){
        message(paste0(third_variable, " is not a valid column in the incoming metadata, 
                       no filtering will take place..."))
      }
      this_metadata = this_metadata %>% 
        dplyr::filter(!!as.symbol(third_variable) %in% third_value)
    }
    
    #get all sample IDs from the metadata if user did not provide sample IDs
    if(is.null(these_sample_ids)){
      sample_ids = this_metadata$sample_id
    }else{
      sample_ids = these_sample_ids
    }
  }
  
  #SCENARIO 3: user provides only sample IDs
  if(!is.null(these_sample_ids) && is.null(this_metadata)){
    message("No metadata provided, only sample IDs...")
    sample_ids = these_sample_ids
  }
  
  #SCENARIO 4: user provides both sample IDs and metadata
  if(!is.null(these_sample_ids) && !is.null(this_metadata)){
    message("Both sample IDs and metadata is provided, the funciton will subset the provided metadata to the provided sample IDs...")
    
    #get samples not in meta, before filtering
    not_in_meta = setdiff(these_sample_ids, this_metadata$sample_id)
    
    #filter the metadata to available sample IDs
    this_metadata = dplyr::filter(this_metadata, sample_id %in% these_sample_ids)
    
    #check if metadata is empty
    if(nrow(this_metadata) == 0){
      stop("Provided sample ID(s) were not found in the provided metadata, try a different sample ID...")
    }
    
    #subset to sample IDs in the metadata
    sample_ids = this_metadata$sample_id
    
    if(length(not_in_meta) > 0){
      message("WARNING! The following sample IDs were not found in the metadata: ")
      print(not_in_meta)
    }
  }
  
  if(return_type == "metadata"){
    return(this_metadata)
  }else if(return_type == "sample_id"){
    return(sample_ids)
  }else if(return_type == "expression_data"){
    if(is.null(this_data)){
      stop("User must provide expression data if the function should subset it to samples of interest...")
    }else{
      this_data = this_data %>% 
        dplyr::select(sample_ids)
      
      #check if the requested sample IDs are actually in the expression data
      expression_samples = colnames(this_data)
      not_in_data = setdiff(sample_ids, expression_samples)
      
      if(length(not_in_data) > 0){
        message("WARNING! The following sample IDs were not found in the expression data:")
        print(not_in_data)
      }
      return(this_data)
    }
    
  }else{
    stop("Valid vlaues are 'metadata', 'sample_id' and 'expression_data'...")
  }
}
 