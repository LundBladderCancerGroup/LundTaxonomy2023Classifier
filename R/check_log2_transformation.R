#' @title Check if RNAseq Expression Data is Log2-Transformed
#' 
#' @description This function checks whether the values in an RNAseq expression data frame have been 
#' log2-transformed. It performs multiple checks, including the range of values, the presence of zeros, 
#' and the presence of negative values, to provide a robust assessment.
#'
#' @details
#' The function uses the following criteria to determine if the data is log2-transformed:
#' - If the maximum value exceeds 20, the data is likely not log2-transformed.
#' - If the minimum value is less than -5, the data is likely not log2-transformed.
#' - If the data contains zeros, it may be log2-transformed but without pseudocounts.
#' - If the range and distribution of values are consistent with log2-transformed data, 
#'   the function concludes that the data is likely log2-transformed.
#'
#' @param expression_df A data frame containing RNAseq expression data.
#' @param print_histogram Set to TRUE to draw histogram with expression values, default is FALSE.
#'
#' @return A list with information about log2 transformation.
#' 
#' @examples
#' check_log2_transformation(expression_df = sjodahl_2017)
#' 
#' @export
#' 
check_log2_transformation = function(expression_df, 
                                     print_histogram = FALSE){
  
  #flatten the data frame into a numeric vector
  values <- as.numeric(as.matrix(expression_df))
  
  #remove NA values (if any)
  values <- values[!is.na(values)]
  
  #initialize a results list
  results <- list()
  
  #check for zeros
  results$contains_zeros <- any(values == 0)
  
  #calculate the range of values
  results$min_value <- min(values)
  results$max_value <- max(values)
  
  #calculate the mean and median
  results$mean_value <- mean(values)
  results$median_value <- median(values)
  
  #check for negative values
  results$contains_negative <- any(values < 0)
  
  #perform checks to determine if the data is log2-transformed
  if (results$max_value > 20) {
    results$log2_transformed <- FALSE
    results$message <- "The data is NOT log2-transformed. The maximum value is too large."
  } else if (results$min_value < -5) {
    results$log2_transformed <- FALSE
    results$message <- "The data is NOT log2-transformed. The minimum value is too small (likely not log2-transformed)."
  } else if (results$contains_zeros) {
    results$log2_transformed <- TRUE
    results$message <- "The data is likely log2-transformed but contains zeros (check if pseudocounts were added)."
  } else {
    results$log2_transformed <- TRUE
    results$message <- "The data is likely log2-transformed. The range and distribution are consistent with log2-transformed data."
  }
  
  #return the results
  return(results)
}
