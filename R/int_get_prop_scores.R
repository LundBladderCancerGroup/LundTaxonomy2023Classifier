#' @title Get Proportional Scores
#' 
#' @description Internal function called by `plot_hm_scores` when `proportional_scores = TRUE`.
#' Not meant for out-of-package usage.
#'
#' @details This function normalizes specific prediction scores so that each value is a fraction of 
#' the total score for that row. It also calculates the immune/stroma ratio.
#'
#' @param these_predictions A list containing a data frame named `scores` with prediction scores. 
#' The data frame should have row names as sample IDs and columns corresponding to different cell types and scores.
#'
#' @return A data frame with the original scores, normalized scores, and the calculated 
#' immune/stroma ratio. The row names are the sample IDs.
#'
#' @import dplyr tibble
#' 
#' @examples
#' #load packages
#' library(dplyr)
#' 
#' #get prediction calls
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                         impute = TRUE)
#'                                         
#' #transform data into proportional signature scores
#' prop_scores = int_get_prop_scores(these_predictions = sjodahl_predicted)
#' 
int_get_prop_scores = function(these_predictions = NULL){
  
  #get names for columns to normalize
  norm_signatures <- c("b_cells", "t_cells", "t_cells_cd8", "nk_cells", "cytotoxicity_score", 
                       "neutrophils", "monocytic_lineage", "macrophages", "m2_macrophage", "myeloid_dendritic_cells", 
                       "endothelial_cells", "fibroblasts", "smooth_muscle")
  
  #subset scores that are to be normalized
  norm_scores = these_predictions$scores %>% 
    select(all_of(norm_signatures))
  
  #normalize each row so that each value is a fraction of the total score for that row
  norm_scores <- norm_scores %>%
    mutate(across(everything(), ~ . / rowSums(across(everything()), na.rm = TRUE)))
  
  #set sample id for joining
  norm_scores <- norm_scores %>%
    rownames_to_column(var = "sample_id")
  
  #subset original scores object and set sample ID
  original_scores <- these_predictions$scores %>%
    rownames_to_column(var = "sample_id")
  
  #merge the data frames on the sample_id column
  updated_scores <- original_scores %>%
    left_join(norm_scores, by = "sample_id", suffix = c("", ".new"))
  
  #get the column names for updateing og scores (exclude sample_id)
  columns_to_update <- colnames(norm_scores)[-1]
  
  #remove the temporary columns and restore row names
  updated_scores <- updated_scores %>%
    select(-ends_with(".new")) %>%
    column_to_rownames(var = "sample_id")
  
  #calculate stroma/immune ratio
  updated_scores$immune_stroma_ratio = updated_scores$immune141_up - updated_scores$stromal141_up
  
  return(updated_scores)
}
