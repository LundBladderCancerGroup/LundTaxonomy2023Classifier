#' @title Merge subUro Matrix.
#'
#' @description Merge prediction score matrices from two classifiers.
#'
#' @description Internal function called by [LundTax2023Classifier::lundtax_calc_sigscore()]. 
#' Not meant for out of package use. This function merges the prediction score matrices from the 
#' 5-class and 7-class (UroA,UroB,UroC) classifiers into 1 unique score matrix.
#'
#' @param score_matrix1 Prediction score matrix from the 7 class classifier.
#' @param score_matrix2 Prediction score matrix from the 5 class classifier.
#' @param row.names Rownames (sample names) for both matrices.
#'
#' @return  Merged matrix including scores for 8 classes (Uro, UroA, UroB, UroC, 
#' GU, BaSq, Mes, ScNE)
#'
int_merge_suburo_matrix = function(score_matrix1,
                                   score_matrix2,
                                   row.names = list(NULL, NULL)){
  #check inputs
  #matrices
  if(!(class(score_matrix1)[1] %in% c("matrix","data.frame") & class(score_matrix2) %in%  c("matrix","data.frame"))[1]){
    stop("Inputs should be score matrices from a multiclassPairs::rule_based_RandomForest object")
  }
  
  #row names
  if(!is.null(row.names) & length(row.names) == 2){
    rownames(score_matrix1) <- row.names[[1]]
    rownames(score_matrix2) <- row.names[[2]]
    name_rows <- rownames(score_matrix2)
  }else{
    stop("Row names for both matrices are needed.")
  }
  
  #create empty matrix
  #number and name of rows
  n_rows <- nrow(score_matrix2)
  final_classes = c("Uro","UroA","UroB","UroC","GU","BaSq","Mes","ScNE")
  final_matrix = matrix(nrow = n_rows, ncol = length(final_classes), dimnames = list(name_rows, final_classes))
  
  for(i in colnames(final_matrix)){
    if(i %in% colnames(score_matrix1) & i %in% colnames(score_matrix2)){
      message("Classes should be just on one of the provided matrices. Values for duplicate classes will be taken from matrix 1.")
    }
    else if(i %in% colnames(score_matrix2)){
      final_matrix[rownames(score_matrix2),i] <- score_matrix2[,i]
    }else if(i %in% colnames(score_matrix1)){
      final_matrix[rownames(score_matrix1),i] <- score_matrix1[,i]
    }
  }
  return(final_matrix)
}
