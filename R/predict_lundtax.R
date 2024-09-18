#' Predict Lund Taxonomy subtypes based on rule-based Random Forest classifiers
#'
#' @param data matrix, data frame or multiclassPairs_object of gene expression values
#' @param subtype_only 
#' @param include_data include data in output (disabled by default)
#' @param include_pred_scores include prediction scores for each sample and class in output (default)
#' @param gene_id specify the type of gene identifier used in the data:
#' - "hgnc_symbol" for HUGO gene symbols
#' - "ensembl_gene_id" for Ensembl gene IDs
#' - "entrezgene" for Entrez IDs
#' Default value is hgnc_symbol
#' @param logTransform h
#' @param adjust hh
#' @param adj_factor default is 5.1431
#' @param ... Additional parameters to be passed to the predict_RF function. If genes are missing in the data, include impute = TRUE here
#' @return
#' Returns a list object including:
#' - Data (optional, not included by default)
#' - Prediction scores for all classes (optional, included by default)
#' - Predicted LundTax class for 7-class system
#' - Predicted LundTax class for 5-class system
#'
#' @details
#' This function uses 2 classifiers to classify the samples: 5-class classifier first  classifies samples into Uro, GU, BaSq, Mes or ScNE.
#' Samples classified as Uro receive a second classification as UroA, B or C by the second classifier
#'
#' @examples
#' # Include data in result
#' results_data <- predict_lundtax(Lund2017,
#'                                 impute = TRUE,
#'                                 include_data = TRUE)
#'
#' @export
#
predict_lundtax <- function(data,
                            subtype_only = FALSE, # include only subtype prediction (no additional scores)
                            include_data = FALSE, # return input data in the results object
                            include_pred_scores = TRUE, # return prediction scores in the results object
                            gene_id = c("hgnc_symbol","ensembl_gene_id")[1],
                            scoring_method = c("ratio","singscore")[1],
                            logTransform = TRUE, # Scores are calculated on log transformed data. If the data is already log transformed, set logTransform to FALSE. If logTranform = TRUE, data will be log2 transformed (log2(data+1)) before calculating the scores
                            adjust = TRUE, # adjust scores by stable genes and adjustment factor
                            adj_factor = 5.1431,# adjustment factor
                            verbose = TRUE,
                            ...)
                                    
  
{
  require(multiclassPairs)
  # Check inputs
  
  ## Data ##
  # Store as dataframe
  if (!(class(data)[1] %in% c("matrix","data.frame","multiclassPairs_object"))) {
    stop("Data should be in one of the following formats: matrix, data.frame, multiclassPairs_object")
  }
  
  if (ncol(data) != length(unique(colnames(data)))) {
    stop("Sample names (column names) should not be duplicated")
  }
  
  if (class(data)[1] == "multiclassPairs_object") {
    D <- data$data$Data
    # Get ref labels
    if (is.null(ref)) {
      ref <- data$data$Labels
    }
  }
  
  if (is.data.frame(data)) {
    D <- as.matrix(data)
  }
  
  if (is.matrix(data)) {
    # D <- as.data.frame(data, stringsAsFactors = FALSE)
    D <- data
  }
  
  # Check gene identifiers ####
  if (gene_id != "hgnc_symbol") {
    
    original_D <- D
    
    # # Testing
    # load("D:/Packages/LundTaxonomy2023Classifier_DEV/gene_info_lund.rda")
    gene_info_lund <- LundTax2023Classifier::gene_info_lund
    
    rownames(gene_info_lund) <- gene_info_lund[[gene_id]]
    int_genes <- rownames(D)[which(rownames(D) %in% gene_info_lund[[gene_id]])]
    rownames(D)[which(rownames(D) %in% gene_info_lund[[gene_id]])] <- gene_info_lund[int_genes,"hgnc_symbol"]
    
  } else {
    original_D <- D
    # change "-" to "_" to avoid errors in RF
    if (TRUE %in% grepl("-",rownames(D))) rownames(D) <- gsub("-","_",rownames(D))
    
    
  }
  
  # Classifier ##
  C <- LundTax2023Classifier::LundTax_RF_5c
  C2 <- LundTax2023Classifier::LundTax_RF_Uro7c
  
  # Results object ##
  
  
  results_suburo <- list(data = original_D,
                         subtype_scores = NULL,
                         predictions_7classes = NULL,
                         predictions_5classes = NULL,
                         scores = NULL)
  
  ## Predict 5 class ###
  
  prediction <- predict_RF(classifier = C,
                           Data = D,
                           verbose = verbose, ...)
  
  # Reorder scores
  pred <- prediction$predictions[,c("Uro","GU","BaSq","Mes","ScNE"), drop=FALSE]
  prediction$predictions <- pred
  prediction$predictions_classes <- colnames(pred)[max.col(replace(pred,is.na(pred),-Inf),ties.method = "first")]
  
  ## Get uro samples
  if ("Uro" %in% prediction$predictions_classes) {
    D_Uro <- D[,which(prediction$predictions_classes == "Uro"), drop = FALSE]
    D_NoUro <- D[,which(prediction$predictions_classes != "Uro"), drop = FALSE]
    
    # Classify suburo if necessary ###
    prediction_suburo <- predict_RF(classifier = C2,
                                    Data = D_Uro,
                                    verbose = verbose, ...)
    
    names_uro <- colnames(D_Uro)
    names_all <- colnames(D)
    
    # Merged score matrix
    score_matrix <- merge_subUro_matrix(score_matrix1 = prediction_suburo$predictions,
                                        score_matrix2 = prediction$predictions,
                                        row.names = list(names_uro,names_all))
  } else { # if there is no Uro
    score_matrix <- cbind("Uro" = prediction$predictions[,"Uro"],
                          "UroA" = NA,
                          "UroB" = NA,
                          "UroC" = NA,
                          "GU" = prediction$predictions[,"GU"],
                          "BaSq" = prediction$predictions[,"BaSq"],
                          "Mes" = prediction$predictions[,"Mes"],
                          "ScNE" = prediction$predictions[,"ScNE"])
    
  }
  
  # Calculate additional scores ##
  if (!subtype_only) {
  all_scores <- score_lundtax(Data = original_D,
                              gene_id = gene_id,
                              scoring_method = scoring_method,
                              logTransform = logTransform,
                              adjust = adjust,
                              adj_factor = adj_factor,
                              verbose = verbose,
                              ...)
  } else {
    all_scores <- NULL
  }
  
  # Collect results ##
  
  # Additional scores
  results_suburo$scores <- all_scores
  
  # Score matrices
  results_suburo$subtype_scores <- score_matrix
  
  score_matrix_suburo <- score_matrix[,2:4, drop = FALSE]
  score_matrix_5c <- score_matrix[,c(1,5:8), drop = FALSE]
  
  # 5 class level
  results_suburo$predictions_5classes <- colnames(score_matrix_5c)[max.col(replace(score_matrix_5c,is.na(score_matrix_5c),-Inf),ties.method = "first")]
  
  # 7 class level
  results_suburo$predictions_7classes <- results_suburo$predictions_5classes
  max_suburo <- colnames(score_matrix_suburo)[max.col(replace(score_matrix_suburo,is.na(score_matrix_suburo),-Inf),ties.method = "first")]
  
  for (i in 1:length(results_suburo$predictions_7classes)) {
    p <- results_suburo$predictions_7classes[i]
    if (p == "Uro") {
      suburo <- max_suburo[i]
      results_suburo$predictions_7classes[i] <- suburo
    }
    
  }
  
  # Score ties ##
  
  # 5 class level
  first5 <- setNames(colnames(score_matrix_5c)[max.col(replace(score_matrix_5c,is.na(score_matrix_5c),-Inf),ties.method = "first")],rownames(score_matrix_5c))
  last5  <- setNames(colnames(score_matrix_5c)[max.col(replace(score_matrix_5c,is.na(score_matrix_5c),-Inf),ties.method = "last")],rownames(score_matrix_5c))
  
  if (sum(first5 != last5)>0) {
    check.ties(first5,last5)
    
  }
  
  # 7 class level
  score_matrix_suburo_ties <- score_matrix_suburo[!is.na(score_matrix_suburo[,1]),]
  first7 <- setNames(colnames(score_matrix_suburo_ties)[max.col(score_matrix_suburo_ties,ties.method = "first")],rownames(score_matrix_suburo_ties))
  last7  <- setNames(colnames(score_matrix_suburo_ties)[max.col(score_matrix_suburo_ties,ties.method = "last")],rownames(score_matrix_suburo_ties))
  
  
  if (sum(first7 != last7)>0) {
    check.ties(first7,last7)
  }
  
  # Final results ##
  names(results_suburo$predictions_7classes) <- colnames(D)
  names(results_suburo$predictions_5classes) <- colnames(D)
  
  if (subtype_only) {
    predictions_suburo <- list(predictions_7classes = results_suburo$predictions_7classes,
                               predictions_5classes = results_suburo$predictions_5classes)
    
    results_suburo_nodata <- list(subtype_scores = results_suburo$subtype_scores,
                                  predictions_7classes = results_suburo$predictions_7classes,
                                  predictions_5classes = results_suburo$predictions_5classes)
  } else {
    predictions_suburo <- list(predictions_7classes = results_suburo$predictions_7classes,
                               predictions_5classes = results_suburo$predictions_5classes,
                               scores = results_suburo$scores)
    
    results_suburo_nodata <- list(subtype_scores = results_suburo$subtype_scores,
                                  predictions_7classes = results_suburo$predictions_7classes,
                                  predictions_5classes = results_suburo$predictions_5classes,
                                  scores = results_suburo$scores)
  }
  
  if (include_data & include_pred_scores) {
    result <- results_suburo
  } else if (include_data == FALSE & include_pred_scores) {
    result <- results_suburo_nodata
  } else if (include_data == FALSE & include_pred_scores == FALSE) {
    result <- predictions_suburo
  }
  
}



