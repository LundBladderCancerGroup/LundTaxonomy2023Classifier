
#' ReadData
#' Function to read the data and Labels
#' @export
ReadData <- function(Data,
                     Labels,
                     Platform = NULL,
                     verbose = TRUE) {

  # check the input Data format
  if (!is.data.frame(Data) &
      !is.matrix(Data) &
      class(Data)[1] != "ExpressionSet") {
    stop("Bad format for the input Data...should be:
         matrix, data.frame, or ExpressionSet")
  }

  if (is.data.frame(Data)) {
    Data_tmp <- Data
  }

  if (is.matrix(Data)) {
    Data_tmp <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # Biobase package is needed
    if(!requireNamespace("Biobase", quietly = TRUE)){
      message("ExpressionSet is used and 'Biobase' package from Bioconductor is needed!")
      stop("Visit their website or install Biobase package using:
      if (!requireNamespace('BiocManager', quietly = TRUE)) {
      install.packages('BiocManager')
      }
      BiocManager::install('Biobase')", call. = FALSE)
    } else {
      requireNamespace("Biobase")
    }

    # extract the expression matrix from the ExpressionSet
    Data_tmp <- as.data.frame(Biobase::exprs(Data), stringsAsFactors = FALSE)

    # if labels are not provided then give the available variables in Eset
    if (!hasArg(Labels)) {
      message(capture.output(cat("Phenotype data has these variables:",
                                 Biobase::varLabels(Data),
                                 fill = TRUE)))
      stop("input a vector with same length of samples number
           or select one of these variable for Labels")
    }

    # extract the Labels - in case it is stored in the ExpressionSet
    if (is.character(Labels) & length(Labels) == 1) {

      if (Labels %in% Biobase::varLabels(Data)) {
        Labels_tmp <- as.character(Biobase::pData(Data)[, Labels])

      } else {
        message(capture.output(cat("Phenotype data has these variables:",
                                   Biobase::varLabels(Data),
                                   fill = TRUE)))
        stop("Labels variable is not found in the phenotype data of your ExpressionSet")
      }
    }

    # get the input Labels vector as it is
    if ((is.character(Labels) | is.factor(Labels)) & length(Labels) != 1) {
      Labels_tmp <- as.character(Labels)

      if (length(Labels_tmp) != ncol(Data_tmp)) {
        message("Number of samples: ", ncol(Data_tmp))
        message("Labels length: ", length(Labels_tmp))
        stop("Labels vector length are not equal to samples in data")
      }
    }

    # if user input platform name or vector
    if (!is.null(Platform)) {
      # extract the Labels - in case it is stored in the ExpressionSet
      if (is.character(Platform) & length(Platform) == 1) {

        if (Platform %in% Biobase::varLabels(Data)) {
          Platform_tmp <- as.character(Biobase::pData(Data)[, Platform])

        } else {
          message(capture.output(
            cat("Phenotype data has these variables:",
                Biobase::varLabels(Data), fill = TRUE)))
          stop("Platform variable is not found in the phenotype
               data of your ExpressionSet")
        }
      }

      # get the input Platform vector as it is
      if ((is.character(Platform) |
           is.factor(Platform)) &
          length(Platform) != 1) {

        Platform_tmp <- as.character(Platform)

        if (length(Platform_tmp) != ncol(Data_tmp)) {
          message("Number of samples:", ncol(Data_tmp))
          message("Labels length:", length(Platform_tmp))
          stop("Platform vector length are not equal to samples in data")
        }
      }
    }
  }

  # check if rownames is not NULL to avoid error later
  if (is.null(rownames(Data_tmp))) {
    stop("Provide features/genes names as rownames in the Data matrix!")
  }

  # get the input Labels vector as it is
  if ((is.character(Labels) |
       is.factor(Labels)) &
      class(Data)[1] != "ExpressionSet") {

    Labels_tmp <- as.character(Labels)

    if (length(Labels_tmp) != ncol(Data_tmp)) {
      message(paste("Number of samples: ", ncol(Data_tmp)))
      message(paste("Labels length: ", length(Labels_tmp)))
      stop("Labels vector length are not equal to samples in data")
    }
  }

  # get the input Platform vector as it is
  if (!is.null(Platform)) {
    if ((is.character(Platform) |
         is.factor(Platform)) &
        class(Data)[1] != "ExpressionSet") {

      Platform_tmp <- as.character(Platform)

      if (length(Platform_tmp) != ncol(Data_tmp)) {
        message(paste("Number of samples: ", ncol(Data_tmp)))
        message(paste("Labels length: ", length(Platform_tmp)))
        stop("Platform vector length are not equal to samples in data")
      }
    }
  } else {
    Platform_tmp <- NULL
  }

  ###
  # Remove genes with NAs in all samples
  remove_na <- rowSums(is.na(Data_tmp)) == ncol(Data_tmp)

  if (sum(remove_na) > 0) {
    message(paste("These features will be removed because they have NA values
                  in all samples:"))
    message(paste(rownames(Data_tmp)[remove_na], collapse = " "))
    Data_tmp <- Data_tmp[!remove_na, ]
  }

  # Remove genes with NAs in all genes
  remove_na <- colSums(is.na(Data_tmp)) == nrow(Data_tmp)

  if (sum(remove_na) > 0) {
    message(paste("These samples will be removed because they have NA values
                  for all features:"))
    message(paste(colnames(Data_tmp)[remove_na], collapse = " "))
    Data_tmp   <- Data_tmp[, !remove_na]
    Labels_tmp <- Labels_tmp[!remove_na]

    if (!is.null(Platform)) {
      Platform_tmp <- Platform_tmp[!remove_na]
    }
  }

  # print info about the input data and labels
  if (verbose) {
    message("Creating Data object...")
    message("Number of samples: ", ncol(Data_tmp))
    message("Number of genes/features: ", nrow(Data_tmp))
    message(capture.output(cat("Classes:", unique(Labels_tmp), fill = TRUE)))

    if (!is.null(Platform)) {
      message(capture.output(cat("Platforms/studies:",
                                 unique(Platform_tmp),
                                 fill = TRUE)))
    } else {
      message("Platforms/studies: NULL")
    }
  }

  # if any labels are NAs then stop
  if (any(is.na(as.character(Labels_tmp)))) {
    stop("NAs are not allowed in labels!")
  }

  if (any(is.na(as.character(Platform_tmp)))) {
    stop("NAs are not allowed in Platform!")
  }


  # give warning message if the gene names have "-"
  # This will give errors in RF models and Boruta and ranger
  if (length(grep(x = rownames(Data_tmp), pattern = "-")) > 0) {
    message("Gene names in the data have '-' symbol! This may generate errors during the training process of random forest! It is recommended to change these '-' to '_' or '.'")
  }

  if (length(grep(x = rownames(Data_tmp), pattern = ",")) > 0) {
    message("Gene names in the data have ',' symbol! This may generate errors during the training process of random forest! It is recommended to change these ',' to '_' or '.'")
  }

  # create the object
  object <- list(
    data = list(Data=Data_tmp,
                Labels=Labels_tmp,
                Platform=Platform_tmp)
  )
  class(object) <- "multiclassPairs_object"

  return(object)
}

#' predict_RF
#' Predict sample class based on gene pair-based random forest classifier
#'
#' @import ranger
#' @import rdist
#' @export

predict_RF <- function(classifier,
                       Data,
                       impute = FALSE,
                       impute_reject = 0.67,
                       impute_kNN = 5, # for kNN
                       verbose = TRUE) {

  # check the object class
  if (!class(Data)[1] %in% c("multiclassPairs_object",
                             "ExpressionSet",
                             "data.frame",
                             "matrix")) {
    stop("Data should be class:
    matrix/data.frame/ExpressionSet/multiclassPairs_object from ReadData function!")
  }

  # check classifier object
  if (class(classifier)[1] != "rule_based_RandomForest") {
    stop("classifier should be rule_based_RandomForest object from train_RF function!")
  }

  if (!is.numeric(impute_reject) |
      !length(impute_reject) == 1 |
      any(impute_reject >= 1) |
      any(impute_reject <= 0)) {
    stop("impute_reject argument should be a number between 0 and 1!")
  }

  # get the data matrix
  if (is.data.frame(Data)) {
    D <- Data
  }

  if (is.matrix(Data)) {
    D <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # Biobase package is needed
    if(!requireNamespace("Biobase", quietly = TRUE)){
      message("ExpressionSet is used and 'Biobase' package from Bioconductor is needed!")
      stop("Visit their website or install Biobase package using:
      if (!requireNamespace('BiocManager', quietly = TRUE)) {
      install.packages('BiocManager')
      }
      BiocManager::install('Biobase')", call. = FALSE)
    } else {
      requireNamespace("Biobase")
    }

    # extract the expression matrix from the ExpressionSet
    D <- as.data.frame(Biobase::exprs(Data), stringsAsFactors = FALSE)
  }

  if (class(Data)[1]  ==  "multiclassPairs_object") {
    D <- Data$data$Data
  }

  # extract the genes and the rules
  genes <- classifier$RF_scheme$genes
  rules <- classifier$RF_scheme$rules

  # check if all genes are in the data
  if (any(!genes %in% rownames(D))) {

    if (verbose){
      message("These genes are not found in the data:")
      message(capture.output(cat(genes[!genes %in% rownames(D)])))
      message("Gene names should as rownames and sample names as columns!")
      message("Check the genes in classifier object to see all the needed genes.")
      message("Check if '-' or ',' symbols in the gene names in your data. You may need to change it to '_' or '.'")
    }

    if (impute == FALSE) {
      stop("All genes should be in the data with no NA values! Or you can turn impute argument to TRUE to impute missed genes to the closest class for each sample!")
    }

    if (impute == TRUE & verbose) {
      message("Missed genes will be imputed to the closest class for each sample!")
    }
  }

  # create empty matrix for the data
  complete <- data.frame(matrix(data = NA,
                                nrow = length(genes),
                                ncol = ncol(D),
                                dimnames = list(genes, colnames(D))),
                         check.names = FALSE,
                         stringsAsFactors = FALSE)

  # fill it with data
  found <- genes[genes %in% rownames(D)]
  complete[found,colnames(D)] <- D[found,]

  # Remove genes with NAs
  if (sum(!complete.cases(complete)) > 0) {
    if (verbose){
      message("These genes have NAs:")
      message(paste(rownames(complete)[!complete.cases(complete)],
                    collapse = " "))
    }

    if (impute == FALSE) {
      message("Turn impute to TRUE to impute NAs to the closest class for each sample with NAs!")
      stop("Gene which is used in the classifier should not have NAs!")
    }

    if (impute == TRUE & verbose) {
      message("These genes will be imputed to the closest class for each sample with NAs")
    }
  }

  # produce the binary matrix
  binary <- complete[rules$gene1, , drop=FALSE] < complete[rules$gene2, , drop=FALSE]
  rownames(binary) <- paste0(rules$gene1,
                             "__",
                             rules$gene2)


  #Impute if needed
  if (impute) {

    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    # get the mode values from the training data - stored in the classifier object
    # mode_df <- classifier$RF_scheme$mode
    TrainingMatrix_df <- classifier$RF_scheme$TrainingMatrix

    # to store the index for the samples to be removed
    # due to lack of a lot of rules
    to_remove_sam <- c()

    for(i in 1:ncol(binary)){

      # get which rules are missed in this sample
      is_na <- is.na(binary[,i])

      # if everything is OK then go to the next sample
      if (sum(is_na) == 0) {
        next
      }

      # skip the sample if it misses >0.67 of rules
      if (sum(is_na) > (nrow(binary)*impute_reject)) {
        to_remove_sam <- c(to_remove_sam,i)
        next()
      } else {
        # give warning if the sample misses >0.5 of the rules
        if (sum(is_na) > (nrow(binary)*0.5) & verbose) {
          message("More than the half of the rules need imputation for this sample:")
          message(colnames(binary)[i])
          message("This could affect the prediction accuracy for this sample!")
        }
      }

      # remove the NAs before find the dist
      ok_rules <- names(which(is_na==FALSE))
      impute_rules <- names(which(is_na == TRUE)) #Pontus

      dist_mat <- rdist::cdist(TrainingMatrix_df[, ok_rules],
                               t(binary[ok_rules, i, drop=FALSE]),
                               metric="jaccard") #Pontus (Uses the package "rdist")

      binary[impute_rules, i] <- apply(TrainingMatrix_df[head(order(dist_mat[,1]),
                                                              impute_kNN),
                                                         impute_rules,
                                                         drop=FALSE],
                                       2, getmode)

      # sam      <- binary[ok_rules,i, drop=FALSE]

      # dist_mat <- as.matrix(dist(t(cbind(sam, mode_df[ok_rules,])),
      #                            method="binary"))
      #
      # # remove the first because it is the sample itself
      # closest  <- which.min(dist_mat[-1,1])
      # closest  <- names(closest)[1]
      #
      # # get the rules those need imputation for this sample
      # impute_rules <- names(which(is_na==TRUE))
      #
      # # get the mode values as imputations
      # binary[impute_rules, i] <- mode_df[impute_rules, closest]
    }

    # tell the user that we skipped these samples
    if (length(to_remove_sam)>0) {
      message("#####")
      message(length(to_remove_sam),
              " sample(s) with missed values and passed the impute_reject cutoff, because of that these sample(s) were removed from the prediction:")
      message(paste0(colnames(binary)[to_remove_sam],
                     collapse = " "))
      message("#####")

      binary  <- binary[,-to_remove_sam, drop=FALSE]
    }

    if (any(dim(binary)== 0)) {
      stop("No samples left!")
    }
  }

  # predict by original ranger function
  results <- predict(classifier[[1]]$RF_classifier,
                     data = t(binary))

  # give the prediction the sample names
  if (is.matrix(results$predictions)) {
    rownames(results$predictions) <- colnames(binary)

    # get the highest score
    pred <- as.data.frame(results$predictions, stringsAsFactors = FALSE)

    # get the prediction labels
    results$predictions_classes <- colnames(pred)[max.col(pred,
                                                          ties.method = "first")]

    names(results$predictions_classes) <- rownames(results$predictions)

    # to generate warnings if there is ties
    first <- colnames(pred)[max.col(pred,
                                    ties.method = "first")]
    last  <- colnames(pred)[max.col(pred,
                                    ties.method = "last")]
    if (sum(first != last)>0) {
      message(paste("Score ties were found in", sum(first != last),
                    "out of",nrow(pred),"samples in the data",
                    collapse = " "))

    }
  }

  if (is.factor(results$predictions)) {
    names(results$predictions) <- colnames(binary)
  }
  #
  return(results)
}


#' Print method for the RF classifier
#' @export
print.rule_based_RandomForest <- function(x, ...) {
  # print info about the input data and labels
  cat("multiclassPairs - Rule based Random Forest classifier\n")

  # print what is not empty in the other slots
  for (y in names(x)) {
    if (all(sapply(x[[y]], is.null))) {
      cat("  Object contains: NULL\n")
    } else {
      cat("  Object contains:\n")
      for (i in names(x[[y]])) {
        if (!is.null(x[[y]][[i]])) {
          cat("     ","-",i)

          if (i == "genes") {
            cat(":", length(x[[y]][[i]]),"genes\n")
          }

          if (i == "rules") {
            cat(":", nrow(x[[y]][[i]]),"rules\n")
          }

          if (i == "calls") {
            cat(": ")
            cat(gsub(capture.output(cat(capture.output(x[[y]][[i]]))),
                     pattern = paste0(c(",", ",     "), collapse = "|"),
                     replacement = ",\n            "))
          }

          if (i == "RF_classifier"| i== "boruta"  | i=="TrainingMatrix") {
            #| i== "mode"
            cat("\n")
          }

        } else { # if NULL
          cat("     ","-",i,": NULL\n")
        }
        message()
      }
    }
  }
}



#' Merge prediciton score amtrices from two classifiers
#'
#' This function merges the prediction score matrices from the 5-class and 7-class (UroA,UroB,UroC) classifiers into 1 unique score matrix
#' @param score_matrix1 prediction score matrix from the 7 class classifier
#' @param score_matrix2 prediction score matrix from the 5 class classifier
#' @param row.names rownames (sample names) for both matrices
#' @return Merged matrix including scores for 8 classes (Uro, UroA, UroB, UroC, GU, BaSq, Mes, ScNE)
#'
#' @export
merge_subUro_matrix <- function(score_matrix1, # Score matrix from the 7-class (Uro) classifier
                                score_matrix2, # Score matrix from the 5-class classifier
                                row.names = list(NULL,NULL)
)
{
  # Check inputs
  # Matrices
  if (!(class(score_matrix1)[1] %in%  c("matrix","data.frame") & class(score_matrix2) %in%  c("matrix","data.frame"))[1]) {
    stop("Inputs should be score matrices from a rule_based_RandomForest object")
  }

  # Row names
  if (!is.null(row.names) & length(row.names) == 2) {
    rownames(score_matrix1) <- row.names[[1]]
    rownames(score_matrix2) <- row.names[[2]]
    name_rows <- rownames(score_matrix2)
  } else {
    stop("Row names for both matrices are needed.")
  }

  # Create empty matrix
  # Number and name of rows
  n_rows <- nrow(score_matrix2)
  final_classes=c("Uro","UroA","UroB","UroC","GU","BaSq","Mes","ScNE")
  final_matrix <- matrix(nrow = n_rows, ncol = length(final_classes), dimnames = list(name_rows, final_classes))

  for (i in colnames(final_matrix)) {
    if (i %in% colnames(score_matrix1) & i  %in%  colnames(score_matrix2)) {
      message("Classes should be just on one of the provided matrices. Values for duplicate classes will be taken from matrix 1.")
    }
    else if (i  %in%  colnames(score_matrix2)) {
      final_matrix[rownames(score_matrix2),i] <- score_matrix2[,i]
    } else if (i  %in%  colnames(score_matrix1)) {
      final_matrix[rownames(score_matrix1),i] <- score_matrix1[,i]
    }
  }
  return(final_matrix)
}

#' Predict Lund Taxonomy subtypes based on rule-based Random Forest classifiers
#'
#' @param data matrix, data frame or multiclassPairs_object of gene expression values
#' @param include_data include data in output (disabled by default)
#' @param include_scores include prediciton scores for each sample and class in output (default)
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
#'
#' @examples data(Lund2017)
#' predict_LundTax2023(Lund2017)
#' @examples
#' # Include data in result
#' data(Lund2017)
#' predict_LundTax2023(Lund2017,
#' include_data = TRUE)
#'
#' @examples
#' # Imputation
#' data(Lund2017)
#' predict_LundTax2023(Lund2017,
#' impute = TRUE)
#'
#' @export
#

predict_LundTax2023 <- function(data,
                                include_data=FALSE, # return input data in the results object
                                include_scores=TRUE, # return prediction scores in the results object
                                ...)

{
  # Check inputs

  ## Data ##
  # Store as dataframe
  if (!(class(data)[1] %in% c("matrix","data.frame","multiclassPairs_object"))) {
    stop("Data should be in one of the following formats: matrix, data.frame, multiclassPairs_object")
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

  # Classifier ##
  C <- LundTax2023::LundTax_RF_5c
  C2 <- LundTax2023::LundTax_RF_Uro7c

  # # testing
  # C <- RF_classifier_IHHK_500_500_10_d5
  # C2 <- RF_classifier_URO_IHHK_500_500_10_d5_final


  # new results object

  results_suburo <- list(data = D,
                         scores = NULL,
                         predictions_7classes = NULL,
                         predictions_5classes = NULL)

  prediction <- predict_RF(classifier = C,
                           Data = D,
                           verbose = TRUE, ...)

  ## get uro samples
  if ("Uro" %in% prediction$predictions_classes) {
    D_Uro <- D[,which(prediction$predictions_classes == "Uro")]
    D_NoUro <- D[,which(prediction$predictions_classes != "Uro")]

    # classify suburo
    prediction_suburo <- predict_RF(classifier = C2,
                                    Data = D_Uro,
                                    verbose = TRUE, ...)

    names_uro <- colnames(D_Uro)
    names_all <- colnames(D)

    # Merged score matrix
    score_matrix <- merge_subUro_matrix(score_matrix1 = prediction_suburo$predictions,
                                        score_matrix2 = prediction$predictions,
                                        row.names = list(names_uro,names_all))
  } else {
    score_matrix <- prediction$predictions
  }

  results_suburo$scores <- score_matrix

  score_matrix_suburo <- score_matrix[,-1]
  score_matrix_5c <- score_matrix[,c(1,5:8)]

  results_suburo$predictions_7classes <- colnames(score_matrix_suburo)[max.col(replace(score_matrix_suburo,is.na(score_matrix_suburo),-Inf))]
  results_suburo$predictions_5classes <- colnames(score_matrix_5c)[max.col(replace(score_matrix_5c,is.na(score_matrix_5c),-Inf))]

  names(results_suburo$predictions_7classes) <- colnames(D)
  names(results_suburo$predictions_5classes) <- colnames(D)

  predictions_suburo <- list(predictions_7classes = results_suburo$predictions_7classes,
                             predictions_5classes = results_suburo$predictions_5classes)

  results_suburo_nodata <- list(scores = results_suburo$scores,
                                predictions_7classes = results_suburo$predictions_7classes,
                                predictions_5classes = results_suburo$predictions_5classes)


  if (include_data & include_scores) {
    result <- results_suburo
  } else if (include_data == FALSE & include_scores) {
    result <- results_suburo_nodata
  } else if (include_data == FALSE & include_scores == FALSE) {
    result <- predictions_suburo
  }

}


