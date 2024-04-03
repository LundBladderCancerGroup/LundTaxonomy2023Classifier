
#' Classifier as a 'rule_based_RandomForest' object. Predicts samples as one of
#' the 5 main Lund Taxonomy molecular subtypes, Uro, GU, BaSq, Mes, or ScNE
#' Object includes the final RF classifier, the used genes and rules in the final model,
#' the Boruta results, and the training matrix. The training matrix is a binary matrix
#' containing the rule values for the training data and it is used for
#' imputation purposes during the prediction if values are missing in the sample.
#' This object was generated using the multiclassPairs package
#'
#' @references \url{https://github.com/NourMarzouka/multiclassPairs}
"LundTax_RF_5c"

#' Classifier as a 'rule_based_RandomForest' object. Predicts samples as one of
#' the 3 Uro subclasses, UroA, UroB, or UroC
#' Object includes the final RF classifier, the used genes and rules in the final model,
#' the Boruta results, and the training matrix. The training matrix is a binary matrix
#' containing the rule values for the training data and it is used for
#' imputation purposes during the prediction if values are missing in the sample.
#' This object was generated using the multiclassPairs package
#'
#' @references \url{https://github.com/NourMarzouka/multiclassPairs}
"LundTax_RF_Uro7c"

#' Gene expression data derived from the Sjödahl et. al. (2017) cohort
#'
#' Matrix of RMA normalized and ComBat adjusted gene expression values for 15697 genes with
#' HGNC symbols and 301 samples
#'
#' @references \url{https://onlinelibrary.wiley.com/doi/full/10.1002/path.4886}
"Lund2017"

#' Colors used for Lund Taxonomy subtypes
"lund_colors"

#' Signature table

"signatures"

#' Gene IDs for classification
"gene_info_classifier"

#' Gene IDs for plotting
"gene_info_heatmap"

#' LundTax2023
#'
#' This packages implements a Random Forest rule-based single-sample predictor that classifies
#' transcriptomic samples into the 5 (or 7, including subclasses) Lund Taxonomy molecular subtypes.
#' The final classifier is composed of two separate predictors applied sequentially: first a sample
#' is classified as one of the 5 main classes (Uro, GU, BaSq, Mes or ScNE), and then samples classified as
#' Uro are subclassified into UroA, UroB or UroC by a second predictor
"_PACKAGE"
