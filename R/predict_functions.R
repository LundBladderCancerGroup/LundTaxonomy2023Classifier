
#' Merge prediction score matrices from two classifiers
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


#' Report score ties
#'
#' Checks score ties in the prediction scores and prints a message indicating
#' the sample where the tie occurred, the two subtypes with the tied scores and the subtype that is reported in the output object
#' @param first predictions when setting ties.method = "first"
#' @param last predictions when setting ties.method = "last"
#' @return Message indicating sample with score tie and the highest scored and reported subtypes
#'
#' @export
check.ties <- function(first, last) {
  for (tie in which(first != last)) {
    sample <- names(first)[tie]
    subtype1 <- first[tie]
    subtype2 <- last[tie]
    message(paste0("Sample ",sample, ": tie between ", subtype1, " and ", subtype2, "\nOutput prediction is ", subtype1))
  }
}



#' Calculate proliferation and progression scores
#'
#' @param Data matrix or data frame of gene expression values
#' @param logTransform if TRUE, log transform data. Set to FALSE is data is already log transformed
#' @param gene_id specify the type of gene identifier used in the data:
#' - "hgnc_symbol" for HUGO gene symbols
#' - "ensembl_gene_id" for Ensembl gene IDs
#' - "entrezgene" for Entrez IDs
#' Default value is hgnc_symbol
#' @param variable score to calculate: proliferation or progression
#' @return
#' Returns values of the calculated score for each sample.
#'
#'
#' @examples
#' # Calculate proliferation score, hgnc symbols 
#' results <- ratio_score(Lund2017, variable = "proliferation", gene_id = "hgnc_symbol")
#' 
#' @examples
#' # Calculate progression score, hgnc symbols 
#' results <- ratio_score(Lund2017, variable = "progression", gene_id = "hgnc_symbol")
#'                                          impute = TRUE)
#'
#' @export
#
ratio_score <- function(Data,
                        variable = c("proliferation", "progression", NULL)[3],
                        logTransform = FALSE,
                        gene_id = c("ensembl_gene_id", "hgnc_symbol")[1],
                        method = c("ratio","singscore")[1]
)
{
  # Data must be a matrix in log2 transformed format, with sample as column and genes as rows
  if (!class(Data)[1] %in% c("data.frame","matrix")) {
    stop("Data must be in dataframe or matrix format.")
  }
  D <- Data
  
  if (logTransform) {
    D <- log2(D+1)
  }
  
  if (is.null(variable)) {
    stop("Variable must be one of the following: proliferation, progression")
  }
  
  
  # Load signatures
  
  signatures <- LundTax2023Classifier::signatures
  
  if (variable == "proliferation") {
    proliferation_signature <- signatures$proliferation
    # When included in the package
    # proliferation_signature <- LundTax2023Classifier::signatures$proliferation
    
    up_genes <- proliferation_signature[proliferation_signature$signature == "LateCellCycle",gene_id]
    down_genes <- proliferation_signature[proliferation_signature$signature == "EarlyCellCycle",gene_id]
    
    
    int_up_genes <- intersect(rownames(D), up_genes)
    int_down_genes <- intersect(rownames(D), down_genes)
    
    diff_genes <- length(c(up_genes,down_genes)) - length(c(int_up_genes,int_down_genes))
    
    if(diff_genes > 0) {
      message(paste0("Proliferation score: ", diff_genes, "/",length(c(up_genes,down_genes)), " genes are missing from the data."))
      print(setdiff(c(up_genes,down_genes), c(int_up_genes,int_down_genes)))
    }
    
    up_genes <- int_up_genes
    down_genes <- int_down_genes
  }
  
  if (variable == "progression") {
    
    progression_signature <- signatures$progression
    
    # When included in the package
    # progression_signature <- LundTax2023Classifier::signatures$progression
    
    up_genes <- progression_signature[progression_signature$direction_in_prog == "Up",gene_id]
    down_genes <- progression_signature[progression_signature$direction_in_prog == "Down",gene_id]
    
    int_up_genes <- intersect(rownames(D), up_genes)
    int_down_genes <- intersect(rownames(D), down_genes)
    
    diff_genes <- length(c(up_genes,down_genes)) - length(c(int_up_genes,int_down_genes))
    
    if (diff_genes > 0) {
      message(paste0("Progression score: ",diff_genes, "/",length(c(up_genes,down_genes)), " genes are missing from the data."))
      print(setdiff(c(up_genes,down_genes), c(int_up_genes,int_down_genes)))
    }
    
    up_genes <- int_up_genes
    down_genes <- int_down_genes
    
  }
  
  if (method == "ratio") {
    
    # Ratio of ranks #
    rank_data <- apply(D,2,rank)
    rank_data <- rank_data/nrow(rank_data)
    
    median_UPgenes <- apply(rank_data[up_genes,,drop=F],2,median)
    median_DOWNgenes <- apply(rank_data[down_genes,,drop=F],2,median)
    
    median_UP_DOWN <- median_UPgenes/median_DOWNgenes
    
    r_score <- data.frame(Score=median_UP_DOWN,row.names = colnames(D))
    
  } else if (method == "singscore") { # I was testing both methods but I think we decided to keep only the ratio to avoid using an extra package
    # So this part could be removed
    
    # Singscore #
    require(singscore)
    rankData <-  rankGenes(D)
    r_score <- simpleScore(rankData, upSet = up_genes, downSet = down_genes, centerScore = FALSE)
    
  }
  
  return(r_score)
  
}


#' Calculate immune, infiltration (141 UP) and prostate scores
#'
#' @param Data matrix or data frame of gene expression values
#' @param variable score to calculate: proliferation or progression
#' @param logTransform if TRUE, log transform data. Set to FALSE is data is already log transformed
#' @param gene_id specify the type of gene identifier used in the data:
#' - "hgnc_symbol" for HUGO gene symbols
#' - "ensembl_gene_id" for Ensembl gene IDs
#' - "entrezgene" for Entrez IDs
#' Default value is hgnc_symbol
#' @return
#' Returns values of the calculated score for each sample.
#'
#'
#' @examples
#' # Calculate immune score, hgnc symbols 
#' results <- single_score(Lund2017, variable = "immune", gene_id = "hgnc_symbol")
#' 
#' @examples
#' # Calculate progression score, hgnc symbols 
#' results <- single_score(Lund2017, variable = "prostate", gene_id = "hgnc_symbol")
#'                                          impute = TRUE)
#'
#' @export
#
single_score <- function(Data,
                         variable = c("immune", "score141up", "prostate", NULL)[4],
                         logTransform = FALSE,
                         gene_id = c("ensembl_gene_id", "hgnc_symbol")[1],
                         adjust = TRUE,
                         adj_factor = 5.1431
) 
{
  # Data must be a matrix in log2 transformed format, with sample as column and genes as rows
  if (!class(Data)[1] %in% c("data.frame","matrix")) {
    stop("Data must be in dataframe or matrix format.")
  }
  D <- Data
  
  if (logTransform) {
    D <- log2(D+1)
  }
  
  if (is.null(variable)) {
    stop("Variable must be one of the following: immune, score141up, prostate")
  }
  
  if (!(gene_id %in% c("hgnc_symbol","ensembl_gene_id"))) {
    stop("Gene ID must be one of the following: 'hgnc_symbol' or 'ensembl_gene_id'")
  }
  else if (gene_id != "hgnc_symbol") {
    
    # # Testing
    # load("D:/Packages/LundTaxonomy2023Classifier_DEV/gene_info_lund.rda")
    # load("C:/Users/earam/LBCG/gene_info_lund.rda")
    gene_info_lund <- LundTax2023Classifier::gene_info_lund
    
    rownames(gene_info_lund) <- gene_info_lund[[gene_id]]
    int_genes <- rownames(D)[which(rownames(D) %in% gene_info_lund[[gene_id]])]
    rownames(D)[which(rownames(D) %in% gene_info_lund[[gene_id]])] <- gene_info_lund[int_genes,"hgnc_symbol"]
    
  }
  gene_id <- "hgnc_symbol"
  
  # Load signatures
  
  # Testing 
  # load("C:/Users/earam/LBCG/signatures.rda")
  
  signatures <- LundTax2023Classifier::signatures
  
  if (variable == "immune") {
    
    s <- signatures$immune[,c(gene_id, "signature"), drop = F]
    
    # When included in the package
    # LundImmune <- LundTax2023Classifier::signatures$immune
    
    genes_immune <- unique(c(s[[gene_id]]))
    genes_immune_int <- intersect(rownames(D),genes_immune)
    
    diff_genes <- length(genes_immune) - length(genes_immune_int)
    
    if (diff_genes > 0) {
      message(paste0("Immune scores: ", diff_genes, "/",length(genes_immune)," genes are missing from the data."))
      print(setdiff(genes_immune, genes_immune_int))
    }
    
    
  }
  
  if (variable == "score141up") {
    
    
    Immune141_UP <- signatures$signatures_plot[which(signatures$signatures_plot$signature == "Immune141_UP."),,drop = F]
    Stromal141_UP <- signatures$signatures_plot[which(signatures$signatures_plot$signature == "Stromal141_UP."),,drop = F]
    
    colnames(Immune141_UP)[2] <- "hgnc_symbol"
    colnames(Stromal141_UP)[2] <- "hgnc_symbol"
    
    genes141 <- unique(c(Immune141_UP$hgnc_symbol, Stromal141_UP$hgnc_symbol))
    genes141_int <- intersect(rownames(D),genes141)
    
    diff_genes <- length(genes141) - length(genes141_int)
    
    if (diff_genes > 0) {
      message(paste0("Infiltration scores: ", diff_genes, "/",length(genes141)," genes are missing from the data."))
      print(setdiff(genes141,genes141_int))
    }
    
    # Updated package
    # signatures_plot <- LundTax2023Classifier::signatures$signatures_plot
    # Immune141_UP <- signatures_plot[which(signatures_plot$signature == "Immune141_UP."),]
    # Stromal141_UP <- signatures_plot[which(signatures_plot$signature == "Stromal141_UP."),]
    
    s <- rbind(Immune141_UP,Stromal141_UP)
    
    
  }
  
  if (variable == "prostate") {
    # score_results <- as.data.frame(matrix(nrow = ncol(Data),
    #                                       ncol = 1,
    #                                       dimnames = list(colnames(Data), "ProstateScore")))
    
    s <- signatures$prostate[,gene_id, drop = F]
    s$signature <- "Prostate"
    
    genes_prostate_int <- intersect(rownames(Data),s[[gene_id]])
    
    diff_genes <- length(s[[gene_id]]) - length(genes_prostate_int)
    
    if(diff_genes > 0) {
      message(paste0("Prostate score: ", diff_genes, "/",length(s[[gene_id]])," genes are missing from the data."))
      print(setdiff(s[[gene_id]],genes_prostate_int))
    }
    
  }
  
  # Results object
  score_results <- as.data.frame(matrix(nrow = ncol(Data),
                                        ncol = length(unique(s$signature)),
                                        dimnames = list(colnames(Data), unique(s$signature))))
  
  for (i in unique(s$signature)) {
    genes_signature_int <- intersect(rownames(D),s[s$signature == i, gene_id])
    res <- unlist(lapply(1:ncol(D),function(x) {mean(D[genes_signature_int,x])}))
    score_results[,i] <- res
  }
  
  if (adjust) {
    
    StableGenes <- signatures$stable_genes
    # StableGenes <- LundTax2023Classifier::signatures$stable_genes
    
    stable_genes_int <- intersect(rownames(D),StableGenes[,gene_id])
    
    
    score_results <- do.call("rbind",lapply(1:nrow(score_results),function(x){
      (score_results[x,]/mean(D[stable_genes_int,x]))*adj_factor
    }))
  }
  
  return(score_results)
  
}

# Immune Proportions #####

calculate_immune_proportions <- function(immune_results) {
  immune_proportions <- t(apply(immune_results,1,function(x){x/sum(x)}))
  colnames(immune_proportions) <- paste0(colnames(immune_proportions)," Proportion")
  return(immune_proportions)
}


# Grade predictor ##########

predict_grade <- function(Data, grade_predictor,
                          gene_id = c("ensembl_gene_id", "hgnc_symbol")[1],
                          ...) {
  
  # Data must be a matrix in log2 transformed format, with sample as column and genes as rows
  if (!class(Data)[1] %in% c("data.frame","matrix")) {
    stop("Data must be in dataframe or matrix format.")
  }
  D <- Data
  
  if (class(grade_predictor)[1] != "rule_based_RandomForest") {
    stop("Classifier must be a rule_based_RandomForest object.")
  }
  
  if (!(gene_id %in% c("hgnc_symbol","ensembl_gene_id"))) {
    stop("Gene ID must be one of the following: 'hgnc_symbol' or 'ensembl_gene_id'")
  } else if (gene_id != "ensembl_gene_id") {
    
    # Testing
    # load("D:/UROSCANSEQ_2024/Analysis/02.New_data/GradeClassifier/gene_info_grade_classifiers.RData")
    # load("C:/Users/earam/LBCG/gene_info_lund.rda")
    gene_info_lund <- LundTax2023Classifier::gene_info_lund
    rownames(gene_info_lund) <- gene_info_lund[[gene_id]]
    int_genes <- rownames(D)[which(rownames(D) %in% gene_info_lund[[gene_id]])]
    rownames(D)[which(rownames(D) %in% gene_info_lund[[gene_id]])] <- gene_info_lund[int_genes,"ensembl_gene_id"]
    
  }
  
  
  require(multiclassPairs)
  grade_results <- predict_RF(classifier = grade_predictor, Data = D, ...)
  
  
  
  # if (ncol(D) == 1) cat(paste0("Prediction: ", grade_results$predictions_classes, "\n","Score: ", grade_results$predictions[,2], "\n"))
  
  return(grade_results)
}


# Function to calculate all scores #

lund_scores <- function(Data, # Input data. Data must be a matrix in log2 transformed format, with sample as column and genes as rows
                        gene_id = c("hgnc_symbol", "ensembl_gene_id")[2], # gene IDs
                        scoring_method = c("ratio", "singscore")[1], # Method to calculate the proliferation score
                        threshold_prostate = 4, # Gene expression threshold to flag a sample as possible prostate
                        threshold_progression = 0.58, #  threshold to flag a sample as high risk of progression
                        logTransform = TRUE, # Scores are calculated on log transformed data. If the data is already log transformed, set logTransformed to FALSE. If logTranform = TRUE, data will be log2 transformed (log2(data+1)) before calculating the scores
                        adjust = FALSE,
                        adj_factor = 5.1431,
                        verbose = FALSE,
                        ... # arguments to pass to the ranger functions (add impute = TRUE here if gene are missing)
)
{
  require(multiclassPairs)
  # Check data
  # Data must be a matrix in log2 transformed format, with sample as column and genes as rows
  if (!class(Data)[1] %in% c("data.frame","matrix")) {
    stop("Data must be in dataframe or matrix format.")
  }
  
  D <- Data
  # Gene symbols
  if (!gene_id %in% c("hgnc_symbol", "ensembl_gene_id")) {
    stop("gene_id must be one of: 'hgnc_symbol', 'ensembl_gene_id'")
  }
  
  # Log transform #
  if (logTransform) {
    D <- log2(D+1)
  }
  
  # Apply #
  
  # Proliferation #
  # Debug
  if (verbose) print("Calculating Proliferation")
  results_proliferation <- ratio_score(Data = D,
                                       variable = "proliferation",
                                       gene_id = gene_id,
                                       method = scoring_method)
  # Grade #
  # WHO 1999 (G3 vs G1/2)
  classifier_GRADE3 <- LundTax2023Classifier::classifier_GRADE3
  
  # WHO 2004/2016 (HG vs LG)
  classifier_HG <- LundTax2023Classifier::classifier_HG
  
  if (verbose) print("Calculating G3 score")
  results_g3 <- predict_grade(Data = D,
                              gene_id = gene_id,
                              grade_predictor = classifier_GRADE3,
                              ...)
  
  if (verbose) print("Calculating HG score")
  results_hg <- predict_grade(Data = D,
                              gene_id = gene_id,
                              grade_predictor = classifier_HG,
                              ...)
  
  # Progression #
  if (verbose) print("Calculating Progression score")
  score_progression <- ratio_score(Data = D,
                                     variable = "progression",
                                     gene_id = gene_id,
                                     method = scoring_method)
  results_progression <- ifelse(score_progression$Score >= threshold_progression, "HR", "LR")
  
  # Prostate #
  
  if (verbose) print("Calculating Prostate score")
  scores_prostate <- single_score(Data = D,
                                  variable = "prostate",
                                  logTransform = FALSE,
                                  gene_id = gene_id,
                                  adjust = adjust,
                                  adj_factor = adj_factor)
  
  results_prostate <- as.numeric(scores_prostate >= threshold_prostate)
  
  # Immune
  
  if (verbose) print("Calculating Immune scores")
  results_immune <- single_score(Data = D,
                                 variable = "immune",
                                 logTransform = FALSE,
                                 gene_id = gene_id,
                                 adjust = adjust,
                                 adj_factor = adj_factor)
  
  # 141 UP
  if (verbose) print("Calculating 141UP scores")
  scores141up <- single_score(Data = D,
                              variable = "score141up",
                              logTransform = FALSE,
                              gene_id = gene_id,
                              adjust = adjust,
                              adj_factor = adj_factor)
                            
  
  # Immune Proportion
  immune_proportions <- calculate_immune_proportions(results_immune)
  
  # Add 141UP to immune results
  results_immune <- cbind(Immune141_UP=scores141up$Immune141_UP,
                          results_immune[,1:10,drop=FALSE],
                          Stromal141_UP=scores141up$Stromal141_UP,
                          results_immune[,11:13,drop=FALSE])
  
  
  
  # Merge_scores
  
  if (verbose) print("Merging scores")
  merge_scores <- cbind(ProliferationScore=results_proliferation$Score,
                        MolecularGradeWHO1999=results_g3$predictions_classes,
                        MolecularGradeWHO1999_score=results_g3$predictions[,"G3"],
                        MolecularGradeWHO2016=results_hg$predictions_classes,
                        MolecularGradeWHO2016_score=results_hg$predictions[,"HG"],
                        ProgressionScore = score_progression$Score,
                        ProgressionRisk = results_progression,
                        ProstateScore = scores_prostate,
                        PossibleProstate = ifelse(results_prostate == 1, "YES", "NO"),
                        results_immune,
                        immune_proportions
  )
  
  
  
}




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
#'
#' @examples
#' results <- predict_LundTax2023(Lund2017)
#' @examples
#' # Include data in result
#' results_data <- predict_LundTax2023(Lund2017,
#'                                     impute = TRUE,
#'                                     include_data = TRUE)
#'
#' @export
#
predict_LundTax2023 <- function(data,
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
  all_scores <- lund_scores(Data = original_D,
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



