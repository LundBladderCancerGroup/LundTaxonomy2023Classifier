#' Plot heatmap for classification results
#'
#' This function plots a heatmap including genes and signatures of interest, with prediction results and scores on top
#' @param results_object description
#' @param data results_object
#' @param title heatmap title
#' @param annotation use 5 or 7 class annotation
#' @param plot_scores plot prediction scores for each class
#' @param show_ann_legend show annotation legend (Lund classes)
#' @param show_hm_legend show heatmap legend
#' @param set_order set sample order. Samples are split by subtype, and order within each subtype. By default, samples are order by late/early cell cycle ratio (low to high)
#' @param ann_heigh annotation height in cm (default = 6)
#' @param font.size font size (default = 8)
#' @param norm normalize the data into Z-scaled values (default TRUE)
#' @param gene_id specify the type of gene identifier used in the data:
#' - "hgnc_symbol" for HUGO gene symbols
#' - "ensembl_gene_id" for Ensembl gene IDs
#' - "entrezgene" for Entrez IDs
#' Default value is hgnc_symbol
#' @return Draws heatmap and silently returns the sample order
#'
#'@examples
#'# Including data in results object
#' results <- predict_LundTax2023(Lund2017,
#'                               include_data = TRUE)
#' plot_signatures(results)
#'
#'@examples
#'# Loading data separately
#'results <- predict_LundTax2023(Lund2017)
#'plot_signatures(results, data = Lund2017)
#'
#'# 5 class annotation and Without prediction scores
#' plot_signatures(results,
#'                data = Lund2017,
#'                annotation = "5 classes",
#'                ann_height = 0.5,
#'                plot_scores = FALSE)
#'
#'# Plot and get sample order
#' sample_order <- plot_signatures(results,
#'                                data = Lund2017)
#'print(head(sample_order))
#'
#'# Save to pdf
#'pdf("heatmap_example.pdf", width = 15, height = 10)
#'plot_signatures(results,
#'                data = Lund2017,
#'                ann_height = 0.5,
#'                plot_scores = FALSE)
#'dev.off()
#'
#'@export


plot_signatures <- function(results_object,
                            data = NULL,
                            title = "",
                            gene_id = c("hgnc_symbol","ensembl_gene_id","entrezgene")[1],
                            annotation = c("5 classes", "7 classes")[2],
                            plot_scores = TRUE,
                            show_ann_legend = FALSE,
                            show_hm_legend = FALSE,
                            set_order = NULL,
                            ann_height = 6,
                            font.size = 8,
                            norm = TRUE
) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) | !requireNamespace("circlize", quietly = TRUE)) {
    stop("The ComplexHeatmap and circlize packages must be installed to use this function")
  }

  suppressPackageStartupMessages(require("ComplexHeatmap",quietly = TRUE))


  ### Data ######

  if (is.null(results_object) | !is.list(results_object)) {
    stop("Input should be the result of applying predict_LundTax2023")
  } else if (is.list(results_object) & "data" %in% names(results_object)) {
    D <- results_object$data
    score_matrix <- results_object$subtype_scores
    pred_labels5 <- results_object$predictions_5classes
    pred_labels7 <- results_object$predictions_7classes
  } else if (is.list(results_object) & !"data" %in% names(results_object)) {
    if (is.null(data)) {
      stop("Data is missing. Include it in results object by running predict_LundTax2023 with include_data = TRUE or provide it in the data argument.")
    } else if (!(class(data)[1] %in%  c("matrix","data.frame"))) {
      stop("Data should be in matrix or data.frame format")
    } else if (class(data)[1] %in%  c("matrix","data.frame")) {
      D <- data[,names(results_object$predictions_7classes)]
      score_matrix <- results_object$subtype_scores
      pred_labels5 <- results_object$predictions_5classes
      pred_labels7 <- results_object$predictions_7classes
    }
  }

  ## Legend ##
  if (!is.null(show_ann_legend)) {
    show_legend = show_ann_legend
  }

  ## Scale ##
  if (norm) {

    D_norm<-as.matrix(D)
    D_norm <- scale(t(D))
    D_norm <- t(D_norm)

  } else {
    D_norm <- as.matrix(D)
  }

  ######## Heatmaps ########

  # Gene signatures #
  signatures <- LundTax2023Classifier::signatures
  signatures_plot <- signatures$signatures_plot

  genes_to_plot <- list(Early_CC=c(signatures$proliferation[which(signatures$proliferation$signature == "EarlyCellCycle"),2]),
                        Late_CC=c(signatures$proliferation[which(signatures$proliferation$signature == "LateCellCycle"),2]),
                        Late_Early=NULL,
                        UroDiff=c("PPARG","FOXA1","GATA3","ELF3"),
                        UPKs=c("UPK1A","UPK1B","UPK2","UPK3A","KRT20"),
                        Circuit=c("FGFR3","CCND1","E2F3","RB1","CDKN2A"),
                        Circuit_score=NULL,
                        FGFR3=c(signatures_plot[which(signatures_plot$signature == "FGFR3."),2]),
                        BaSq=c("KRT5","KRT14","FOXA1","GATA3"),
                        BaSq_ratio=NULL,
                        Keratinization=c(signatures_plot[which(signatures_plot$signature == "Keratinization_QTC."),2]),
                        Adhesion=c("EPCAM","CDH1","CDH3"),
                        MYC=c("MYCL","MYCN","MYC"),
                        ERBB=c("EGFR","ERBB2","ERBB3"),
                        ERBB_score=NULL,
                        ScNE=c("CHGA","SYP","ENO2"),
                        Immune141_UP=c(signatures_plot[which(signatures_plot$signature == "Immune141_UP."),2]),
                        Stromal141_UP=c(signatures_plot[which(signatures_plot$signature == "Stromal141_UP."),2]),
                        Immune141_UP_score=NULL,
                        Stromal141_UP_score=NULL)
  # Lund colors #
  lund_colors <- LundTax2023Classifier::lund_colors
  # Gene IDs #

  if (!(gene_id %in% c("hgnc_symbol","ensembl_gene_id","entrezgene"))) {
    stop("Gene ID must be one of the following: 'hgnc_symbol','ensembl_gene_id' or 'entrezgene'")
  }
  else if (gene_id != "hgnc_symbol") {
   
    gene_info_lund <- LundTax2023Classifier::gene_info_lund
    rownames(gene_info_lund) <- gene_info_lund[[gene_id]]
    int_genes <- rownames(D_norm)[which(rownames(D_norm) %in% gene_info_lund[[gene_id]])]
    rownames(D_norm)[which(rownames(D_norm) %in% gene_info_lund[[gene_id]])] <- gene_info_lund[int_genes,"hgnc_symbol"]

  }

  ## hm1 -> late/early cell cycle #####

  genes_early <- genes_to_plot$Early_CC[which(genes_to_plot$Early_CC %in% rownames(D_norm))]
  genes_late <- genes_to_plot$Late_CC[which(genes_to_plot$Late_CC %in% rownames(D_norm))][1:10]
  genes_cc <- na.omit(c(genes_early,genes_late))

  # Check if genes from both early and late are present
  if (length(genes_early) != 0 & !all(is.na(genes_late))) {
    # Row split for the heatmap
    row_split <- c(rep("Early",length(genes_early)),
                   rep("Late",length(genes_late)))
    # Row title
    row_title_cc <- c("Late Cell Cycle", "Early Cell Cycle")

    # Late and Early scores
    late_score <- apply(D_norm[intersect(rownames(D_norm),genes_to_plot$Late_CC),], 2, median)
    early_score <- apply(D_norm[intersect(rownames(D_norm),genes_to_plot$Early_CC),], 2, median)

    # Ratio
    late_early <- late_score-early_score

    genes_to_plot$Late_Early <- late_early
    col_fun_cc <- circlize::colorRamp2(c(quantile(late_early, 0.05),median(late_early),quantile(late_early, 0.95)),
                                       c("blue","white", "red"))

    # Order samples by late_early cell cycle
    sample_order <- order(late_early)

    # Heatmap annotation
    col = list(Predictions = lund_colors$lund_colors,
               late_early = col_fun_cc)
    ha1b = HeatmapAnnotation(late_early = genes_to_plot$Late_Early,
                             annotation_name_side = "left",
                             simple_anno_size = unit(4, "mm"),
                             simple_anno_size_adjust = TRUE,
                             col = col,
                             show_legend = show_ann_legend,
                             border = TRUE,
                             annotation_name_gp = gpar(fontsize = 8))


  } else if (length(genes_early) == 0 & !all(is.na(genes_late))) { # If early cell cycle is missing
    message("All genes from the early cell cycle signature are missing.\nSamples will not be ordered by cell cycle.")
    # skip_cellcycle <- TRUE
    sample_order <- NULL
    ha1b <- NULL
    row_split <- NULL
    row_title_cc <- "Late cell cycle"
  } else if (length(genes_early) != 0 & all(is.na(genes_late))) { # If late cell cycle is missing
    message("All genes from the late cell cycle signature are missing.\nSamples will not be ordered by cell cycle.")
    # skip_cellcycle <- TRUE
    sample_order <- NULL
    ha1b <- NULL
    row_split <- NULL
    row_title_cc <- "Early cell cycle"
  } else { # Both are missing
    message("All genes from the late/early cell cycle signatures are missing.\nSamples will not be ordered by cell cycle.")
    # skip_cellcycle <- TRUE
    sample_order <- NULL
    ha1b <- NULL
    row_split <- NULL
    row_title_cc <- NULL
  }

  # Subtype annotations #

  # 7 classes
  if (annotation == "7 classes") {

    # Predictions
    pred_lab <- pred_labels7

    # Column split
    split <- factor(pred_lab,levels=c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))

    # Score plots
    if (plot_scores) {
      bar0 = anno_barplot(as.numeric(score_matrix[,"Uro"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["Uro"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar1 = anno_barplot(as.numeric(score_matrix[,"UroA"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["UroA"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar2 = anno_barplot(as.numeric(score_matrix[,"UroB"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["UroB"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar3 = anno_barplot(as.numeric(score_matrix[,"UroC"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["UroC"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar4 = anno_barplot(as.numeric(score_matrix[,"GU"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["GU"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar5 = anno_barplot(as.numeric(score_matrix[,"BaSq"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["BaSq"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar6 = anno_barplot(as.numeric(score_matrix[,"Mes"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["Mes"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
      bar7 = anno_barplot(as.numeric(score_matrix[,"ScNE"]), ylim = c(0, 1), gp = gpar(fill=lund_colors$lund_colors["ScNE"], border=NA, col=NA), bar_width = 1, height = unit(6, "mm"))
    } else {
      bar0 = NULL
      bar1 = NULL
      bar2 = NULL
      bar3 = NULL
      bar4 = NULL
      bar5 = NULL
      bar6 = NULL
      bar7 = NULL
    }

    # Colors #
    col = list(Predictions = lund_colors$lund_colors)


    # Annotation first heatmap #
    ha1 = HeatmapAnnotation(Predictions = pred_lab,
                            annotation_name_side = "left",
                            col = col,
                            Uro = bar0,
                            GU = bar4, BaSq = bar5, Mes = bar6,
                            ScNE = bar7,
                            UroA = bar1, UroB = bar2, UroC = bar3,
                            na_col = "gray83",
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            show_legend = show_ann_legend,
                            border = TRUE,
                            height = unit(ann_height, "cm"),
                            annotation_name_gp = gpar(fontsize = 8))


    # 5 classes
  } else if (annotation == "5 classes") {

    pred_lab <- pred_labels5

    # Column split
    split <- factor(pred_lab,levels=c("Uro","GU","BaSq","Mes","ScNE"))


    # Score plots
    if (plot_scores) {
      bar1 = anno_barplot(as.numeric(score_matrix[,"Uro"]), ylim = c(0, 1),gp = gpar(fill=lund_colors$lund_colors["UroA"],border=NA,col=NA),bar_width = 1,height = unit(6, "mm"))
      bar2 = anno_barplot(as.numeric(score_matrix[,"GU"]), ylim = c(0, 1),gp = gpar(fill=lund_colors$lund_colors["GU"],border=NA,col=NA),bar_width = 1,height = unit(6, "mm"))
      bar3 = anno_barplot(as.numeric(score_matrix[,"BaSq"]), ylim = c(0, 1),gp = gpar(fill=lund_colors$lund_colors["BaSq"],border=NA,col=NA),bar_width = 1,height = unit(6, "mm"))
      bar4 = anno_barplot(as.numeric(score_matrix[,"Mes"]), ylim = c(0, 1),gp = gpar(fill=lund_colors$lund_colors["Mes"],border=NA,col=NA),bar_width = 1,height = unit(6, "mm"))
      bar5 = anno_barplot(as.numeric(score_matrix[,"ScNE"]), ylim = c(0, 1),gp = gpar(fill=lund_colors$lund_colors["ScNE"],border=NA,col=NA),bar_width = 1,height = unit(6, "mm"))
    } else {
      bar1 = NULL
      bar2 = NULL
      bar3 = NULL
      bar4 = NULL
      bar5 = NULL
    }

    # Colors #
    col = list(Predictions = lund_colors$lund_colors)

    # Annotation first heatmap #
    ha1 = HeatmapAnnotation(Predictions = pred_lab,
                            annotation_name_side = "left",
                            col = col,
                            Uro = bar1, GU = bar2,
                            BaSq = bar3, Mes = bar4,
                            ScNE = bar5,
                            na_col = "gray83",
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            show_legend = show_ann_legend,
                            border = TRUE,
                            height = unit(ann_height, "cm"),
                            annotation_name_gp = gpar(fontsize = 8))

  }


  # Heatmap colors
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

  # Heatmap 1
  hm1 = Heatmap(D_norm[genes_cc,,drop=FALSE],
                top_annotation = ha1,
                bottom_annotation = ha1b,
                name="hm1_cc",
                col = col_fun,
                column_split = split,
                row_split = row_split,
                row_title = row_title_cc,
                row_title_rot = 0,
                cluster_row_slices = FALSE,
                cluster_column_slices = FALSE,
                cluster_columns = FALSE,
                column_order = sample_order,
                row_names_side = "left",
                show_column_names = FALSE,
                show_row_names = FALSE,
                show_row_dend = FALSE,
                border = TRUE,
                row_names_gp = gpar(fontsize = font.size),
                row_title_gp = gpar(fontsize = 7),
                clustering_distance_rows= "spearman",
                clustering_method_rows = "ward.D2",
                column_title = title,
                height = 8,
                border_gp = gpar(lwd=0.3),
                show_heatmap_legend = FALSE,

  )

  ## hm2 -> Urodiff #########

  genes_ud <- genes_to_plot$UroDiff
  genes_ud <- genes_ud[which(genes_ud %in% rownames(D_norm))]

  hm2 = Heatmap(D_norm[genes_ud,,drop=FALSE],
                name="hm2_ud",
                col = col_fun,
                column_split = split,
                cluster_row_slices = FALSE,
                cluster_column_slices = FALSE,
                cluster_columns = FALSE,
                column_order = sample_order,
                row_names_side = "left",
                show_column_names = FALSE,
                show_row_names = TRUE,
                show_row_dend = FALSE,
                border = TRUE,
                row_names_gp = gpar(fontsize = font.size),
                row_title_gp = gpar(fontsize = 7),
                cluster_rows = FALSE,
                column_title = title,
                border_gp = gpar(lwd=0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90
  )

  ## hm3 -> UROPLAKINS #########

  genes_upk <- genes_to_plot$UPKs
  genes_upk <- genes_upk[which(genes_upk %in% rownames(D_norm))]

  hm3 = Heatmap(D_norm[genes_upk,,drop=FALSE],
                name="hm3_upk",
                col = col_fun,
                column_split = split,
                cluster_row_slices = FALSE,
                cluster_column_slices = FALSE,
                cluster_columns = FALSE,
                column_order = sample_order,
                row_names_side = "left",
                show_column_names = FALSE,
                show_row_names = TRUE,
                show_row_dend = FALSE,
                border = TRUE,
                row_names_gp = gpar(fontsize = font.size),
                row_title_gp = gpar(fontsize = 7),
                cluster_rows = FALSE,
                column_title = title,
                border_gp = gpar(lwd=0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90
  )

  ## hm4 -> Circuit score #########

  genes_circ <- genes_to_plot$Circuit
  genes_circ <- genes_circ[which(genes_circ %in% rownames(D_norm))]

  if (length(genes_circ) == 5) { # Check if all genes are present

    ### CIRCUIT SCORE
    circuit_score <- apply(D_norm, 2, function(col) sum(col[c("RB1", "FGFR3", "CCND1")]) - sum(col[c("E2F3", "CDKN2A")]))

    genes_to_plot$Circuit_score <- circuit_score

    col_fun_circ <- circlize::colorRamp2(c(quantile(circuit_score, 0.10),median(circuit_score),quantile(circuit_score, 0.90)),
                                         c("blue","white", "red"))

    col = list(circuit_score = col_fun_circ)

    ha4 = HeatmapAnnotation(circuit_score = genes_to_plot$Circuit_score,
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            annotation_name_side = "left",
                            col = col,
                            show_legend = FALSE,
                            border = TRUE,
                            annotation_name_gp = gpar(fontsize = 8))
    hm4 = Heatmap(D_norm[genes_circ,,drop=FALSE],
                  name="hm4_circ",
                  bottom_annotation = ha4,
                  col = col_fun,
                  column_split = split,
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  cluster_columns = FALSE,
                  column_order = sample_order,
                  row_names_side = "left",
                  show_column_names = FALSE,
                  show_row_names = TRUE,
                  show_row_dend = FALSE,
                  border = TRUE,
                  row_names_gp = gpar(fontsize = font.size),
                  row_title_gp = gpar(fontsize = 7),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd=0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90
    )
  } else if (length(genes_circ) != 0) { # if only some genes are missing
    # col_fun_circ <- NULL
    ha4 <- NULL
    message("Genes involved in the circuit score are missing.\nCircuit score will not be calculated.")
    hm4 = Heatmap(D_norm[genes_circ,,drop=FALSE],
                  name="hm4_circ",
                  bottom_annotation = ha4,
                  col = col_fun,
                  column_split = split,
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  cluster_columns = FALSE,
                  column_order = sample_order,
                  row_names_side = "left",
                  show_column_names = FALSE,
                  show_row_names = TRUE,
                  show_row_dend = FALSE,
                  border = TRUE,
                  row_names_gp = gpar(fontsize = font.size),
                  row_title_gp = gpar(fontsize = 7),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd=0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90
    )
  } else {
    hm4 <- NULL
  }


  ## hm 4_1 TP63 ##
  if ("TP63" %in% rownames(D_norm)) {
    tp63 <- D_norm["TP63",,drop=FALSE]
    rownames(tp63) <- c("TP63")

    hm4_1 = Heatmap(tp63,
                    name="tp63",
                    col = col_fun,
                    column_split = split,
                    cluster_row_slices = FALSE,
                    cluster_column_slices = FALSE,
                    cluster_columns = FALSE,
                    column_order = sample_order,
                    row_names_side = "left",
                    show_column_names = FALSE,
                    show_row_names = TRUE,
                    show_row_dend = FALSE,
                    border = TRUE,
                    row_names_gp = gpar(fontsize = font.size),
                    row_title_gp = gpar(fontsize = 7),
                    column_title = title,
                    border_gp = gpar(lwd=0.3),
                    show_heatmap_legend = FALSE,
                    row_title_rot = 90
    )
  } else {
    hm4_1 <- NULL
  }

  ## hm5 -> FGFR3 #########
  ## FGFR3 signature + circuit score annotation

  genes_fgfr3 <- genes_to_plot$FGFR3
  genes_fgfr3 <- genes_fgfr3[which(genes_fgfr3 %in% rownames(D_norm))]

  if (length(genes_fgfr3) != 0) {

    hm5 = Heatmap(D_norm[genes_fgfr3,,drop=FALSE],
                  name="hm5_fgfr3",
                  col = col_fun,
                  column_split = split,
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  cluster_columns = FALSE,
                  column_order = sample_order,
                  row_names_side = "left",
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  show_row_dend = FALSE,
                  row_title = "FGFR3\nsignature",
                  row_title_rot = 360,
                  border = TRUE,
                  row_names_gp = gpar(fontsize = font.size),
                  row_title_gp = gpar(fontsize = 7),
                  clustering_distance_rows= "pearson",
                  # clustering_distance_columns=  "pearson",
                  clustering_method_rows = "ward.D2",
                  # clustering_method_columns = "ward.D2",
                  column_title = title,
                  height = 4,
                  border_gp = gpar(lwd=0.3),
                  show_heatmap_legend = FALSE,

    )
  } else {
    hm5 <- NULL
  }

  ## hm6 -> BA/SQ RATIO #########
  ## baSq genes and ratio

  genes_basq <- genes_to_plot$BaSq
  genes_basq <- genes_basq[which(genes_basq %in% rownames(D_norm))]

  if (length(genes_basq) == 4) { # Check if all genes are present
    ###   BA/SQ RATIO
    basq_ratio <- apply(D_norm, 2, function(col) sum(col[c("KRT5", "KRT14", "FOXA1")]) - sum(col[c("GATA3")]))

    genes_to_plot$BaSq_ratio <- basq_ratio
    col_fun_basq <- circlize::colorRamp2(c(quantile(basq_ratio, 0.10),median(basq_ratio),quantile(basq_ratio, 0.90)),
                                         c("blue","white", "red"))

    col = list(BaSq_ratio = col_fun_basq)

    ha6 = HeatmapAnnotation(BaSq_ratio = genes_to_plot$BaSq_ratio,
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            annotation_name_side = "left",
                            col=col,
                            show_legend = FALSE,
                            border = TRUE,
                            annotation_name_gp= gpar(fontsize = 8))
    hm6 = Heatmap(D_norm[genes_basq,,drop=FALSE],
                  name="hm6_basq",
                  bottom_annotation = ha6,
                  col = col_fun,
                  column_split = split,
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  cluster_columns = FALSE,
                  column_order = sample_order,
                  row_names_side = "left",
                  show_column_names = FALSE,
                  show_row_names = TRUE,
                  show_row_dend = FALSE,
                  border = TRUE,
                  row_names_gp = gpar(fontsize = font.size),
                  row_title_gp = gpar(fontsize = 7),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd=0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90)

  } else if (length(genes_basq) != 0) { # if only some genes are missing
    ha6 <- NULL
    message("Genes involved in the BaSq score are missing.\nBaSq ratio will not be calculated.")
    hm6 = Heatmap(D_norm[genes_basq,,drop=FALSE],
                  name="hm6_basq",
                  bottom_annotation = ha6,
                  col = col_fun,
                  column_split = split,
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  cluster_columns = FALSE,
                  column_order = sample_order,
                  row_names_side = "left",
                  show_column_names = FALSE,
                  show_row_names = TRUE,
                  show_row_dend = FALSE,
                  border = TRUE,
                  row_names_gp = gpar(fontsize = font.size),
                  row_title_gp = gpar(fontsize = 7),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd=0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90)
  } else { # If all are missing
    hm6 <- NULL
    message("Genes involved in the BaSq score are missing.\nBaSq ratio will not be calculated.")
    col_fun_basq <- NULL
  }

  ## hm7 -> KERATINIZATION #########
  ## baSq genes and ratio

  genes_krt <- genes_to_plot$Keratinization
  genes_krt <- genes_krt[which(genes_krt %in% rownames(D_norm))]

  if (length(genes_krt) != 0) {
    hm7 = Heatmap(D_norm[genes_krt,,drop=FALSE],
                  name="hm7_krt",
                  col = col_fun,
                  column_split = split,
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  cluster_columns = FALSE,
                  column_order = sample_order,
                  row_names_side = "left",
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  show_row_dend = FALSE,
                  row_title = "Keratinization\nsignature",
                  row_title_rot = 0,
                  border = TRUE,
                  row_names_gp = gpar(fontsize = font.size),
                  row_title_gp = gpar(fontsize = 7),
                  clustering_distance_rows= "pearson",
                  # clustering_distance_columns=  "pearson",
                  clustering_method_rows = "ward.D2",
                  # clustering_method_columns = "ward.D2",
                  column_title = title,
                  height = 4,
                  border_gp = gpar(lwd=0.3),
                  show_heatmap_legend = FALSE,
    )
  } else {
    hm7 <- NULL
  }


  ## hm8 -> ADHESION #########

  genes_ad <- genes_to_plot$Adhesion
  genes_ad <- genes_ad[which(genes_ad %in% rownames(D_norm))]

  hm8 = Heatmap(D_norm[genes_ad,,drop=FALSE],
                name="hm8_ad",
                col = col_fun,
                column_split = split,
                cluster_row_slices = FALSE,
                cluster_column_slices = FALSE,
                cluster_columns = FALSE,
                column_order = sample_order,
                row_names_side = "left",
                show_column_names = FALSE,
                show_row_names = TRUE,
                show_row_dend = FALSE,
                border = TRUE,
                row_names_gp = gpar(fontsize = font.size),
                row_title_gp = gpar(fontsize = 7),
                cluster_rows = FALSE,
                column_title = title,
                border_gp = gpar(lwd=0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90
  )


  ## hm9 -> MYC #########

  genes_myc <- genes_to_plot$MYC
  genes_myc <- genes_myc[which(genes_myc %in% rownames(D_norm))]

  hm9 = Heatmap(D_norm[genes_myc,,drop=FALSE],
                name="hm9_myc",
                col = col_fun,
                column_split = split,
                cluster_rows = FALSE,
                cluster_row_slices = FALSE,
                cluster_column_slices = FALSE,
                cluster_columns = FALSE,
                column_order = sample_order,
                row_names_side = "left",
                show_column_names = FALSE,
                show_row_names = TRUE,
                show_row_dend = FALSE,
                border = TRUE,
                row_names_gp = gpar(fontsize = font.size),
                row_title_gp = gpar(fontsize = 7),
                column_title = title,
                border_gp = gpar(lwd=0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90
  )


  ## hm10 -> ERBB SCORE #########

  genes_erbb <- genes_to_plot$ERBB
  genes_erbb <- genes_erbb[which(genes_erbb %in% rownames(D_norm))]

  if (length(genes_erbb) == 3) { # Check if all genes are present
    ###   ERBB SCORE
    erbb_score <- apply(D_norm, 2, function(col) sum(col[c("EGFR")]) - sum(col[c("ERBB2","ERBB3")]))

    genes_to_plot$ERBB_score <- erbb_score
    col_fun_erbb <- circlize::colorRamp2(c(quantile(erbb_score, 0.10),median(erbb_score),quantile(erbb_score, 0.90)),
                                         c("blue","white", "red"))

    col = list(ERBB_score = col_fun_erbb)

    ha10 = HeatmapAnnotation(ERBB_score = genes_to_plot$ERBB_score,
                             simple_anno_size = unit(4, "mm"),
                             simple_anno_size_adjust = TRUE,
                             annotation_name_side = "left",
                             col=col,
                             annotation_legend_param = list(title = "Scores"),
                             show_legend = FALSE,
                             border = TRUE,
                             annotation_name_gp= gpar(fontsize = 8))
    hm10 = Heatmap(D_norm[genes_erbb,,drop=FALSE],
                   name="hm10_erbb",
                   bottom_annotation = ha10,
                   col = col_fun,
                   column_split = split,
                   cluster_row_slices = FALSE,
                   cluster_rows = FALSE,
                   cluster_column_slices = FALSE,
                   cluster_columns = FALSE,
                   column_order = sample_order,
                   row_names_side = "left",
                   show_column_names = FALSE,
                   show_row_names = TRUE,
                   show_row_dend = FALSE,
                   border = TRUE,
                   row_names_gp = gpar(fontsize = font.size),
                   row_title_gp = gpar(fontsize = 7),
                   clustering_distance_rows= "pearson",
                   # clustering_distance_columns=  "pearson",
                   clustering_method_rows = "ward.D2",
                   # clustering_method_columns = "ward.D2",
                   column_title = title,
                   border_gp = gpar(lwd=0.3),
                   show_heatmap_legend = FALSE,
                   row_title_rot = 90
    )

  } else if (length(genes_erbb) != 0) { # If only some genes are missing
    message("Genes involved in the ERBB score are missing.\nERBB score will not be calculated.")

    ha10 <- NULL
    hm10 = Heatmap(D_norm[genes_erbb,,drop=FALSE],
                   name="hm10_erbb",
                   bottom_annotation = ha10,
                   col = col_fun,
                   column_split = split,
                   cluster_row_slices = FALSE,
                   cluster_rows = FALSE,
                   cluster_column_slices = FALSE,
                   cluster_columns = FALSE,
                   column_order = sample_order,
                   row_names_side = "left",
                   show_column_names = FALSE,
                   show_row_names = TRUE,
                   show_row_dend = FALSE,
                   border = TRUE,
                   row_names_gp = gpar(fontsize = font.size),
                   row_title_gp = gpar(fontsize = 7),
                   clustering_distance_rows= "pearson",
                   # clustering_distance_columns=  "pearson",
                   clustering_method_rows = "ward.D2",
                   # clustering_method_columns = "ward.D2",
                   column_title = title,
                   border_gp = gpar(lwd=0.3),
                   show_heatmap_legend = FALSE,
                   row_title_rot = 90
    )
  } else {
    message("Genes involved in the ERBB score are missing.\nERBB score will not be calculated.")
    hm10 <- NULL
  }

  ## hm11 -> ScNE #########
  # + scores for Stromal & immune infiltration
  genes_immune <- intersect(rownames(D_norm),genes_to_plot$Immune141_UP)
  genes_stroma <- intersect(rownames(D_norm),genes_to_plot$Stromal141_UP)

  # Immune
  if (length(genes_immune) != 0) { # Check if all genes are missing
    immune_score <- apply(D_norm[genes_immune,],2,median)

    genes_to_plot$Immune141_UP_score <- immune_score

    col_fun_immune <- circlize::colorRamp2(c(quantile(immune_score, 0.10),median(immune_score),quantile(immune_score, 0.90)),
                                           c("blue","white","red"))
  } else {
    message("Genes involved in the Immune141_UP score are missing.\nImmune score will not be calculated.")
    immune_score <- NULL
    col_fun_immune <- NULL
  }

  # Stroma
  if (length(genes_stroma) != 0) { # Check if all genes are missing

    stromal_score <- apply(D_norm[genes_stroma,],2,median)
    genes_to_plot$Stromal141_UP_score <- stromal_score
    col_fun_stromal <- circlize::colorRamp2(c(quantile(stromal_score, 0.10),median(stromal_score),quantile(stromal_score, 0.90)),
                                            c("blue","white","red"))
  } else {
    message("Genes involved in the Stromal141_UP score are missing.\nStromal score will not be calculated.")
    stromal_score <- NULL
    col_fun_stromal <- NULL
  }


  # ScNE genes
  col = list(Immune141_UP = col_fun_immune,
             Stromal141_UP = col_fun_stromal)

  ha11 = HeatmapAnnotation(Immune141_UP = genes_to_plot$Immune141_UP_score,
                           Stromal141_UP = genes_to_plot$Stromal141_UP_score,
                           simple_anno_size = unit(4, "mm"),
                           simple_anno_size_adjust = TRUE,
                           annotation_name_side = "left",
                           col=col,
                           show_legend = FALSE,
                           border = TRUE,
                           annotation_name_gp= gpar(fontsize = 8))
  genes_ne <- genes_to_plot$ScNE
  genes_ne <- genes_ne[which(genes_ne %in% rownames(D_norm))]

  hm11 = Heatmap(D_norm[genes_ne,,drop=FALSE],
                 name="hm11",
                 col = col_fun,
                 column_split = split,
                 bottom_annotation = ha11,
                 cluster_row_slices = FALSE,
                 cluster_column_slices = FALSE,
                 cluster_columns = FALSE,
                 column_order = sample_order,
                 row_names_side = "left",
                 show_column_names = FALSE,
                 show_row_names = TRUE,
                 show_row_dend = FALSE,
                 border = TRUE,
                 row_names_gp = gpar(fontsize = font.size),
                 row_title_gp = gpar(fontsize = 7),
                 cluster_rows = FALSE,
                 column_title = title,
                 border_gp = gpar(lwd=0.3),
                 show_heatmap_legend = show_hm_legend,
                 row_title_rot = 90
  )


  final_hm = draw(hm1 %v% hm6 %v% hm7 %v% hm4 %v% hm4_1 %v% hm2 %v% hm3 %v% hm5 %v% hm8 %v% hm9 %v% hm10 %v% hm11); hm_sample_order <- column_order(final_hm@ht_list$hm1_cc)
  invisible(hm_sample_order)

}
