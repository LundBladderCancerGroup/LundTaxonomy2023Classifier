#' @title Plot Heatmap Signatures.
#'
#' @description Plot heatmap for classification results.
#'
#' @details This function plots a heatmap including genes and signatures of 
#' interest, with prediction results and scores on top.
#' 
#' @param these_predictions Required parameter, should be the output from 
#' [LundTax2023Classifier::lundtax_predict_sub()].
#' @param this_data Expression data used for predictions. 
#' Required if the output from [LundTax2023Classifier::lundtax_predict_sub()] is run with 
#' include_data = FALSE (default).
#' @param gene_id Specify the type of gene identifier used in `this_data`. Accepted values are; 
#' hgnc_symbol (default) or ensembl_gene_id.
#' @param subtype_annotation Can be one of the following; "5 _class" (default) or "7_class" 
#' annotation.
#' @param norm Boolean parameter. Set to TRUE (default) to normalize the data into Z-scaled values.
#' @param plot_scores Boolean parameter. Set to TRUE (default) to plot prediction scores for each 
#' class.
#' @param show_ann_legend Boolean parameter, set to TRUE to show annotation legend (Lund classes). 
#' Default is FALSE.
#' @param show_hm_legend Boolean parameter, set to TRUE to show heatmap legend, default is FALSE.
#' @param ann_height Plotting parameter, optional. Annotation height in cm. Default = 8.
#' @param title Plotting parameter. The title for the generated heatmap. Deafult is "My Plot".
#' @param plot_width This parameter controls the width in inches. Default is 14 (4200 pixels at 300 
#' PPI).
#' @param plot_height This parameter controls the height in pixels. Default is 10 (3000 pixels at 
#' 300 PPI)
#' @param plot_font_size Optional parameter to control the size of the font in the generated 
#' heatmap. Note, the title of the plot will always be twice that of the set value here 
#' (default = 10).
#' @param plot_font_row_size Optional parameter to control the size of the font in the generated 
#' heatmap. Note, the title of the plot will always be twice that of the set value here 
#' (default = 8).
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can control the dimensions with `plot_width` and `plot_height`.
#' 
#' @return Draws heatmap and silently returns the sample order.
#' 
#' @import ComplexHeatmap circlize dplyr
#' @importFrom stats na.omit quantile
#' @importFrom grDevices dev.off pdf png
#' 
#' @export 
#'
#' @examples
#' \dontrun{ 
#' #example 1 including data in results object
#' #run predictor on the bundled expression data
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017,
#'                                         include_data = TRUE, 
#'                                         impute = TRUE)
#'

#'
#'
#' #example 2 - 5 class annotation and Without prediction scores
#' plot_hm_signatures(these_predictions = sjodahl_predicted, 
#'                    subtype_annotation = "5_class",
#'                    ann_height = 0.5,
#'                    plot_scores = FALSE)
#'}
#'
plot_hm_signatures = function(these_predictions = NULL,
                              this_data = NULL,
                              gene_id = "hgnc_symbols",
                              subtype_annotation = "5_class",
                              norm = TRUE,
                              plot_scores = TRUE,
                              show_ann_legend = FALSE,
                              show_hm_legend = FALSE,
                              ann_height = 8,
                              title = "My Plot",
                              plot_width = 14,
                              plot_height = 10,
                              plot_font_size = 10,
                              plot_font_row_size = 8,
                              out_path = NULL,
                              out_format = "png"){
  
  #check incoming data and parameter combinations
  if(is.null(these_predictions) | !is.list(these_predictions)){
    stop("Input should be the result of applying predict_lundtax...")
  }else if(is.list(these_predictions) & "data" %in% names(these_predictions)){
    this_data = these_predictions$data
    score_matrix = these_predictions$subtype_scores
    pred_labels5 = these_predictions$predictions_5classes
    pred_labels7 = these_predictions$predictions_7classes
  }else if(is.list(these_predictions) & !"data" %in% names(these_predictions)){
    if(is.null(this_data)){
      stop("Data is missing. Include it in results object by running predict_lundtax with include_data = TRUE or provide it in the this_data argument...")
    }else if(!(class(this_data)[1] %in%  c("matrix","data.frame"))){
      stop("Data should be in matrix or data.frame format")
    }else if(class(this_data)[1] %in%  c("matrix","data.frame")){
      this_data = this_data[,names(these_predictions$predictions_7classes)]
      score_matrix = these_predictions$subtype_scores
      pred_labels5 = these_predictions$predictions_5classes
      pred_labels7 = these_predictions$predictions_7classes
    }
  }
  
  #convert ensembl gene id to hgnc symbols
  if(gene_id == "ensembl_gene_id"){
    
    #convert the rownames to first column
    mutated_data = dplyr::as_tibble(this_data, 
                                    rownames = "ensembl_gene_id")
    
    #left join with gene list to get ensembl IDs
    mutated_data = dplyr::left_join(gene_list, 
                                    mutated_data, 
                                    by = "ensembl_gene_id")
    
    #remove duplicated rows
    mutated_data = dplyr::filter(mutated_data, 
                                 duplicated(hgnc_symbol) == FALSE)
    
    #convert the first column back to rownames
    row.names(mutated_data) = unique(mutated_data$hgnc_symbol)
    mutated_data[1:2] = NULL
    
    #convert back to expected name
    this_data = mutated_data
  }

  #scale and set this_data as matrix
  if(norm){
    this_data = as.matrix(this_data)
    this_data = scale(t(this_data))
    this_data = t(this_data)
  }else{
    this_data = as.matrix(this_data)
  }
  
  #gene signatures for plotting
  genes_to_plot = list(Early_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "EarlyCellCycle"), 1]),
                       Late_CC = c(signatures$proliferation[which(signatures$proliferation$signature == "LateCellCycle"), 1]),
                       Late_Early = NULL,
                       UroDiff = c("PPARG", "FOXA1", "GATA3", "ELF3"),
                       UPKs = c("UPK1A", "UPK1B", "UPK2", "UPK3A", "KRT20"),
                       Circuit = c("FGFR3", "CCND1", "E2F3", "RB1", "CDKN2A"),
                       Circuit_score = NULL,
                       FGFR3 = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "FGFR3"), 1]),
                       BaSq = c("KRT5", "KRT14", "FOXA1", "GATA3"),
                       BaSq_ratio = NULL,
                       Keratinization = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "Keratinization_QTC"), 1]),
                       Adhesion = c("EPCAM", "CDH1", "CDH3"),
                       MYC = c("MYCL", "MYCN", "MYC"),
                       ERBB = c("EGFR", "ERBB2", "ERBB3"),
                       ERBB_score = NULL,
                       ScNE = c("CHGA", "SYP", "ENO2"),
                       Immune141_UP = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "Immune141_UP"), 1]),
                       Stromal141_UP = c(signatures$signatures_plot[which(signatures$signatures_plot$signature == "Stromal141_UP"), 1]),
                       Immune141_UP_score = NULL,
                       Stromal141_UP_score = NULL)
  
  #create wrapper function for subtype score plots
  bar_anno = function(plot_scores = plot_scores,
                      subtype = NULL){
    
    if(plot_scores){
      this_bar = anno_barplot(as.numeric(score_matrix[,subtype]), ylim = c(0, 1), gp = gpar(fill = lund_colors$lund_colors[subtype], border = NA, col = NA), bar_width = 1, height = unit(6, "mm"))
    }else{
      this_bar = NULL
    }
    return(this_bar)
  }

  #plotting
  ##heatmap annotations - subtype predictions with scores
  #7 classes
  if(subtype_annotation == "7_class"){
    
    #predictions
    pred_lab = pred_labels7
    
    #column split
    split = factor(pred_lab ,levels = c("UroA" ,"UroB", "UroC", "GU", "BaSq", "Mes", "ScNE"))
  
    #score plots
    bar0 = bar_anno(plot_scores = plot_scores, subtype = "Uro")
    bar1 = bar_anno(plot_scores = plot_scores, subtype = "UroA")
    bar2 = bar_anno(plot_scores = plot_scores, subtype = "UroB")
    bar3 = bar_anno(plot_scores = plot_scores, subtype = "UroC")
    bar4 = bar_anno(plot_scores = plot_scores, subtype = "GU")
    bar5 = bar_anno(plot_scores = plot_scores, subtype = "BaSq")
    bar6 = bar_anno(plot_scores = plot_scores, subtype = "Mes")
    bar7 = bar_anno(plot_scores = plot_scores, subtype = "ScNE")
    
    #colors
    col = list(Predictions = lund_colors$lund_colors)
    
    #heatmap annotations
    ha1 = HeatmapAnnotation(Predictions = pred_lab,
                            annotation_name_side = "left",
                            col = col,
                            Uro = bar0,
                            UroA = bar1,
                            UroB = bar2,
                            UroC = bar3,
                            GU = bar4, 
                            BaSq = bar5, 
                            Mes = bar6,
                            ScNE = bar7,
                            na_col = "gray83", 
                            gap = unit(2, "mm"),
                            simple_anno_size = unit(10, "mm"),
                            simple_anno_size_adjust = TRUE,
                            show_legend = show_ann_legend,
                            border = TRUE,
                            height = unit(ann_height, "cm"),
                            annotation_name_gp = gpar(fontsize = plot_font_size))

  #5 classes
  }else if(subtype_annotation == "5_class"){
  
    #predictions
    pred_lab = pred_labels5
  
    #column split
    split = factor(pred_lab, levels = c("Uro", "GU", "BaSq", "Mes", "ScNE"))
  
    #score plots
    bar1 = bar_anno(plot_scores = plot_scores, subtype = "Uro")
    bar2 = bar_anno(plot_scores = plot_scores, subtype = "GU")
    bar3 = bar_anno(plot_scores = plot_scores, subtype = "BaSq")
    bar4 = bar_anno(plot_scores = plot_scores, subtype = "Mes")
    bar5 = bar_anno(plot_scores = plot_scores, subtype = "ScNE")

    #colors
    col = list(Predictions = lund_colors$lund_colors)

    #draw heatmap for subtype predictions
    ha1 = HeatmapAnnotation(Predictions = pred_lab,
                            annotation_name_side = "left",
                            col = col,
                            Uro = bar1, 
                            GU = bar2,
                            BaSq = bar3, 
                            Mes = bar4,
                            ScNE = bar5,
                            na_col = "gray83",
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            show_legend = show_ann_legend,
                            border = TRUE,
                            height = unit(ann_height, "cm"),
                            annotation_name_gp = gpar(fontsize = plot_font_size))
  }else{
    stop("Please provide a valid subtype annotation (5 classes or 7 classes)...")
  }
  
  ##heatmap 1 - late/early cell cycle
  #get genes in provided data for downstream filtering steps
  these_genes = row.names(this_data)
  genes_early = intersect(genes_to_plot$Early_CC, these_genes)
  genes_late = intersect(genes_to_plot$Late_CC, these_genes)

  #combine genes
  genes_cc = na.omit(c(genes_early, genes_late))

  #check if genes from both early and late are present
  if(length(genes_early) != 0 & !all(is.na(genes_late))){
    
    #row split for the heatmap
    row_split = c(rep("Early", length(genes_early)),
                  rep("Late", length(genes_late)))
    #row title
    row_title_cc = c("Late Cell Cycle", "Early Cell Cycle")

    #late and Early scores
    late_score = apply(this_data[intersect(rownames(this_data),genes_to_plot$Late_CC),], 2, median)
    early_score = apply(this_data[intersect(rownames(this_data),genes_to_plot$Early_CC),], 2, median)

    #ratio
    late_early = late_score - early_score

    #add genes to genes_to_plot object
    genes_to_plot$Late_Early = late_early
    
    #create color palette for late/early
    col_fun_cc = circlize::colorRamp2(c(quantile(late_early, 0.05),
                                        median(late_early),
                                        quantile(late_early, 0.95)),
                                      c("blue","white", "red"))

    #order samples by late_early cell cycle
    sample_order = order(late_early)

    #heatmap annotation
    col = list(Predictions = lund_colors$lund_colors,
               late_early = col_fun_cc)
    
    #draw annotations track, late/early
    ha1b = HeatmapAnnotation(late_early = genes_to_plot$Late_Early,
                             annotation_name_side = "left",
                             simple_anno_size = unit(4, "mm"),
                             simple_anno_size_adjust = TRUE,
                             col = col,
                             show_legend = show_ann_legend,
                             border = TRUE,
                             annotation_name_gp = gpar(fontsize = plot_font_row_size))

  }else if(length(genes_early) == 0 & !all(is.na(genes_late))){ #if early cell cycle is missing
    message("All genes from the early cell cycle signature are missing.
            \nSamples will not be ordered by cell cycle...")
    sample_order = NULL
    ha1b = NULL
    row_split = NULL
    row_title_cc = "Late cell cycle"
  }else if(length(genes_early) != 0 & all(is.na(genes_late))){ #if late cell cycle is missing
    message("All genes from the late cell cycle signature are missing.
            \nSamples will not be ordered by cell cycle...")
    sample_order = NULL
    ha1b = NULL
    row_split = NULL
    row_title_cc = "Early cell cycle"
  }else{ #both are missing
    message("All genes from the late/early cell cycle signatures are missing.
            \nSamples will not be ordered by cell cycle...")
    sample_order = NULL
    ha1b = NULL
    row_split = NULL
    row_title_cc = NULL
  }

  #heatmap colors
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

  #draw heatmap- late/early
  hm1 = Heatmap(this_data[genes_cc,, drop = FALSE],
                top_annotation = ha1,
                bottom_annotation = ha1b,
                name = "hm1_cc",
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
                row_names_gp = gpar(fontsize = plot_font_row_size),
                row_title_gp = gpar(fontsize = plot_font_size),
                clustering_distance_rows = "spearman",
                clustering_method_rows = "ward.D2",
                column_title = title,
                height = 8,
                border_gp = gpar(lwd = 0.3),
                show_heatmap_legend = FALSE)

  ##heatmap 2 - urodiff
  genes_ud = genes_to_plot$UroDiff
  genes_ud = genes_ud[which(genes_ud %in% rownames(this_data))]

  hm2 = Heatmap(this_data[genes_ud,,drop = FALSE],
                name = "hm2_ud",
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
                row_names_gp = gpar(fontsize = plot_font_row_size),
                row_title_gp = gpar(fontsize = plot_font_size),
                cluster_rows = FALSE,
                column_title = title,
                border_gp = gpar(lwd=0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90)

  ##heatmap 3 - uroplakins
  genes_upk = genes_to_plot$UPKs
  genes_upk = genes_upk[which(genes_upk %in% rownames(this_data))]

  hm3 = Heatmap(this_data[genes_upk,,drop = FALSE],
                name = "hm3_upk",
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
                row_names_gp = gpar(fontsize = plot_font_row_size),
                row_title_gp = gpar(fontsize = plot_font_size),
                cluster_rows = FALSE,
                column_title = title,
                border_gp = gpar(lwd = 0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90)

  ##heatmap 4 - circuit score
  genes_circ = genes_to_plot$Circuit
  genes_circ = genes_circ[which(genes_circ %in% rownames(this_data))]

  #check if all genes are present
  if(length(genes_circ) == 5){

    #calculate circuit sore
    circuit_score = apply(this_data, 2, function(col) sum(col[c("RB1", "FGFR3", "CCND1")]) - sum(col[c("E2F3", "CDKN2A")]))
  
    #add scores to genes_to_plot object
    genes_to_plot$Circuit_score = circuit_score

    #generate color palette
    col_fun_circ = circlize::colorRamp2(c(quantile(circuit_score, 0.10),
                                          median(circuit_score),
                                          quantile(circuit_score, 0.90)),
                                         c("blue","white", "red"))

    col = list(circuit_score = col_fun_circ)

    #create annotation track for heatmap 4
    ha4 = HeatmapAnnotation(circuit_score = genes_to_plot$Circuit_score,
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            annotation_name_side = "left",
                            col = col,
                            show_legend = FALSE,
                            border = TRUE,
                            annotation_name_gp = gpar(fontsize = plot_font_row_size))

    #draw heatmap 4 - circuit score
    hm4 = Heatmap(this_data[genes_circ,, drop = FALSE],
                  name = "hm4_circ",
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
                  row_names_gp = gpar(fontsize = plot_font_row_size),
                  row_title_gp = gpar(fontsize = plot_font_size),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd = 0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90)
  
  }else if(length(genes_circ) != 0){
    ha4 = NULL
    message("Genes involved in the circuit score are missing.
            \nCircuit score will not be calculated...")

    #draw heatmap 4 - circuit score
    hm4 = Heatmap(this_data[genes_circ,,drop = FALSE],
                  name = "hm4_circ",
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
                  row_names_gp = gpar(fontsize = plot_font_row_size),
                  row_title_gp = gpar(fontsize = plot_font_size),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd = 0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90)
  }else{
    message("No genes involved in the circuit score are found.
            \nCircuit score will not be calculated...")
    hm4 = NULL
  }

  ##heatmap 4_1 - TP63
  if("TP63" %in% rownames(this_data)){
    tp63 = this_data["TP63",, drop = FALSE]
    rownames(tp63) = c("TP63")

    #draw heatm 4_1 - TP63
    hm4_1 = Heatmap(tp63,
                    name = "tp63",
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
                    row_names_gp = gpar(fontsize = plot_font_row_size),
                    row_title_gp = gpar(fontsize = plot_font_size),
                    column_title = title,
                    border_gp = gpar(lwd = 0.3),
                    show_heatmap_legend = FALSE,
                    row_title_rot = 90)
  }else{
    message("TP63 was not found in the expresison matrix.
            \nNo heatmap will be generated for this gene...")
    hm4_1 = NULL
  }

  ##heatmap 5 - FGFR3
  genes_fgfr3 = genes_to_plot$FGFR3
  genes_fgfr3 = genes_fgfr3[which(genes_fgfr3 %in% rownames(this_data))]

  if(length(genes_fgfr3) != 0){
    hm5 = Heatmap(this_data[genes_fgfr3,, drop = FALSE],
                  name = "hm5_fgfr3",
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
                  row_names_gp = gpar(fontsize = plot_font_row_size),
                  row_title_gp = gpar(fontsize = plot_font_size),
                  clustering_distance_rows= "pearson",
                  clustering_method_rows = "ward.D2",
                  column_title = title,
                  height = 4,
                  border_gp = gpar(lwd = 0.3),
                  show_heatmap_legend = FALSE)
  }else{
    message("FGFR3 signature genes was not found in the expresison matrix.
             \nNo heatmap will be generated for this gene...")
    hm5 = NULL
  }

  ##heatmap 6 - Ba/Sq ratio
  genes_basq = genes_to_plot$BaSq
  genes_basq = genes_basq[which(genes_basq %in% rownames(this_data))]

  #check if all genes are present
  if(length(genes_basq) == 4) {
    
    #calculate BaSq ratio
    basq_ratio = apply(this_data, 2, function(col) sum(col[c("KRT5", "KRT14", "FOXA1")]) - sum(col[c("GATA3")]))

    genes_to_plot$BaSq_ratio = basq_ratio
    
    #generate color palette
    col_fun_basq = circlize::colorRamp2(c(quantile(basq_ratio, 0.10),
                                          median(basq_ratio),
                                          quantile(basq_ratio, 0.90)),
                                        c("blue", "white", "red"))

    col = list(BaSq_ratio = col_fun_basq)

    #draw annotation track for heatmap 6 - BaSq ratio
    ha6 = HeatmapAnnotation(BaSq_ratio = genes_to_plot$BaSq_ratio,
                            simple_anno_size = unit(4, "mm"),
                            simple_anno_size_adjust = TRUE,
                            annotation_name_side = "left",
                            col = col,
                            show_legend = FALSE,
                            border = TRUE,
                            annotation_name_gp= gpar(fontsize = plot_font_row_size))
    
    #draw heatmap 6 - BaSq ratio
    hm6 = Heatmap(this_data[genes_basq,, drop = FALSE],
                  name = "hm6_basq",
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
                  row_names_gp = gpar(fontsize = plot_font_row_size),
                  row_title_gp = gpar(fontsize = plot_font_size),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd = 0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90)

  }else if(length(genes_basq) != 0){
    #draw annotation track for heatmap 6 - BaSq ratio
    ha6 = NULL
    message("Genes involved in the BaSq score are missing.
            \nBaSq ratio will not be calculated...")
    
    #draw heatmap 6 - BaSq ratio
    hm6 = Heatmap(this_data[genes_basq,,drop = FALSE],
                  name = "hm6_basq",
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
                  row_names_gp = gpar(fontsize = plot_font_row_size),
                  row_title_gp = gpar(fontsize = plot_font_size),
                  cluster_rows = FALSE,
                  column_title = title,
                  border_gp = gpar(lwd = 0.3),
                  show_heatmap_legend = FALSE,
                  row_title_rot = 90)
  }else{
    message("Genes involved in the BaSq score are missing.
            \nNo heatmap forBaSq ratio will be generated...")
    hm6 = NULL
    col_fun_basq = NULL
  }

  ##heatmap 7 - keratinization
  genes_krt = genes_to_plot$Keratinization
  genes_krt = genes_krt[which(genes_krt %in% rownames(this_data))]

  if(length(genes_krt) != 0){
    hm7 = Heatmap(this_data[genes_krt,,drop = FALSE],
                  name = "hm7_krt",
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
                  row_names_gp = gpar(fontsize = plot_font_row_size),
                  row_title_gp = gpar(fontsize = plot_font_size),
                  clustering_distance_rows= "pearson",
                  clustering_method_rows = "ward.D2",
                  column_title = title,
                  height = 4,
                  border_gp = gpar(lwd = 0.3),
                  show_heatmap_legend = FALSE)
  }else{
    message("Genes involved in the Keratinization signatures are missing.
            \nNo heatmap for keratinization signature ratio will be generated...")
    hm7 = NULL
  }

  ##heatmap 8 - adhesion
  genes_ad = genes_to_plot$Adhesion
  genes_ad = genes_ad[which(genes_ad %in% rownames(this_data))]

  hm8 = Heatmap(this_data[genes_ad,,drop = FALSE],
                name = "hm8_ad",
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
                row_names_gp = gpar(fontsize = plot_font_row_size),
                row_title_gp = gpar(fontsize = plot_font_size),
                cluster_rows = FALSE,
                column_title = title,
                border_gp = gpar(lwd = 0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90)

  ##heatmap 9 - MYC
  genes_myc = genes_to_plot$MYC
  genes_myc = genes_myc[which(genes_myc %in% rownames(this_data))]

  hm9 = Heatmap(this_data[genes_myc,,drop = FALSE],
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
                row_names_gp = gpar(fontsize = plot_font_row_size),
                row_title_gp = gpar(fontsize = plot_font_size),
                column_title = title,
                border_gp = gpar(lwd = 0.3),
                show_heatmap_legend = FALSE,
                row_title_rot = 90)


  ##heatmap 10 - ERBB score
  genes_erbb = genes_to_plot$ERBB
  genes_erbb = genes_erbb[which(genes_erbb %in% rownames(this_data))]

  #check if all genes are present
  if(length(genes_erbb) == 3){
    
    #calculate erbb score
    erbb_score = apply(this_data, 2, function(col) sum(col[c("EGFR")]) - sum(col[c("ERBB2", "ERBB3")]))

    #add scores back to the genes_to_plot object
    genes_to_plot$ERBB_score = erbb_score
    
    #generate color palette
    col_fun_erbb = circlize::colorRamp2(c(quantile(erbb_score, 0.10),
                                          median(erbb_score),
                                          quantile(erbb_score, 0.90)),
                                         c("blue", "white", "red"))

    col = list(ERBB_score = col_fun_erbb)

    #draw annotation track for ERBB score heatmap (heatmap 10)
    ha10 = HeatmapAnnotation(ERBB_score = genes_to_plot$ERBB_score,
                             simple_anno_size = unit(4, "mm"),
                             simple_anno_size_adjust = TRUE,
                             annotation_name_side = "left",
                             col=col,
                             annotation_legend_param = list(title = "Scores"),
                             show_legend = FALSE,
                             border = TRUE,
                             annotation_name_gp = gpar(fontsize = plot_font_row_size))
    
    #draw heatmap 10, ERBB scores
    hm10 = Heatmap(this_data[genes_erbb,,drop = FALSE],
                   name = "hm10_erbb",
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
                   row_names_gp = gpar(fontsize = plot_font_row_size),
                   row_title_gp = gpar(fontsize = plot_font_size),
                   clustering_distance_rows= "pearson",
                   clustering_method_rows = "ward.D2",
                   column_title = title,
                   border_gp = gpar(lwd = 0.3),
                   show_heatmap_legend = FALSE,
                   row_title_rot = 90)

  }else if(length(genes_erbb) != 0){ 
    message("Some genes involved in the ERBB score are missing.
            \nERBB score will not be calculated.")
    
    #draw annotation track for ERBB score heatmap (heatmap 10)
    ha10 = NULL
    
    #draw heatmap 10, ERBB scores
    hm10 = Heatmap(this_data[genes_erbb,,drop = FALSE],
                   name = "hm10_erbb",
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
                   row_names_gp = gpar(fontsize = plot_font_row_size),
                   row_title_gp = gpar(fontsize = plot_font_size),
                   clustering_distance_rows= "pearson",
                   clustering_method_rows = "ward.D2",
                   column_title = title,
                   border_gp = gpar(lwd = 0.3),
                   show_heatmap_legend = FALSE,
                   row_title_rot = 90)
  }else{
    message("Genes involved in the ERBB scores are missing.
            \nNo heatmap associated with ERBB will be generated...")
    hm10 = NULL
  }

  ##heatmap 11 - ScNE, scores for Stromal & immune infiltration
  genes_immune = intersect(rownames(this_data),genes_to_plot$Immune141_UP)
  genes_stroma = intersect(rownames(this_data),genes_to_plot$Stromal141_UP)

  #immune
  if(length(genes_immune) != 0){
    
    #calculate immune score
    immune_score = apply(this_data[genes_immune,],2,median)

    #add immune scores back to genes_to_plot object
    genes_to_plot$Immune141_UP_score = immune_score

    #generate color palette
    col_fun_immune = circlize::colorRamp2(c(quantile(immune_score, 0.10),
                                            median(immune_score),
                                            quantile(immune_score, 0.90)),
                                          c("blue", "white", "red"))
  }else{
    message("Genes involved in the Immune141_UP score are missing.
            \nImmune score will not be calculated.")
    immune_score = NULL
    col_fun_immune = NULL
  }

  #stroma
  if(length(genes_stroma) != 0){

    #calculate stromal score
    stromal_score <- apply(this_data[genes_stroma,],2,median)
    
    #add stromal scores back to genes_to_plot object
    genes_to_plot$Stromal141_UP_score = stromal_score
    
    #generate color palette
    col_fun_stromal = circlize::colorRamp2(c(quantile(stromal_score, 0.10),
                                             median(stromal_score),
                                             quantile(stromal_score, 0.90)),
                                           c("blue", "white", "red"))
  }else{
    message("Genes involved in the Stromal141_UP score are missing.
            \nStromal score will not be calculated.")
    stromal_score = NULL
    col_fun_stromal = NULL
  }

  col = list(Immune141_UP = col_fun_immune,
             Stromal141_UP = col_fun_stromal)

  #draw annotation track for heatmap 11 - SnNE and immune/stromal scores
  ha11 = HeatmapAnnotation(Immune141_UP = genes_to_plot$Immune141_UP_score,
                           Stromal141_UP = genes_to_plot$Stromal141_UP_score,
                           simple_anno_size = unit(4, "mm"),
                           simple_anno_size_adjust = TRUE,
                           annotation_name_side = "left",
                           col = col,
                           show_legend = FALSE,
                           border = TRUE,
                           annotation_name_gp = gpar(fontsize = plot_font_row_size))
  
  #ScNE genes
  genes_ne = genes_to_plot$ScNE
  genes_ne = genes_ne[which(genes_ne %in% rownames(this_data))]

  #draw heatmap 11 - ScNE, immune scores and stromal scores
  hm11 = Heatmap(this_data[genes_ne,,drop = FALSE],
                 name = "hm11",
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
                 row_names_gp = gpar(fontsize = plot_font_row_size),
                 row_title_gp = gpar(fontsize = plot_font_size),
                 cluster_rows = FALSE,
                 column_title = title,
                 border_gp = gpar(lwd = 0.3),
                 show_heatmap_legend = show_hm_legend,
                 row_title_rot = 90)
  
  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, title, "_heatmap_scores.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, title, "heatmap_scores.png"),
          width = plot_width,
          height = plot_height,
          units = "in",
          res = 300,
          pointsize = 10,
          bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
    }
  }else{
    message("No out_path provided, the function will return the heatmap within your R session...")
  }

  #combine all heatmaps to the final object
  final_hm = draw(hm1 %v% hm6 %v% hm7 %v% hm4 %v% hm4_1 %v% hm2 %v% hm3 %v% hm5 %v% hm8 %v% hm9 %v% hm10 %v% hm11); hm_sample_order <- column_order(final_hm@ht_list$hm1_cc)
  invisible(hm_sample_order)
  dev.off()
}
