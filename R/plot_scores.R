#' @title Plot Signature Scores.
#'
#' @description Build a heatmap with scores retrieved with the `predict_lundtax` function.
#'
#' @details Construct and export (pdf or png) a highly customizable heatmap visualizing prediction
#' scores for each sample and class, predicted with `predict_lundtax`.
#' This function depends on Complexheatmap. It is also possible to return a data frame with
#' prediction scores in a tidy format. To do so, set `return_scores = TRUE`. For a greater explanation
#' on how to use the function, see parameter descriptions and examples.
#'
#' @param these_predictions A list with a data frame object called scores. Returned with `predict_lundtax`.
#' @param out_path Optional, set path to export plot. If not provided, tidy version of incoming
#' scores in data frame format will be returned (`return_scores` will be auto-defaulted to TRUE).
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' if pdf, the returned pdf will be in A4 (horizontal) format, if png is specified the user can
#' control the dimensions, see `plot_width` and `plot_height`.
#' @param return_scores Set to TRUE to return prediction scores in a tidy format. Default is FALSE.
#' @param to_xlsx Boolean parameter, set to TRUE to export score data frame in xlsx format. Default is FALSE.
#' If set to TRUE, the spreadsheet will be saved to the same path as the heatmap.
#' @param title Required parameter. Heatmap title, will also be pasted to the exported file(s) as well
#' as a new column in the scores data frame under cohort.
#' @param hm_split Optional parameter for controlling how the data is split into different groups.
#' If not provided, the function will split on subtype.
#' @param hm_cluster Boolean parameter, set to TRUE to cluster the rows (default is FALSE).
#' @param plot_anno_legend Expects a vector with TRUE/FALSE (7 in total), thsi decides what legends
#' will be on the final heatmap. Default is to only show the legend for the subtypes.
#' @param plot_hm_legend Boolean parameter. Set to TRUE to show heatmap legend. Default is FALSE.
#' @param plot_width This parameter controls the width in inches. Default is 14(4200 pixels at 300 PPI).
#' @param plot_height This parameter controls the height in inches.Default is 6(1800 pixels at 300 PPI)
#' @param plot_font_size Optional parameter to control the size of the font in the generated heatmap.
#' Note, the title of the plot will always be twice that of the set value here (default = 12).
#' @param verbose Set to TRUE for debugging purposes. Default is FALSE.
#'
#' @return Data frame with prediction score for each sample and class, if return_scores = TRUE.
#' Otherwise, nothing.
#'
#' @import ComplexHeatmap ggplot2 circlize
#' @importFrom stats median quantile
#' @importFrom grDevices dev.off pdf png
#' @importFrom openxlsx write.xlsx
#' 
#'
#' @export
#'
#' @examples
#' \dontrun{
#' my_predictions = lundtax_predict_sub(these_predictions = sjodahl_2017, 
#'                                      gene_id = "hgnc_symbol", 
#'                                      impute = TRUE, 
#'                                      adjust = TRUE)
#' 
#' plot_scores(these_predictions = my_predictions,
#'             out_path = "../",
#'             out_format = "pdf",
#'             title = "Lund2017")
#'}
#'
plot_scores = function(these_predictions = NULL,
                       out_path = NULL,
                       out_format = "png",
                       return_scores = FALSE,
                       to_xlsx = FALSE,
                       title = NULL,
                       hm_split = NULL,
                       hm_cluster = FALSE,
                       plot_anno_legend = NULL,
                       plot_hm_legend = FALSE,
                       plot_width = 14,
                       plot_height = 6,
                       plot_font_size = 12,
                       verbose = TRUE){

  #deal with nonsensical parameter combinations
  if(to_xlsx){
    return_scores = TRUE
  }

  if(is.null(title)){
    stop("Please specify a title for your plot with `title`...")
  }

  #set the default legend properties
  if(is.null(plot_anno_legend)){
    plot_anno_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  }else{
    plot_anno_legend = plot_anno_legend
    }

  #check the type of incoming data
  if(is.null(these_predictions) | !is.list(these_predictions)){
    stop("Input should be the result of applying predict_lundtax")
  }

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

    #set sample order
    sample_order <- order(these_predictions$scores$proliferation_score)

    #set split (top annotation track of heatmap)
    if(!is.null(hm_split)){
      split = hm_split
      }else{
        split <- factor(these_predictions$predictions_7classes,
                        levels = c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))
     }

    #get immune names
    immune_names <- c("immune141_up", "b_cells", "t_cells",
                      "t_cells_cd8", "nk_cells",
                      "cytotoxicity_score", "neutrophils",
                      "monocytic_lineage", "macrophages",
                      "m2_macrophage", "myeloid_dendritic_cells")

    #get stromal names
    stromal_names <- c("stromal141_up", "endothelial_cells",
                       "fibroblasts", "smooth_muscle")

    #set colours
    #proliferation
    col_fun_proliferation =
      circlize::colorRamp2(c(quantile(these_predictions$scores$proliferation_score, 0.05),
                             median(these_predictions$scores$proliferation_score),
                             quantile(these_predictions$scores$proliferation_score, 0.95)),
                           c("#21908CFF","white", "#B63679FF"))

    #progression
    col_fun_progression =
      circlize::colorRamp2(c(quantile(these_predictions$scores$progression_score, 0.05),
                             median(these_predictions$scores$progression_score),
                             quantile(these_predictions$scores$progression_score, 0.90)),
                           c("#FAEBDDFF","#A11A5BFF", "#4C1D4BFF"))

    #create colour object
    colour_obj = list(lund_subtype = lund_colors$lund_colors,
                      proliferation_score = col_fun_proliferation,
                      molecular_grade_who_1999 = c("G1_2" = "white", "G3" = "black"),
                      molecular_grade_who_2016 = c("HG" = "black", "LG" = "white"),
                      progression_score = col_fun_progression,
                      progression_risk = c("HR" = "#A11A5BFF", "LR" = "#FAEBDDFF"))

    #plotting
    #build heatmap annotation (top)
    hm <- HeatmapAnnotation(lund_subtype = split,
                            proliferation_score = these_predictions$scores$proliferation_score,lund_subtype = split,
                            molecular_grade_who_1999 = these_predictions$scores$molecular_grade_who_1999,
                            molecular_grade_who_2016 = these_predictions$scores$molecular_grade_who_2016,
                            progression_score = these_predictions$scores$progression_score,
                            progression_risk = these_predictions$scores$progression_risk,
                            annotation_name_side = "left",
                            show_legend = plot_anno_legend,
                            annotation_name_gp = gpar(fontsize = plot_font_size),
                            border = TRUE,
                            col = colour_obj,
                            column_title = NULL,
                            row_title = NULL)

    #immune scores heatmap
    hm_immune_scores <- Heatmap(t(scale(these_predictions$scores[,immune_names, drop = FALSE])),
                                top_annotation = hm,
                                column_order = sample_order,
                                height = unit(5*ncol(these_predictions$scores[,immune_names,drop = FALSE]), "mm"),
                                width = unit(0.5*nrow(these_predictions$scores[,immune_names,drop = FALSE]), "mm"),
                                name = "Immune Scores",
                                border = TRUE,
                                column_split = split,
                                cluster_rows = hm_cluster,
                                show_heatmap_legend = plot_hm_legend,
                                row_names_side = "left",
                                row_names_gp = gpar(fontsize = plot_font_size),
                                show_column_names = FALSE,
                                column_title = NULL,
                                row_title = NULL)

    #stroma scores heatmap
    hm_stroma_scores <- Heatmap(t(scale(these_predictions$scores[,stromal_names, drop = FALSE])),
                                column_order = sample_order,
                                height = unit(5*ncol(these_predictions$scores[,stromal_names, drop = FALSE]), "mm"),
                                width = unit(0.5*nrow(these_predictions$scores[,stromal_names,drop = FALSE]), "mm"),
                                name = "Stroma Scores",
                                border = TRUE,
                                column_split = split,
                                cluster_rows = hm_cluster,
                                show_heatmap_legend = plot_hm_legend,
                                row_names_side = "left",
                                row_names_gp = gpar(fontsize = plot_font_size),
                                show_column_names = FALSE,
                                column_title = NULL,
                                row_title = NULL)

    #combine heatmaps
    hm_combined = draw(hm_immune_scores %v% hm_stroma_scores,
                       column_title = title,
                       column_title_gp = gpar("fontface", fontsize = (plot_font_size*2)));

    hm_sample_order <- column_order(hm_combined@ht_list$`Stroma Scores`)

    dev.off()

    if(verbose){
      message(paste0("Heatmap exported to ", out_path, title, "_heatmap_scores.", out_format))
      }

  }else{
    message("No output path provided, returning the score data frame in tidy format...")
    return_scores = TRUE
    }

  if(return_scores == TRUE){

    #subset scores data frame from list
    my_scores = these_predictions$scores
    my_scores = as.data.frame(my_scores)

    #add new column with cohort information
    my_scores$cohort = title

    #convert correct columns to factors
    factor_cols <- c('molecular_grade_who_1999' ,'molecular_grade_who_2016',
                     'progression_risk', 'cohort')

    my_scores[,factor_cols] <- lapply(my_scores[,factor_cols] , factor)

    if(to_xlsx){
      my_scores <- tibble::rownames_to_column(my_scores, "SampleID")
      write.xlsx(my_scores, paste0(out_path, title, "_scores.xlsx"))
      if(verbose){
        message(paste0("Prediction scores exported as xlsx to ", out_path, title, "_scores.xlsx"))
        }
    }
    
    dev.off()

    #return data frame
    return(my_scores)

    if(verbose){
      message("Prediction Scores for each sample and class sucessfully generated!")
      }
    }
  }
