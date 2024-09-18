#' @title Plot Signature Scores.
#'
#' @description Build a heatmap with scores retrieved with the predict_LundTax2023 function.
#'
#' @details Construct and export (pdf or png) a highly customizable heatmap visualizing prediction
#' scores for each sample and class, predicted with `predict_lundtax`.
#' This function depends on Complexheatmap. It is also possible to return a data frame with
#' prediction scores in a tidy format. To do so, set `return_scores = TRUE`. For a greater explanation
#' on how to use the function, see parameter descriptions and examples.
#'
#' @param this_data A list with a data frame object called scores. Returned with `predict_lundtax`.
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
#' @param plot_width If `out_format = "png"`, this parameter controls the width in pixels.
#' Default is 14 inches (4200 pixels at 300 PPI).
#' @param plot_height If `out_format = "png"`, this parameter controls the height in pixels.
#' Default is 6 inches (1800 pixels at 300 PPI)
#' @param plot_font_size Optional parameter to control the size of the font in the generated heatmap.
#' Note, the title of the plot will always be twice that of the set value here (default = 12).
#' @param verbose Set to TRUE for debugging purposes. Default is FALSE.
#'
#' @return Data frame with prediction score for each sample and class, if return_scores = TRUE.
#' Otherwise, nothing.
#'
#' @import ComplexHeatmap ggplot2 circlize openxlsx grDevices utils grid
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #load pacakges
#' library(ComplexHeatmap, ggplot2, grid, circlize, openxlsx, grDevices, utils)
#' 
#' my_predictions = predict_lundtax(this_data = sjodahl_2017, 
#'                                  gene_id = "hgnc_symbol", 
#'                                  impute = TRUE, 
#'                                  adjust = TRUE)
#' 
#' plot_scores(this_data = my_predictions,
#'             out_path = "../",
#'             out_format = "pdf",
#'             title = "Lund2017")
#'}
#'
plot_scores = function(this_data = NULL,
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
  if(is.null(this_data) | !is.list(this_data)){
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
    sample_order <- order(this_data$scores$ProliferationScore)

    #set split (top annotation track of heatmap)
    if(!is.null(hm_split)){
      split = hm_split
      }else{
        split <- factor(this_data$predictions_7classes,
                        levels = c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))
     }

    #get immune names
    immune_names <- c("Immune141_UP", "B-cells", "T-cells",
                      "T-cells CD8+", "NK-cells",
                      "Cytotoxicity Score", "Neutrophils",
                      "Monocytic lineage", "Macrophages",
                      "M2 macrophage", "Myeloid Dendritic Cells")

    #get stromal names
    stromal_names <- c("Stromal141_UP", "Endothelial cells",
                       "Fibroblasts", "Smooth muscle")

    #set colours
    #proliferation
    col_fun_proliferation =
      circlize::colorRamp2(c(quantile(this_data$scores$ProliferationScore, 0.05),
                             median(this_data$scores$ProliferationScore),
                             quantile(this_data$scores$ProliferationScore, 0.95)),
                           c("#21908CFF","white", "#B63679FF"))

    #progression
    col_fun_progression =
      circlize::colorRamp2(c(quantile(this_data$scores$ProgressionScore, 0.05),
                             median(this_data$scores$ProgressionScore),
                             quantile(this_data$scores$ProgressionScore, 0.90)),
                           c("#FAEBDDFF","#A11A5BFF", "#4C1D4BFF"))

    #create colour object
    colour_obj = list(Lund = lund_colors$lund_colors,
                      ProliferationScore = col_fun_proliferation,
                      MolecularGradeWHO1999 = c("G1_2"="white","G3"="black"),
                      MolecularGradeWHO2016 = c("HG"="black","LG"="white"),
                      ProgressionScore = col_fun_progression,
                      ProgressionRisk = c("HR" = "#A11A5BFF", "LR" = "#FAEBDDFF"),
                      PossibleProstate = c("NO" = "#eae5eb", "YES" = "#ee82ee"))

    #plotting
    #build heatmap annotation (top)
    hm <- HeatmapAnnotation(Lund = split,
                            ProliferationScore = this_data$scores$ProliferationScore,Lund = split,
                            MolecularGradeWHO1999 = this_data$scores$MolecularGradeWHO1999,
                            MolecularGradeWHO2016 = this_data$scores$MolecularGradeWHO2016,
                            ProgressionScore = this_data$scores$ProgressionScore,
                            ProgressionRisk = this_data$scores$ProgressionRisk,
                           #PossibleProstate = this_data$scores$PossibleProstate,
                            annotation_name_side = "left",
                            show_legend = plot_anno_legend,
                            annotation_name_gp = gpar(fontsize = plot_font_size),
                            border = TRUE,
                            col = colour_obj,
                            column_title = NULL,
                            row_title = NULL)

    #immune scores heatmap
    hm_immune_scores <- Heatmap(t(scale(this_data$scores[,immune_names, drop = FALSE])),
                                top_annotation = hm,
                                column_order = sample_order,
                                height = unit(5*ncol(this_data$scores[,immune_names,drop = FALSE]), "mm"),
                                width = unit(0.5*nrow(this_data$scores[,immune_names,drop = FALSE]), "mm"),
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
    hm_stroma_scores <- Heatmap(t(scale(this_data$scores[,stromal_names, drop = FALSE])),
                                column_order = sample_order,
                                height = unit(5*ncol(this_data$scores[,stromal_names, drop = FALSE]), "mm"),
                                width = unit(0.5*nrow(this_data$scores[,stromal_names,drop = FALSE]), "mm"),
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
    my_scores = this_data$scores
    my_scores = as.data.frame(my_scores)

    #add new column with cohort information
    my_scores$Cohort = title

    #convert correct columns to factors
    factor_cols <- c('MolecularGradeWHO1999' ,'MolecularGradeWHO2016',
                     'ProgressionRisk', 'PossibleProstate', 'Cohort')

    my_scores[,factor_cols] <- lapply(my_scores[,factor_cols] , factor)

    if(to_xlsx){
      my_scores <- tibble::rownames_to_column(my_scores, "SampleID")
      write.xlsx(my_scores, paste0(out_path, title, "_scores.xlsx"))
      if(verbose){
        message(paste0("Prediction scores exported as xlsx to ", out_path, title, "_scores.xlsx"))
        }
      }

    #return data frame
    return(my_scores)

    if(verbose){
      message("Prediction Scores for each sample and class sucessfully generated!")
      }
    }
  }
