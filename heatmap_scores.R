#' @title Heatmap Scores.
#'
#' @description Build a heatmap with scores retrieved with the predict_LundTax2023 function.
#'
#' @details Construct and export (pdf or png) a highly customizable heatmap visualizing scores
#' predicted with `predict_LundTax2023`. This function depends on Complexheatmap.
#' For a greater detail of what the function can do, see examples and parameter descriptions.
#'
#' @param this_data A list with a data frame object called scores. Returned with `predict_LundTax2023`.
#' @param out_path Optional, set path to export plot. If not provided, tidy version of incoming
#' scores in data frame format will be returned.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' if pdf, the returned pdf will be in A4 (horizontal) format, if png is specified the user can
#' control the dimensions, see `plot_width` and `plot_height`.
#' @param title Heatmap title, will also be pasted to the exported file.3
#' @param hm_split Optional parameter for controlling how the data is split into different groups.
#' If not provided, the function will split on subtype.
#' @param hm_cluster Boolean parameter, set to TRUE to cluster the rows (default is FALSE).
#' @param plot_anno_legend Expects a vector with TRUE/FALSE (7 in total), thsi decides what legends
#' will be on the final heatmap. Default is to only show the legend for the subtypes.
#' @param plot_hm_legend Boolean parameter. Set to TRUE to show heatmap legend. Default is FALSE.
#' @param plot_width If `out_format = "png"`, this parameter controls the width in pixels. Default is 4300 pixels
#' @param plot_height If `out_format = "png"`, this parameter controls the height in pixels. Default is 1600 pixels
#' @param plot_font_size Optional parameter to control the size of the font in the generated heatmap.
#' Note, the title of the plot will always be twice that of the set value here (default = 12).
#'
#' @return Nothing
#'
#' @import ComplexHeatmap ggplot2 circlize
#'
#' @export
#'
#' @examples
#' results = predict_LundTax2023(data = Lund2017,
#'                              include_data = TRUE,
#'                              adjust = FALSE,
#'                              impute = TRUE,
#'                              gene_id = "hgnc_symbol")
#'
#' heatmap_scores(this_data = results,
#'                out_path = "../",
#'                out_format = "pdf",
#'                title = "Lund2017")
#'
heatmap_scores = function(this_data = NULL,
                          out_path = NULL,
                          out_format = "png",
                          title = NULL,
                          hm_split = NULL,
                          hm_cluster = FALSE,
                          plot_anno_legend = NULL,
                          plot_hm_legend = FALSE,
                          plot_width = 4300,
                          plot_height = 1600,
                          plot_font_size = 12){

  if(is.null(plot_anno_legend)){
    plot_anno_legend = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  }else{
    plot_anno_legend = plot_anno_legend
    }

  #check the type of incoming data
  if(is.null(this_data) | !is.list(this_data)){
    stop("Input should be the result of applying predict_LundTax2023")
  }

  if(!is.null(out_path)){
    if(is.null(title)){
      stop("Please specify a title for your plot with `title`...")
    }

    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, title, "_heatmap_scores.pdf"),
          paper = "a4r")
    }else if(out_format == "png"){
        png(paste0(out_path, title, "heatmap_scores.png"),
            width = plot_width,
            height = plot_height,
            units = "px",
            res = 300,
            pointsize = 10,
            bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
      }
  }else{
    message("No output path provided, returning the score data frame in tidy format...")
    return(as.data.frame(this_data$scores)) #TODO: tidy...
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
                          PossibleProstate = this_data$scores$PossibleProstate,
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

  }

heatmap_scores(this_data = bowden_geTMM_scores_out,
               title = "NEWPLasdasdOTdsds",
               out_path = "../../figs_scratch/new/",
               plot_font_size = 13,
               out_format = "png",
               hm_cluster = FALSE,
               plot_hm_legend = FALSE)
