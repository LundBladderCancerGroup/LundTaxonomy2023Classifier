#' @title Violin Plot Subtype Scores
#'
#' @description Visualize subtype prediction score for a set subtype in a violin plot.
#'
#' @details Take the output from `lundtax_predict_sub` and return a violin plot where each violin 
#' represent the distribution of the, for that class, subtype prediction score. Set the subtype of 
#' desire with `this_subtype`. The subtype can be one of the subtypes included in the LundTax subtype 
#' classification nomenclature. Parameters starting with the plot prefix are available for the user 
#' to gain more control of the returned plot. The user can also return the wrangled data used in the 
#' plot by setting `return_data = TRUE`.
#'
#' @param these_predictions Required parameter, should be the output from 
#' [LundTax2023Classifier::lundtax_predict_sub()].
#' @param this_subtype Required parameter. Should be one of the set subtype classes from the
#'  LundTax nomenclature.
#' @param plot_title Required parameter, if `out_path` is specified. plot title, will also be pasted to 
#' the exported file.
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can further specify the dimensions of the returned plot with `plot_width` and `plot_height`.
#' @param plot_width This parameter controls the width in inches. 
#' Default is 4 (1200 pixels at 300 PPI).
#' @param plot_height This parameter controls the height in inches. 
#' Default is 4 (1200 pixels at 300 PPI).
#' @param plot_adjust A multiplicate bandwidth adjustment. This makes it possible to adjust the 
#' bandwidth while still using the a bandwidth estimator. For example, adjust = 1/2 means use half 
#' of the default bandwidth.
#' @param plot_scale if "area" (default), all violins have the same area (before trimming the tails). 
#' If "count", areas are scaled proportionally to the number of observations. If "width", all 
#' violins have the same maximum width.
#' @param plot_trim If TRUE (default), trim the tails of the violins to the range of the data. If
#' FALSE, don't trim the tails.
#' @param return_data Set to TRUE to return tidy data used by the plotting function. Default is FALSE.
#'
#' @return Nothing.
#' 
#' @import ggplot2 dplyr reshape2
#' @rawNamespace import(gridExtra, except = combine)
#'
#' @export
#'
#' @examples
#' my_predictions = lundtax_predict_sub(these_predictions = sjodahl_2017, 
#'                                      gene_id = "hgnc_symbol", 
#'                                      impute = TRUE, 
#'                                      adjust = TRUE)
#'                                      
#' uro_scores = plot_subscore_violin(these_predictions = my_predictions, 
#'                                   this_subtype = "Uro")
#'                                   
plot_subscore_violin = function(these_predictions, 
                                this_subtype,
                                plot_title = NULL,
                                out_path = NULL,
                                out_format = "png",
                                plot_width = 4,
                                plot_height = 4,
                                plot_adjust = 2,
                                plot_scale = "width",
                                plot_trim = FALSE,
                                return_data = FALSE){
  
  #subset scores
  these_scores = data.frame(these_predictions$subtype_scores)
  
  #get subtypes based on subtype class
  if(this_subtype %in% c("Uro", "GU", "BaSq", "Mes", "ScNE")){
    these_subtypes = data.frame(these_predictions$predictions_5classes)
    subtype_classes = c("Uro", "GU", "BaSq", "Mes", "ScNE")
  }else if(this_subtype %in% c("UroA", "UroB", "UroC")){
    these_subtypes = data.frame(these_predictions$predictions_7classes)
    subtype_classes = c("UroA", "UroB", "UroC")
  }

  #convert rowname to columns
  these_scores = tibble::rownames_to_column(these_scores, "sample_id")
  these_subtypes = tibble::rownames_to_column(these_subtypes, "sample_id")

  #rename column
  colnames(these_subtypes)[2] = "subtype"

  #subset consensus subtypes
  these_samples = dplyr::filter(these_subtypes, subtype == this_subtype)

  #replce NAs with zero
  these_scores[is.na(these_scores)] <- 0
  
  #melt data frame
  melted_scores = melt(data = these_scores,
                       id = "sample_id",
                       variable.name = "subtype",
                       measure.vars = c("Uro",
                                        "UroA",
                                        "UroB",
                                        "UroC",
                                        "GU",
                                        "BaSq",
                                        "Mes",
                                        "ScNE"))
  
  #convert to factor
  melted_scores$subtype = as.factor(melted_scores$subtype) 
  
  #subset on class
  class_melted = melted_scores %>%  
    dplyr::filter(subtype %in% subtype_classes)

  #get sample IDs for each subtype
  this_melted = dplyr::filter(class_melted, sample_id %in% these_samples$sample_id)
  
  #get colour for the subtype to be plotted
  subtype_col = lund_colors$lund_colors[this_subtype]
  
  if(return_data){
    message("No plot generated, returning tidy data instead...")
    return(this_melted)
  }

  #plot subtype scores
  my_plot = this_melted %>% 
    arrange(value) %>%    
    mutate(subtype = factor(subtype, levels = c("Uro", "GU", "BaSq", "Mes", "ScNE", "UroA", "UroB", "UroC"))) %>% 
    ggplot(aes(x = subtype, y = value, fill = subtype), show.legend = TRUE) + 
    geom_violin(scale = plot_scale, trim = plot_trim, color = NA, adjust = plot_adjust) +
    scale_fill_manual(values = lund_colors$lund_colors, drop = FALSE) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    ggtitle(label = plot_title) +
    ylab(this_subtype) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks = seq(0, 1, by = 0.5)) +
    theme(legend.position = "none", 
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
          axis.line.x = element_blank(), 
          axis.title.y = element_blank())
  
  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, plot_title, "_subscore_viol.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, plot_title, "_subscore_viol.png"),
          width = plot_width,
          height = plot_height,
          units = "in",
          res = 300,
          pointsize = 10,
          bg = "white")
    }else{
      stop("Enter a valid output format (pdf or png)...")
    }
    print(my_plot)
    dev.off()
    message(paste0("Plot exported to ", out_path, plot_title, "_subscore_viol.", out_format))
  }else{
    return(my_plot) 
  }
}
