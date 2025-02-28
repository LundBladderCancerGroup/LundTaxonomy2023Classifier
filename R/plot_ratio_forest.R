#' @title Forest Plot
#'
#' @description Construct a forest plot using Cox model or generalized linnear model for a set of 
#' signature scores.
#'
#' @details Depending on the user defines `stat_plot`, this function internally calls `get_survival` 
#' or `get_glm`, which internally calls `int_prediction_wrangler`. The user can also provide the 
#' plotting data with `this_data` directly. If `plot_data` is set to `odds_ratio` the function will 
#' visualize odds ratio for one selected categorical variable. The user sets this variable with `categorical_factor`.
#' This should be a valid column name from either the metadata provided with `these_samples_metadata` 
#' or from the signature score data subset. The user can then set what numerical columns that are to 
#' be tested for association. This is controlled with `predictor_columns`. The input for this 
#' parameter should match column names corresponding to numeric variables in the metadata or signature 
#' score data subset. Alternatively, the function can also visualize hazard ratios calculated with 
#' `get_surv`. To do so, set `plot_data = "hazard_ratio"`. If so, the user also needs to point the 
#' function to the corresponding survival columns in the provided metadata. This is done with `surv_time` 
#' and `surv_event`. To return a plot for just one subtype, specify the subtype and class with 
#' `this_subtype` and `subtype_class`. To return a plot for more than one subtype, it is recommended 
#' to provide the plot data directly with `this_data`.
#'
#' @param these_predictions Required parameter if `this_data` is not provided. Should be output from 
#' [LundTax2023Classifier::lundtax_predict_sub()].
#' @param these_samples_metadata Required parameter if `this_data`is not provided. Metadata associated
#' with he prediction output. Also possible for the user to provide a metadata subset with samples 
#' of interest, the return will be restricted to the samples within the specified group
#' @param this_data Required parameter if `these_predictions` and `these_samples_metadata` is not 
#' provided. Should be the output from `get_survival` or `get_glm`, depending on the plot type set with `stat_plot`.
#' @param stat_plot Required parameter. Sets the plot type, for hazard ratio, set to `hazard_ratio`. 
#' For odds ratio, set to `odds_ratio`.
#' @param this_subtype Optional Specify subtype of interest. Leave as NULL to not separate statistics on subtype.
#' @param subtype_class Can be one of the following; 5_class or 7_class. Default is 5_class.
#' @param scale Optional parameter. A numeric value to scale the numeric scores. If provided, all 
#' numeric scores will be multiplied by this value.
#' @param bin_scores Boolean parameter. Set to TRUE to bin the numeric scores into discrete bins. Default is FALSE.
#' @param n_bins Optional parameter. The number of bins to use when binning numeric scores. Default is 10.
#' @param surv_time Required parameter if `stat_plot` is set to `hazard_ratio`, should be the name 
#' of the column in the metadata with survival time. Should be of value numeric.
#' @param surv_event Required parameter if `stat_plot` is set to `hazard_ratio`, should be the name 
#' of the column in the metadata with survival event. Should be of value factor, with two levels.
#' @param categorical_factor Required parameter if `this_glm`is not provided. This should be the categorical
#' variable that is intended for testing. In addition, this should also be a variable of type factor, 
#' with exactly 2 levels.
#' @param predictor_columns Optional, should be a vector with column names, either from the provided 
#' metadata or signature score object, to be tested for. If not provided, the function will subset 
#' data to the signature scores returned with `lundtax_predict_sub`.
#' @param exclude_columns Optional argument, specify columns you wish to exclude from the standard 
#' predictor columns. Note, this parameter is only validated if predictor_columns is NULL (default).
#' @param significant_p Numeric parameter for flagging significant p values. Default is 0.05.
#' @param sample_id_col Optional parameter. Allows the user to manually specify the name of a column with sample ID.
#' @param row_to_col Optional parameter, set to TRUE to convert row names in metadata to a new column 
#' called sample_id. Default is FALSE.
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can further specify the dimensions of the returned plot with `plot_width` and `plot_height`.
#' @param file_name Optional, if plot is being saved to disk, specify the file name for the file. 
#' @param plot_title Title for plot.
#' @param plot_subtitle Subtitle for plot.
#' @param plot_subtitle Caption for plot.
#' @param plot_width This parameter controls the width in inches.  Default is 8 (2400 pixels at 300 PPI).
#' @param plot_height This parameter controls the height in inches. Default is 8 (2400 pixels at 300 PPI).
#' @param plot_order Optional parameter for setting the order of scores in the returned plot. 
#' Only applies if `arrange_plot` is set to TRUE.
#' @param plot_arrange Boolean parameter, if set to TRUE the user can specify the order of the score 
#' levels on the y axis. Default is FALSE.
#' @param return_data Boolean parameter, set to TRUE and return the formatted data used for plotting. 
#' Default is FALSE
#'
#' @return A forest plot as grub object.
#'
#' @import ggplot2 dplyr
#'
#' @export
#'
#' @examples
#' #load pacakges
#' library(dplyr, ggplot2)
#' 
#' #get prediction calls
#' sjodahl_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                         impute = TRUE)
#' 
#' #hazard ratio                                      
#' plot_ratio_forest(these_predictions = sjodahl_predicted,
#'                   these_samples_metadata = sjodahl_2017_meta,
#'                   stat_plot = "hazard_ratio",
#'                   plot_subtitle = "n Samples: 267",
#'                   plot_title = "All 7 class samples",
#'                   subtype_class = "7_class", 
#'                   this_subtype = NULL,
#'                   surv_time = "surv_css_time",
#'                   surv_event = "surv_css_event")
#'
#' #odds ratio
#' plot_ratio_forest(these_predictions = sjodahl_predicted,
#'                   these_samples_metadata = sjodahl_2017_meta,
#'                   stat_plot = "odds_ratio",
#'                   plot_title = "Uro Samples - Adjuvant Chemo",
#'                   plot_subtitle = "n Samples: 121",
#'                   subtype_class = "5_class", 
#'                   this_subtype = "Uro",
#'                   predictor_columns = c("progression_score", 
#'                                         "proliferation_score", 
#'                                         "monocytic_lineage"),
#'                   categorical_factor = "adj_chemo")
#' 
plot_ratio_forest = function(these_predictions = NULL,
                             these_samples_metadata = NULL,
                             this_data = NULL,
                             stat_plot,
                             all_subs = FALSE,
                             this_subtype = NULL,
                             subtype_class = "5_class",
                             scale = NULL,
                             bin_scores = FALSE,
                             n_bins = 10,
                             surv_time = NULL,
                             surv_event = NULL,
                             categorical_factor = NULL,
                             predictor_columns = NULL,
                             exclude_columns = NULL,
                             significant_p = 0.05,
                             sample_id_col = NULL, 
                             row_to_col = FALSE,
                             out_path = NULL,
                             out_format = "png",
                             file_name = "my_plot",
                             plot_title = "My Plot",
                             plot_subtitle = "",
                             plot_caption = "",
                             plot_width = 8,
                             plot_height = 8,
                             plot_order = NULL,
                             plot_arrange = TRUE,
                             return_data = FALSE){
  #checks
  if(missing(stat_plot)){
    stop("User must specify stat_plot, can be one of the following; hazard_ratio or odds_ratio...")
  }
  
  if(length(this_subtype) > 1){
    stop("If you want more than one subtype (but not all), it's recommended to provide these as this_surv")
  }
  
  if(!is.null(this_subtype)){
    if(!this_subtype %in% names(lund_colors$lund_colors)){
      stop("Please check spelling of subtype...")
    }
  }
  
  #run get_survival if user has provided prediction data and not survival object
  #if not all subtypes are requested
  if(!is.null(these_predictions) && is.null(this_data)){
    if(stat_plot == "hazard_ratio"){
      
      this_data = get_survival(these_predictions = these_predictions,
                               these_samples_metadata = these_samples_metadata,
                               subtype_class = subtype_class,
                               predictor_columns = predictor_columns,
                               exclude_columns = exclude_columns, 
                               surv_time = surv_time, 
                               surv_event = surv_event,
                               this_subtype = this_subtype,
                               scale = scale, 
                               bin_scores = bin_scores, 
                               n_bins = n_bins,
                               sample_id_col = sample_id_col, 
                               row_to_col = row_to_col)
      
    }else if(stat_plot == "odds_ratio"){
      
      this_data = get_glm(these_predictions = these_predictions,
                          these_samples_metadata = these_samples_metadata,
                          subtype_class = subtype_class,
                          categorical_factor = categorical_factor,
                          predictor_columns = predictor_columns,
                          exclude_columns = exclude_columns, 
                          this_subtype = this_subtype,
                          scale = scale, 
                          bin_scores = bin_scores, 
                          n_bins = n_bins,
                          sample_id_col = sample_id_col, 
                          row_to_col = row_to_col)
      
    }else{
      stop("stat_plot must be one of the following; hazard_ratio or odds_ratio...")
    }
  }else if(!is.null(these_predictions) && !is.null(this_data)){
    message("Both predictions and plot data are provided, function will plot the provided data...")
  }else if(is.null(these_predictions) && !is.null(this_data)){
    
  }
  
  if(return_data){
    message("No plot generated, returning plot data instead...")
    return(this_data)
  }
  
  #annotate significance of p value
  this_data$significant = ifelse(this_data$p_value < significant_p, "significant", "not significant")
  
  #rename column names for generalization
  names(this_data)[3] = "ratio"
  names(this_data)[4] = "conf_2.5"
  names(this_data)[5] = "conf_97.5"
  
  if(plot_arrange){
    # Define the mapping of old names to new names
    score_name_mapping <- c(
      "proliferation_score" = "Proliferation Score",
      "molecular_grade_who_1999_score" = "Mol. grade (WHO1999)",
      "molecular_grade_who_2022_score" = "Mol. grade (WHO2022)",
      "progression_score" = "Progression Score",
      "immune141_up" = "Immune 141_UP",
      "b_cells" = "B Cells",
      "t_cells" = "T Cells",
      "t_cells_cd8" = "CD8+ T Cells",
      "nk_cells" = "NK Cells",
      "cytotoxicity_score" = "Cytotoxicity Score",
      "neutrophils" = "Neutrophils",
      "monocytic_lineage" = "Monocytic Lineage",
      "macrophages" = "Macrophages",
      "m2_macrophage" = "M2 Macrophages",
      "myeloid_dendritic_cells" = "Myeloid DCs",
      "stromal141_up" = "Stromal 141_UP",
      "endothelial_cells" = "Endothelial Cells",
      "fibroblasts" = "Fibroblasts",
      "smooth_muscle" = "Smooth Muscle"
    )
    
    # Update the score column in this_data
    this_data$score <- score_name_mapping[this_data$score]
    
    if(is.null(plot_order)){
      #create the desired order as a vector
      desired_order <- c("Proliferation Score", "Mol. grade (WHO1999)", "Mol. grade (WHO2022)", "Progression Score", 
                         "Immune 141_UP", "NK Cells", "T Cells", "CD8+ T Cells", "Cytotoxicity Score", "B Cells", 
                         "Myeloid DCs", "Monocytic Lineage", "Macrophages", "M2 Macrophages", "Neutrophils", 
                         "Stromal 141_UP", "Fibroblasts", "Endothelial Cells", "Smooth Muscle")
    }else{
      desired_order = plot_order
    }
    
    # Convert the score column to a factor with the specified levels
    this_data$score <- factor(this_data$score, levels = rev(desired_order))
  }
  
  #build plot
  my_plot = ggplot(this_data, aes(x = ratio, y = score)) +
    geom_point(aes(color = significant), size = 3) +
    geom_errorbarh(aes(xmin = conf_2.5, xmax = conf_97.5, color = significant), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(title = plot_title, subtitle = plot_subtitle,  caption = plot_caption, x = stat_plot, y = "") +
    scale_color_manual(values = c("significant" = "red", "not significant" = "black")) +
    theme(legend.position = "none",
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(color = "black", size = 10),
          axis.ticks.x = element_line(linewidth = 0.4),
          axis.ticks.y = element_line(linewidth = 0.4),
          plot.title = element_text(),
          plot.caption = element_text(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
          axis.line.x = element_blank())
  
  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, file_name, "_", stat_plot, "_forest.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, file_name, "_", stat_plot, "_forest.png"),
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
    message(paste0("Plot exported to ", out_path, file_name, "_", stat_plot, "_forest.", out_format))
  }else{
    return(my_plot) 
  }
}
