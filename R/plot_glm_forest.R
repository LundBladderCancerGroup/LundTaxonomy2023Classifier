#' @title GLM Forest Plot.
#'
#' @description Construct a forest plot using GLM data for a set of signature scores.
#'
#' @details This function internally calls `get_glm_scores`, which internally calls `int_prediction_wrangler`. 
#' This happens if the user does not call the plotting function using the `this_glm` parameter. 
#' Many of the parameters in this function are recycled in the internal function calls to improve 
#' user flexibility. It's possible to return a forest plot for all subtypes within the specified 
#' `subtype_class`. To do so, set `all_subs = TRUE`. This will trigger a loop for which `get_glm_scores` 
#' is internally called for each subtype. To return a plot for just one subtype, set `all_subs = FALSE` 
#' and specify the subtype and class with `this_subtype` and `subtype_class`. To return a plot for more
#' than one subtype, but not all, call get_glm_score for each iteration of subtypes, see last example
#'  in the docs.
#'
#' @param these_predictions Required parameter if `this_glm` is not provided. Should be output from 
#' [LundTax2023Classifier::lundtax_predict_sub()].
#' @param these_samples_metadata Required parameter if `this_glm`is not provided. Metadata associated
#' with he prediction output. Also possible for the user to provide a metadata subset with samples 
#' of interest, the return will be restricted to the samples within the specified group
#' @param this_glm Required parameter if `these_predictions` and `these_samples_metadata` is not 
#' provided (The output from `get_glm_scores`).
#' @param plot_title Title for plot.
#' @param all_subs Boolean, default is FALSE. Set to TRUE to return a forest plot for all subtypes 
#' within the specified class (`subtype_class`).
#' @param subtype_class Can be one of the following; 5_class or 7_class. Default is 5_class.
#' @param cat_variable Required parameter if `this_glm`is not provided. This should be the categorical
#' variable that is intended for testing. In addition, this should also be a variable of type factor, 
#' with exactly 2 levels.
#' @param this_subtype Required parameter if `this_glm` is not provided, see `get_glm_scores`.
#' @param sample_id_col Required parameter if `this_glm` is not provided, see `get_glm_scores`.
#' @param row_to_col Required parameter if `this_glm` is not provided, see `get_glm_scores`.
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can further specify the dimensions of the returned plot with `plot_width` and `plot_height`.
#' @param plot_width This parameter controls the width in inches.  Default is 8 (2400 pixels at 300 PPI).
#' @param plot_height This parameter controls the height in inches. Default is 8 (2400 pixels at 300 PPI).
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
#' #build plot                                        
#' plot_glm_forest(these_predictions = sjodahl_predicted,
#'                 these_samples_metadata = sjodahl_2017_meta,
#'                 plot_title = "Uro, adj_chemo",
#'                 all_subs = FALSE,
#'                 subtype_class = "5_class", 
#'                 this_subtype = "Uro",
#'                 cat_variable = "adj_chemo")
#' 
plot_glm_forest = function(these_predictions = NULL,
                           these_samples_metadata = NULL,
                           this_glm = NULL,
                           plot_title = "My Plot",
                           all_subs = FALSE,
                           subtype_class = "5_class",
                           cat_variable = NULL,
                           this_subtype = NULL,
                           sample_id_col = NULL, 
                           row_to_col = FALSE,
                           out_path = NULL,
                           out_format = "png",
                           plot_width = 8,
                           plot_height = 8,
                           return_data = FALSE){
  #checks
  if(length(this_subtype) > 1){
    stop("If you want more than one subtype (but not all), it's recommended to provide these as this_glm")
  }
  
  if(!is.null(this_subtype)){
    if(!this_subtype %in% names(lund_colors$lund_colors)){
      stop("Please check spelling of subtype...")
    }
  }

  #run get_glm if user has provided prediction data and not glm object
  #if not all subtypes are requested
  if(!is.null(these_predictions) && is.null(this_glm)){
    if(!all_subs){
      this_glm = get_glm_scores(these_predictions = these_predictions,
                                these_samples_metadata = these_samples_metadata,
                                subtype_class = subtype_class,
                                cat_variable = cat_variable,
                                this_subtype = this_subtype,
                                sample_id_col = sample_id_col, 
                                row_to_col = row_to_col)
    }else{
      if(subtype_class == "5_class"){
        uro_glms = get_glm_scores(these_predictions = these_predictions,
                                  these_samples_metadata = these_samples_metadata,
                                  subtype_class = "5_class",
                                  cat_variable = cat_variable,
                                  this_subtype = "Uro",
                                  sample_id_col = sample_id_col, 
                                  row_to_col = row_to_col)

      }else if(subtype_class == "7_class"){
        uroa_glms = get_glm_scores(these_predictions = these_predictions,
                                   these_samples_metadata = these_samples_metadata,
                                   subtype_class = subtype_class,
                                   cat_variable = cat_variable,
                                   this_subtype = "UroA",
                                   sample_id_col = sample_id_col, 
                                   row_to_col = row_to_col)

        urob_glms = get_glm_scores(these_predictions = these_predictions,
                                   these_samples_metadata = these_samples_metadata,
                                   subtype_class = subtype_class,
                                   cat_variable = cat_variable,
                                   this_subtype = "UroB",
                                   sample_id_col = sample_id_col, 
                                   row_to_col = row_to_col)

        uroc_glms = get_glm_scores(these_predictions = these_predictions,
                                   these_samples_metadata = these_samples_metadata,
                                   subtype_class = subtype_class,
                                   cat_variable = cat_variable,
                                   this_subtype = "UroC",
                                   sample_id_col = sample_id_col, 
                                   row_to_col = row_to_col)
      }

      gu_glms = get_glm_scores(these_predictions = these_predictions,
                               these_samples_metadata = these_samples_metadata,
                               subtype_class = subtype_class,
                               cat_variable = cat_variable,
                               this_subtype = "GU",
                               sample_id_col = sample_id_col, 
                               row_to_col = row_to_col)

      basq_glms = get_glm_scores(these_predictions = these_predictions,
                                 these_samples_metadata = these_samples_metadata,
                                 subtype_class = subtype_class,
                                 cat_variable = cat_variable,
                                 this_subtype = "BaSq",
                                 sample_id_col = sample_id_col, 
                                 row_to_col = row_to_col)

      scne_glms = get_glm_scores(these_predictions = these_predictions,
                                 these_samples_metadata = these_samples_metadata,
                                 subtype_class = subtype_class,
                                 cat_variable = cat_variable,
                                 this_subtype = "ScNE",
                                 sample_id_col = sample_id_col, 
                                 row_to_col = row_to_col)

      mes_glms = get_glm_scores(these_predictions = these_predictions,
                                these_samples_metadata = these_samples_metadata,
                                subtype_class = subtype_class,
                                cat_variable = cat_variable,
                                this_subtype = "Mes",
                                sample_id_col = sample_id_col, 
                                row_to_col = row_to_col)

      if(subtype_class == "5_class"){
        this_glm = rbind(uro_glms, gu_glms, basq_glms, scne_glms, mes_glms)
      }else if(subtype_class == "7_class"){
        this_glm = rbind(uroa_glms, urob_glms, uroc_glms, gu_glms, basq_glms, scne_glms, mes_glms)
      }
    }
  }else if(is.null(these_predictions) && !is.null(this_glm)){
    message("Both predictions and GLM object are provided, function will plot the provided GLM...")
  }
  
  if(return_data){
    message("No plot generated, returning plot data instead...")
    return(this_glm)
  }

  #build plot
  my_plot = ggplot(data = this_glm, aes(x = score, y = odds_ratio, ymin = conf_2.5, ymax = conf_97.5, color = subtype)) +
    geom_hline(yintercept = 1, lty = 2) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    coord_flip() +
    xlab("Signature") +
    ylab("Hazard Ratio (97.5 CI)") +
    theme_bw() +
    scale_color_manual(values = lund_colors$lund_colors) +
    labs(title = plot_title) +
    theme(legend.title = element_blank())

  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, plot_title, "_forest_plot.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, plot_title, "_forest_plot.png"),
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
    message(paste0("Plot exported to ", out_path, plot_title, "_ranked_score.", out_format))
  }else{
    return(my_plot) 
  }
}
