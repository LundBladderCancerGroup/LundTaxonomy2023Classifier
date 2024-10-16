#' @title Plot Ranked Score
#'
#' @description Return a point plot with ranked scores for a set variable.
#'
#' @details This function takes predictions returned with 
#' [LundTax2023Classifier::lundtax_predict_sub()], together with a score variable,
#' subtype class and returns a point plot (grob) with ranked scores along the x axis and 
#' scores along the y axis. The returned plot will also color fill the points based on subtype 
#' classification.
#'
#' @param these_predictions Required parameter, if no data provided with `this_data`.
#' Predictions returned with [LundTax2023Classifier::lundtax_predict_sub()].
#' @param this_data Optional parameter, makes it possible for the user to give the plotting function 
#' their own data, preferably retrieved with this function with return_data = TRUE. If provided, 
#' the plotting function will disregard any other parameters and only draw the plot using the data 
#' provided.
#' @param this_score Required parameter, should be one of the numeric columns in the scores data 
#' frame from the list returned with [LundTax2023Classifier::lundtax_predict_sub()].
#' @param this_subtype Optional parameter, lets the user subset the return to a specific subtype.
#' If nopt provided, all subtypes within the selected class will be returned.
#' @param subtype_class Required, one of the following; 5_class or 7_class. Needed for coloring the 
#' points based on subtype classification.
#' @param add_stat Boolean parameter, set to TRUE to add regression lines, p value for each regression. Default is FALSE.
#' @param plot_title Required parameter, if `out_path` is specified. plot title, will also be pasted to 
#' the exported file.
#' @param out_path Optional, set path to export plot.
#' @param out_format Required parameter if `out_path` is specified. Can be "png" (default) or "pdf".
#' The user can further specify the dimensions of the returned plot with `plot_width` and `plot_height`.
#' @param plot_width This parameter controls the width in inches. 
#' Default is 4 (1200 pixels at 300 PPI).
#' @param plot_height This parameter controls the height in inches. 
#' Default is 4 (1200 pixels at 300 PPI).
#' @param return_data Boolean parameter, set to TRUE and return the formatted data used for 
#' plotting. Default is FALSE
#'
#' @return Nothing.
#' 
#' @import ggplot2 dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom ggpubr stat_cor
#'
#' @export
#'
#' @examples
#' sjodahl_2017_predicted = lundtax_predict_sub(this_data = sjodahl_2017, 
#'                                              impute = TRUE)
#'
#' plot_ranked_score(these_predictions = sjodahl_2017_predicted, 
#'                   this_score = "proliferation_score", 
#'                   subtype_class = "7_class", 
#'                   plot_title = "Sjodahl 2017")
#'
plot_ranked_score = function(these_predictions = NULL, 
                             this_data = NULL,
                             this_score = NULL, 
                             this_subtype = NULL,
                             subtype_class = NULL, 
                             add_stat = FALSE,
                             title = NULL,
                             out_path = NULL,
                             out_format = "png",
                             plot_width = 4,
                             plot_height = 4,
                             plot_title = NULL,
                             return_data = FALSE){
  
  if(is.null(this_data)){
    #checks
    if(is.null(these_predictions)){
      stop("these_predictions is missing, should be the return from lund_predict_sub...")
    }
    
    #impute plot title
    if(is.null(plot_title)){
      plot_title = this_score
    }
    
    if(is.null(this_score)){
      stop("this_score is missing, should be a numeric variable in the scores data frame within these_predictions...")
    }else if(!this_score %in% colnames(these_predictions$scores)){
      stop(paste0(this_score, " is non-existing in the scores data frame within these_predictions..."))
    }else if(!is.numeric(these_predictions$scores[,this_score])){
      stop(paste0(this_score, " is not of type numeric..."))
    }
    
    #subset on subtype
    if(subtype_class == "5_class"){
      subtypes = data.frame(these_predictions$predictions_5classes) %>% 
        tibble::rownames_to_column("sample_id")
      if(is.null(this_subtype)){
        this_subtype = c("Uro", "GU", "ScNE", "BaSq", "Mes")
      }else{
        if(!this_subtype %in% c("Uro", "GU", "ScNE", "BaSq", "Mes")){
          stop("The requested subtype is not available in the selected class..")
        }
      }
    }else if(subtype_class == "7_class"){
      subtypes = data.frame(these_predictions$predictions_7classes) %>% 
        tibble::rownames_to_column("sample_id") 
      if(is.null(this_subtype)){
        this_subtype = c("UroA","UroB","UroC", "GU", "ScNE", "BaSq", "Mes")
      }else{
        if(!this_subtype %in% c("UroA","UroB","UroC", "GU", "ScNE", "BaSq", "Mes")){
          stop("The requested subtype is not available in the selected class..")
        }
      }
    }else{
      stop("A valid subtype class must be provided, one of the following; 5_class or 7_class...")
    }
    
    #subset score, left join with subtypes, add rank
    this_data = these_predictions$scores %>% 
      dplyr::select(all_of(this_score)) %>% 
      tibble::rownames_to_column("sample_id") %>% 
      dplyr::left_join(subtypes, by = "sample_id") %>% 
      dplyr::rename(score = 2, subtype = 3) %>%
      dplyr::filter(subtype %in% this_subtype) %>% 
      dplyr::mutate_at(vars(subtype), factor) %>%
      mutate(rank = as.numeric(as.factor(score)))
    
    if(return_data){
      message("No plot generated, returning plot data instead...")
      return(this_data)
    }
  }else{
    
    if(!is.null(this_subtype)){
      message("You have provided your own data, the function will not subset plotting data on the requested subtype...")
    }
    
    if(is.null(this_score)){
      stop("You must specify the type of score with `this_score` in your provided dataset...")
    }
  }
  
  #draw plot
  my_plot = ggplot(data = this_data, aes(x = rank, y = score, color = subtype)) +
    geom_point(size = 2, shape = 16) +
    {if(add_stat) geom_smooth(method = lm, se = FALSE)} +
    {if(add_stat) stat_cor(method = "pearson", label.x = 20)} +
    scale_color_manual(values = lund_colors$lund_colors) + 
    xlab("Index") + 
    ylab(this_score) +
    theme(legend.position = "none", 
          axis.text.y = element_text(color = "black", size = 7),
          axis.ticks.x = element_line(linewidth = 0.4),
          axis.ticks.y = element_line(linewidth = 0.4),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
          axis.line.x = element_blank(), 
          axis.title.y = element_text(angle = 90, colour = "black", size = 10, vjust = 2))
  
  if(!is.null(out_path)){
    #set PDF outputs
    if(out_format == "pdf"){
      pdf(paste0(out_path, title, "_ranked_score.pdf"),
          width = plot_width,
          height = plot_height)
      #set PNG outputs
    }else if(out_format == "png"){
      png(paste0(out_path, title, "_ranked_score.png"),
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
    message(paste0("Plot exported to ", out_path, title, "_ranked_score.", out_format))
  }else{
    return(my_plot) 
  }
}
