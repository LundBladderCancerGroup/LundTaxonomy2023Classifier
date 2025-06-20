% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_hm_scores.R
\name{plot_hm_scores}
\alias{plot_hm_scores}
\title{Plot Signature Scores.}
\usage{
plot_hm_scores(
  these_predictions = NULL,
  proportional_scores = FALSE,
  out_path = NULL,
  out_format = "png",
  return_scores = FALSE,
  to_xlsx = FALSE,
  plot_title = NULL,
  hm_split = NULL,
  subtype_annotation = "5_class",
  hm_cluster = FALSE,
  plot_anno_legend = NULL,
  plot_hm_legend = FALSE,
  plot_width = 14,
  plot_height = 6,
  plot_font_size = 12,
  verbose = TRUE
)
}
\arguments{
\item{these_predictions}{A list with a data frame object called scores. Returned with
\code{\link[=lundtax_predict_sub]{lundtax_predict_sub()}}.}

\item{proportional_scores}{Set to TRUE to transform signature scores into proportions, default is FALSE.}

\item{out_path}{Optional, set path to export plot. If not provided, tidy version of incoming
scores in data frame format will be returned (\code{return_scores} will be auto-defaulted to TRUE).}

\item{out_format}{Required parameter if \code{out_path} is specified. Can be "png" (default) or "pdf".
The user can further specify the dimensions of the returned plot with \code{plot_width} and \code{plot_height}.}

\item{return_scores}{Set to TRUE to return prediction scores in a tidy format. Default is FALSE.}

\item{to_xlsx}{Boolean parameter, set to TRUE to export score data frame in xlsx format. Default
is FALSE.
If set to TRUE, the spreadsheet will be saved to the same path as the heatmap.}

\item{plot_title}{Required parameter. Heatmap title, will also be pasted to the exported file(s) as
well as a new column in the scores data frame under cohort.}

\item{hm_split}{Optional parameter for controlling how the data is split into different groups.
If not provided, the function will split on what is specified within \code{subtype_annotation}.}

\item{subtype_annotation}{Can be one of the following; "5_class" (default) or "7_class"
annotation.}

\item{hm_cluster}{Boolean parameter, set to TRUE to cluster the rows (default is FALSE).}

\item{plot_anno_legend}{Expects a vector with TRUE/FALSE (7 in total), thsi decides what legends
will be on the final heatmap. Default is to only show the legend for the subtypes.}

\item{plot_hm_legend}{Boolean parameter. Set to TRUE to show heatmap legend. Default is FALSE.}

\item{plot_width}{This parameter controls the width in inches.
Default is 14(4200 pixels at 300 PPI).}

\item{plot_height}{This parameter controls the height in inches.
Default is 6(1800 pixels at 300 PPI)}

\item{plot_font_size}{Optional parameter to control the size of the font in the generated heatmap.
Note, the title of the plot will always be twice that of the set value here (default = 12).}

\item{verbose}{Set to TRUE for debugging purposes. Default is FALSE.}
}
\value{
Data frame with prediction score for each sample and class, if return_scores = TRUE.
Otherwise, nothing.
}
\description{
Build a heatmap with scores retrieved with the
\code{\link[=lundtax_predict_sub]{lundtax_predict_sub()}} function.
}
\details{
Construct and export (pdf or png) a highly customizable heatmap visualizing prediction
scores for each sample and class, predicted with \code{\link[=lundtax_predict_sub]{lundtax_predict_sub()}}.
This function depends on Complexheatmap. It is also possible to return a data frame with
prediction scores in a tidy format. To do so, set \code{return_scores = TRUE}. For a greater explanation
on how to use the function, see parameter descriptions and examples.
}
\examples{
\dontrun{
my_predictions = lundtax_predict_sub(these_predictions = sjodahl_2017, 
                                     gene_id = "hgnc_symbol", 
                                     impute = TRUE, 
                                     adjust = TRUE)

plot_hm_scores(these_predictions = my_predictions,
               out_path = "../",
               out_format = "pdf",
               plot_title = "Lund2017")
}

}
