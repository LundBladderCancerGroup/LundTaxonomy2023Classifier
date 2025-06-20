% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/int_ratio_score.R
\name{int_ratio_score}
\alias{int_ratio_score}
\title{Calculate Scores.}
\usage{
int_ratio_score(
  this_data = NULL,
  variable = NULL,
  gene_id = "hgnc_symbol",
  verbose = TRUE
)
}
\arguments{
\item{this_data}{Required parameter. Data frame or matrix with expression values.}

\item{variable}{Required parameter. Input should be one of the following;
proliferation, or progression.}

\item{gene_id}{Specify the type of gene identifier used in \code{this_data}.
Accepted values are; hgnc_symbol (default) or ensembl_gene_id.}

\item{verbose}{A logical value indicating whether processing messages will be
printed or not. Default is TRUE.}
}
\value{
A list with two objects. 1, A data frame with scores for the selected variable.
2, A data frame indicating what genes from the incoming data are missing, based on the expected
genes for signature calculations.
}
\description{
Calculate proliferation and progression scores.
}
\details{
Internal function called by \code{\link[=lundtax_calc_sigscore]{lundtax_calc_sigscore()}}.
Not meant for out of package use. Takes a data frame of matrix with expression values and
calculates scores based on gene expression.
}
