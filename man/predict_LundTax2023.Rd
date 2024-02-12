\name{predict_LundTax2023}
\alias{predict_LundTax2023}
\title{Predict Lund Taxonomy subtypes based on rule-based Random Forest classifiers}
\usage{
predict_LundTax2023(data, include_data = FALSE, include_scores = TRUE, ...)
}
\arguments{
\item{data}{matrix, data frame or multiclassPairs_object of gene expression values}

\item{include_data}{include data in output (disabled by default)}

\item{include_scores}{include prediciton scores for each sample and class in output (default)}

\item{gene_id}{specify the type of gene identifier used in the data. Default value is hgnc_symbol}
\itemize{
\item "hgnc_symbol" for HUGO gene symbols
\item "ensembl_gene_id" for Ensembl gene IDs
\item "entrezgene" for Entrez IDs
}

\item{...}{Additional parameters to be passed to the predict_RF function. If genes are missing in the data, include impute = TRUE here}
}
\value{
Returns a list object including:
\itemize{
\item Data (optional, not included by default)
\item Prediction scores for all classes (optional, included by default)
\item Predicted LundTax class for 7-class system
\item Predicted LundTax class for 5-class system
}
}
\description{
Predict Lund Taxonomy subtypes based on rule-based Random Forest classifiers
}
\details{
This function uses 2 classifiers to classify the samples: 5-class classifier first  classifies samples into Uro, GU, BaSq, Mes or ScNE.
Samples classified as Uro receive a second classification as UroA, B or C by the second classifier
}
\examples{
results <- predict_LundTax2023(Lund2017)
# Include data in result
results_data <- predict_LundTax2023(Lund2017,
                               include_data = TRUE)

# Imputation
# Remove 100 genes from data
missing_genes <- sample(1:nrow(Lund2017),100)
Lund2017_missinggenes <- Lund2017[-missing_genes,]
results_imputation <- predict_LundTax2023(Lund2017_missinggenes,
                                          impute = TRUE)

}
