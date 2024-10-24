\name{plot_signatures}
\alias{plot_signatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot heatmap of classification results
}
\description{
This function plots a heatmap including genes and signatures of interest, with prediction results and scores on top.
}
\usage{
  plot_signatures(results_object,
           data = NULL,
           title = "",
           gene_id = c("hgnc_symbol","ensembl_gene_id","entrezgene")[1],
           annotation = c("5 classes", "7 classes")[2],
           plot_scores = TRUE,
           show_ann_legend = FALSE,
           show_hm_legend = FALSE,
           set_order = NULL,
           ann_height = 6,
           font.size = 8,
           norm = TRUE
  )
}
\arguments{
  \item{results_object}{Results object generated by the predict_LundTax2023 function}

  \item{data}{Matrix, data frame or multiclassPairs_object of gene expression values with genes in rows and samples in columns. One single sample can be classified, but it should also be in matrix format: one column with gene identifiers as rownames. The default gene identifier is HGNC symbols, but they can also be provided in ensembl gene or entrezgene IDs}

  \item{title}{Heatmap title}

  \item{gene_id}{Specify the type of gene identifier used in the data. Default value is hgnc_symbol}
  \itemize{
    \item "hgnc_symbol" for HUGO gene symbols
    \item "ensembl_gene_id" for Ensembl gene IDs
    \item "entrezgene" for Entrez IDs
  }

  \item{annotation}{Use 5 ("5 classes") or 7 class ("7 classes") annotation. Default is "7 classes".}

  \item{plot_scores}{Logical, plot prediction scores for each class, default = TRUE}

  \item{show_ann_legend}{Logical, show annotation legend (Lund classes), default = FALSE}

  \item{show_hm_legend}{Logical, show heatmap legend, fefault = FALSE}

  \item{set_order}{Set sample order. Samples are split by subtype, and ordered within each subtype. By default, samples are order by late/early cell cycle signatures ratio (late/early ratio, low to high)}

  \item{font.size}{Font size, default = 8}

  \item{norm}{Logical, normalize the data into Z-scaled values, default = TRUE}

  \item{ann_heigh}{Annotation height in cm, default = 6}
}

\value{
Draws heatmap and silently returns the sample order.
}
\references{
This function is based on \href{https://github.com/jokergoo/ComplexHeatmap}{complexHeatmap}.

Zuguang Gu, et al., Complex heatmaps reveal patterns and correlations in multidimensional genomic data, Bioinformatics, 2016.
}

\examples{
  # Including data in results object
  results <- predict_LundTax2023(Lund2017,
                                 include_data = TRUE)
  plot_signatures(results)

  # Loading data separately
  results <- predict_LundTax2023(Lund2017)
  plot_signatures(results, data = Lund2017)

  # 5 class annotation and Without prediction scores
  plot_signatures(results,
                  data = Lund2017,
                  annotation = "5 classes",
                  ann_height = 0.5,
                  plot_scores = FALSE)

  # Plot and get sample order
  sample_order <- plot_signatures(results,
                                  data = Lund2017)
  print(head(sample_order))

  # Save to pdf
  pdf("heatmap_example.pdf", width = 15, height = 10)
  plot_signatures(results,
                  data = Lund2017,
                  ann_height = 0.5,
                  plot_scores = FALSE)
  dev.off()
}

