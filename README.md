
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LundTax2023Classifier

<!-- badges: start -->
<!-- badges: end -->

This package implements a Random Forest rule-based single-sample
predictor that classifies gene expression samples into the 5 (or 7,
including subclasses) Lund Taxonomy molecular subtypes. The final
classifier is composed of two separate predictors applied sequentially -
first a sample is classified as one of the 5 main classes (Uro, GU,
BaSq, Mes or ScNE), and then samples classified as Uro are subclassified
into UroA, UroB or UroC by a second predictor. The package includes a
sample dataset (Lund2017) to run the classifier.

## Installation

You can install LundTax2023Classifier from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("LundBladderCancerGroup/LundTaxonomy2023Classifier")
```

## Usage

### Prediction

``` r
predict_LundTax2023(data, include_data = FALSE, include_scores = TRUE, ...)
```

Where `data` is a matrix, data frame or multiclassPairs_object of gene
expression values with genes as HUGO gene symbols in rows and samples in
columns.

`include_data` is a logical value indicating if the gene expression
values should be included in the results object.

`include_scores` is a logical value indicating if the prediction scores
should be included in the results object.

The predict function includes an imputation feature to handle missing
genes in the data. This can be accessed by adding the “impute = TRUE”
argument.

#### Example

``` r
library(LundTax2023Classifier)
results <- predict_LundTax2023(Lund2017)
str(results)
#> List of 3
#>  $ scores              : num [1:301, 1:8] 0.9924 0.0016 0.0682 0.993 0.9808 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...
#>   .. ..$ : chr [1:8] "Uro" "UroA" "UroB" "UroC" ...
#>  $ predictions_7classes: Named chr [1:301] "UroA" "Mes" "GU" "UroA" ...
#>   ..- attr(*, "names")= chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...
#>  $ predictions_5classes: Named chr [1:301] "Uro" "Mes" "GU" "Uro" ...
#>   ..- attr(*, "names")= chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...

# Include data in result
results_data <- predict_LundTax2023(Lund2017,
                               include_data = TRUE)
str(results_data)
#> List of 4
#>  $ data                : num [1:15697, 1:301] 4.25 8.35 6.69 7.26 3.74 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:15697] "A1CF" "A2M" "A2ML1" "A4GALT" ...
#>   .. ..$ : chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...
#>  $ scores              : num [1:301, 1:8] 0.9924 0.0016 0.0682 0.993 0.9808 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...
#>   .. ..$ : chr [1:8] "Uro" "UroA" "UroB" "UroC" ...
#>  $ predictions_7classes: Named chr [1:301] "UroA" "Mes" "GU" "UroA" ...
#>   ..- attr(*, "names")= chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...
#>  $ predictions_5classes: Named chr [1:301] "Uro" "Mes" "GU" "Uro" ...
#>   ..- attr(*, "names")= chr [1:301] "p1404_1.CEL" "2.CEL" "3.CEL" "4.CEL" ...
```

``` r
# Imputation
# Remove 100 genes from data
missing_genes <- sample(1:nrow(Lund2017),100)
Lund2017_missinggenes <- Lund2017[-missing_genes,]
results_imputation <- predict_LundTax2023(Lund2017_missinggenes,
                                          impute = TRUE)
#> These genes are not found in the data:
#> ZNF561 LRRC8A FANCC PRLR PPIE
#> Gene names should as rownames and sample names as columns!
#> Check the genes in classifier object to see all the needed genes.
#> Check if '-' or ',' symbols in the gene names in your data. You may need to change it to '_' or '.'
#> Missed genes will be imputed to the closest class for each sample!
#> These genes have NAs:
#> ZNF561 LRRC8A FANCC PRLR PPIE
#> These genes will be imputed to the closest class for each sample with NAs
```

The classifier returns a list of up to 4 elements: - `data` original
gene expression values. - `scores` matrix containing predictions scores
for 8 classes (Uro, UroA, UroB, UroC, GU, BaSq, Mes and ScNE). -
`predictions_7classes` named vector, with sample names as names and
subtype labels as values. - `predictions_5classes` named vector, with
sample names as names and subtype labels as values.

Both data and scores can be excluded or included from the final output
in the include_data and include_scores parameters, respectively.

### Plotting

A plotting function is included to draw a heatmap showing genes, gene
signatures and scores of interest. This function requires the
[ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) package.

``` r
plot_signatures(
  results_object,
  data = NULL,
  title = "",
  annotation = c("5 classes", "7 classes")[2],
  plot_scores = TRUE,
  show_ann_legend = FALSE,
  show_hm_legend = FALSE,
  set_order = NULL,
  ann_height = 6,
  font.size = 8,
  norm = c("scale", NULL)[1]
)
```

Parameters:

- `results_object` is a list resulting from applying the
  predict_LundTax2023 function

- `data` is is a matrix, data frame or multiclassPairs_object of gene
  expression values with genes as HUGO gene symbols in rows and samples
  in columns. This can be included if the results_object does not
  include the data, and samples should be in the same order as in the
  results object.

- `title` title for the heatmap.

- `annotation` is acharacter indicating if 5 (“5 classes”) or 7 class
  (“7 classes”) annotations should be plotted.

- `plot_scores` is a logical value indicating if the prediction scores
  should be plotted.

\-`show_ann_legend` is a logical value indicating if the annotation
legend should be shown.

- `show_hm_legend` is a logical value indicating if the heatmap legend
  should be shown.

- `set_order` is a logical value indicating if the prediction scores
  should be plotted.

\-`font.size` font size, default is 8.

- `norm` indicates if data should be scaled/Z-normalized. If “NULL”,
  data is plotted as is.

- `ann_heigh` annotation height in cm, default is 6.

#### Example

``` r
# Including data in results object
results <- predict_LundTax2023(Lund2017,
                              include_data = TRUE)
plot_signatures(results)
```

<img src="man/figures/README-heatmap-1.png" width="100%" />

# Session info

``` r
sessionInfo()
#> R version 4.2.1 (2022-06-23 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19045)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] ComplexHeatmap_2.14.0         LundTax2023Classifier_0.0.0.9
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.11         highr_0.10          compiler_4.2.1     
#>  [4] RColorBrewer_1.1-3  iterators_1.0.14    tools_4.2.1        
#>  [7] digest_0.6.33       rdist_0.0.5         clue_0.3-65        
#> [10] evaluate_0.23       lattice_0.22-5      png_0.1-8          
#> [13] rlang_1.1.2         Matrix_1.6-3        foreach_1.5.2      
#> [16] cli_3.6.1           rstudioapi_0.15.0   yaml_2.3.7         
#> [19] parallel_4.2.1      xfun_0.41           fastmap_1.1.1      
#> [22] cluster_2.1.4       ranger_0.16.0       knitr_1.45         
#> [25] GlobalOptions_0.1.2 S4Vectors_0.36.2    IRanges_2.32.0     
#> [28] stats4_4.2.1        GetoptLong_1.0.5    rmarkdown_2.25     
#> [31] codetools_0.2-19    matrixStats_1.1.0   htmltools_0.5.7    
#> [34] BiocGenerics_0.44.0 shape_1.4.6         circlize_0.4.15    
#> [37] colorspace_2.1-0    doParallel_1.0.17   crayon_1.5.2       
#> [40] rjson_0.2.21
```
