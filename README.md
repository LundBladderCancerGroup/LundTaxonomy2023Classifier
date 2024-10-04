# LundTax2023Classifier
A Random Forest rule-based single-sample predictor that classifies urothelial cancer samples to one 
of the the 5 (or 7, including subclasses) Lund Taxonomy molecular subtypes, based on gene expression
data. 

This classifier is composed of two separate predictors applied sequentially - first a sample is 
classified as one of the 5 main classes (Uro, GU, BaSq, Mes or ScNE). Then, samples classified as
Uro are further sub classified into UroA, UroB or UroC by a second predictor. 

Furthermore, the package can also be deployed to visualize the prediction calls and signatures used 
for subtype prediction. For more information on how to use the package, please refer to package 
documentation and vignettes.

## Installation
To install the package, open an R session and run the following command (requires [devtools](https://github.com/r-lib/devtools)):
``` r
devtools::install_github("LundBladderCancerGroup/LundTaxonomy2023Classifier")
```

## Usage
Please refer to the package [documentation]() and [vignettes]()
