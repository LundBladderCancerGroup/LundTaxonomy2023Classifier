#' LundTax2023
#'
#' This packages implements a Random Forest rule-based single-sample predictor that classifies
#' transcriptomic samples into the 5 (or 7, including subclasses) Lund Taxonomy molecular subtypes.
#' The final classifier is composed of two separate predictors applied sequentially: first a sample
#' is classified as one of the 5 main classes (Uro, GU, BaSq, Mes or ScNE), and then samples classified as
#' Uro are subclassified into UroA, UroB or UroC by a second predictor
"_PACKAGE"

#' Classifier Grade 3.
#'
#' Classifier generated with [multiclassPairs::predict_RF()], this dataset is needed for predicting 
#' the grade (1, 2, or 3) and calculating scores for incoming data.
#'
#' A Large rule_based_RandomForest object. A list of 6.
#'
#' \itemize{
#'  \item genes. Genes in ensembl format.
#'  \item rules. A set of rules for the classifier in ensembl format.
#'  \item TrainingMatrix. Binary matrix for the rules in the training data.
#'  \item boruta. Boruta results for the classifier.
#'  \item RF_classifier. Random forest classifier details.
#'  \item calls. Information on how the model was generated.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name classifier_grade3
#' @usage data(classifier_grade3)
#' @format A list of 6.
NULL

#' Classifier High Grade.
#'
#' Classifier generated with [multiclassPairs::predict_RF()], this dataset is needed for predicting 
#' the grade (high grade/lowgrade) and calculating scores for incoming data.
#'
#' A Large rule_based_RandomForest object. A list of 6.
#'
#' \itemize{
#'  \item genes. Genes in ensembl format.
#'  \item rules. A set of rules for the classifier in ensembl format.
#'  \item TrainingMatrix. Binary matrix for the rules in the training data.
#'  \item boruta. Boruta results for the classifier.
#'  \item RF_classifier. Random forest classifier details.
#'  \item calls. Information on how the model was generated.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name classifier_hg
#' @usage data(classifier_hg)
#' @format A list of 6.
NULL

#' Classifier LundTax 5c.
#'
#' Classifier as a 'rule_based_RandomForest' object. Predicts samples as one of the 5 main Lund 
#' Taxonomy molecular subtypes, Uro, GU, BaSq, Mes, or ScNE. Object includes the final RF classifier,
#' the used genes and rules in the final model, the Boruta results, and the training matrix. The 
#' training matrix is a binary matrix containing the rule values for the training data and it is
#' used for imputation purposes during the prediction if values are missing in the sample. This object
#' was generated using the [multiclassPairs::predict_RF()] function.
#'
#' A Large rule_based_RandomForest object. A list of 6.
#'
#' \itemize{
#'  \item genes. Genes in hgnc format.
#'  \item rules. A set of rules for the classifier in hgnc format.
#'  \item TrainingMatrix. Binary matrix for the rules in the training data.
#'  \item boruta. Boruta results for the classifier.
#'  \item RF_classifier. Random forest classifier details.
#'  \item calls. Information on how the model was generated.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name classifier_lundtax_5c
#' @usage data(classifier_lundtax_5c)
#' @format A list of 6.
NULL

#' Classifier LundTax 7c.
#'
#' Classifier as a rule_based_RandomForest object. Predicts samples as one of the 3 Uro subclasses, 
#' UroA, UroB, or UroC. Object includes the final RF classifier, the used genes and rules in the 
#' final model, the Boruta results, and the training matrix. The training matrix is a binary matrix
#' containing the rule values for the training data and it is used for imputation purposes during 
#' the prediction if values are missing in the sample. This object was generated using the
#' [multiclassPairs::predict_RF()] function.
#'
#' A Large rule_based_RandomForest object. A list of 6.
#'
#' \itemize{
#'  \item genes. Genes in hgnc format.
#'  \item rules. A set of rules for the classifier in hgnc format.
#'  \item TrainingMatrix. Binary matrix for the rules in the training data.
#'  \item boruta. Boruta results for the classifier.
#'  \item RF_classifier. Random forest classifier details.
#'  \item calls. Information on how the model was generated.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name classifier_lundtax_7c
#' @usage data(classifier_lundtax_7c)
#' @format A list of 6.
NULL

#' Gene List.
#'
#' Gene annotations in hgnc and ensembl format. Features, for the classifier relevant genes. 
#' For convenience, both formats are available for seamless conversion of IDs, regardless of the 
#' format of the incoming data.
#'
#' A data frame with gene information in different formats.
#'
#' \itemize{
#'  \item ensembl_gene_id. Gene annotation in ensembl format.
#'  \item hgnc_symbol. Gene annotation in HUGO format.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name gene_list
#' @usage data(gene_list)
#' @format A data frame with 1900 rows (genes) and 2 columns (hgnc_symbol and ensembl_gene_id).
NULL

#' Lund Colors.
#'
#' Standardized colors used for Lund Taxonomy subtypes and datasets.
#'
#' A list of 4 with color palettes frequently used in this package.
#'
#' \itemize{
#'  \item lund_colors. Color palette for the LundTax subtypes.
#'  \item lund_colors_transp. Same as `lund_colors` but transparent.
#'  \item stage_colors. Color palette for stages.
#'  \item dataset_colors. Color palette for the LundTax datasets.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name lund_colors
#' @usage data(lund_colors)
#' @format A list of 4.
NULL

#' Signatures
#'
#' Bundles a set of signatures in different molecular processes and the genes associated.
#' This dataset is needed for calculating the signature scores.
#'
#' A list of 6 with different signatures and the associated genes.
#'
#' \itemize{
#'  \item signatures_plot. Genes in hgnc and ensembl format associated with specific signatures.
#'  \item proliferation. Genes in hgnc and ensembl format associated with proliferation.
#'  \item progression. Genes in hgnc and ensembl format associated with progression.
#'  \item prostate. Genes in hgnc and ensembl format associated with prostate cancer.
#'  \item immune. Genes in hgnc and ensembl format associated with immuno process.
#'  \item stable_genes. Genes in hgnc and ensembl format associated with stable genes.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name signatures
#' @usage data(signatures)
#' @format A list of 6.
NULL

#' Sjodahl 2017.
#'
#' Gene expression data derived from the Sjödahl et. al. (2017) cohort. A matrix of RMA normalized
#' and ComBat adjusted gene expression values for 15697 genes (hgnc format) for 267 samples
#'
#' A data frame with gene information in different formats.
#'
#' @docType data
#' @keywords datasets
#' @name sjodahl_2017
#' @usage data(sjodahl_2017)
#' @format A data frame with 15697 rows (gene expressions) and 267 columns (samples).
NULL

#' Sjodahl 2017 Meta.
#'
#' metadata associated with the bundled expresison data.
#'
#' A data frame with meta data for the bundled samples.
#' 
#' \itemize{
#'  \item sample_id. Unique sample ID.
#'  \item cluster_order. Cluster order.
#'  \item age. Age of patient, years.
#'  \item gender. Patient gender.
#'  \item region_cx. Geographical region of CX.
#'  \item turb_grade. Pathological grading of TURB-specimen. 0 = WHO1999 Grade 2, 1 = WHO1999 Grade 3.
#'  \item turb_stage. Pathological staging of TURB-specimen. 0 = Non muscle-invasive, 1 = Muscle invasive.
#'  \item turb_skiv. Presence of keratinization or squamous differentiation in TURB-specimen. 0 = No presence, 1 = presence.
#'  \item turb_cis. Presence of cis in the TURB-specimen. 0 = No presence, 1 = presence.
#'  \item turb_lvi. Presence of lymphovascular invasion in the TURB-specimen. 0 = No presence, 1 = presence.
#'  \item turb_necros. Presence of necrosis in TURB-specimen. 0 = no, 1 = present, 2 = extensive.
#'  \item turb_apoptos. Presence of apoptosis in the TURB-specimen. 0 = No presence, 1 = presence.
#'  \item variant_histology. Variant histology present 1 = yes, 0 = no.
#'  \item cT. Clinical tumor state
#'  \item cN. Clinical node state
#'  \item hydroneph. Presence of hydronephrose. 0 = No presence, 1 = presence.
#'  \item pTCx. Pathological tumor state at cystectomy.
#'  \item pN. Patological node state at cystectomy.
#'  \item adj_chemo. Adjeuvant chemotherapy. 0 = No presence, 1 = presence.
#'  \item surv_os_event. Survival data, 1 = event, 0 = non event. 
#'  \item surv_os_time. Survival data, in months.
#'  \item surv_css_event. Cancer specific survival, 1 = event, 0 = non event. 
#'  \item surv_css_time. Cancer specific survival, in months.
#'  \item surv_pfs_event. Progression-free survival, 1 = event, 0 = non event. 
#'  \item surv_pfs_time. Progression-free survival, in months.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sjodahl_2017_meta
#' @usage data(sjodahl_2017_meta)
#' @format A data frame with 15697 rows (gene expressions) and 267 columns (samples).
NULL
