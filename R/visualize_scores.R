# Plot Lund Scores

# Data
# TPMSTAR665log2
# Adjusted data
load("D:/UROSCANSEQ_2024/Datasets/RNAseq/02.NF_Core314/Merged_Results_AllBatches/FinalPreprocessing/filteredgeTMM572_star_ComBatSeq_killPC6_non_centered.Rdata",verbose=T)

# Full dataset TPM
TPMSTAR <- read.csv("D:/UROSCANSEQ_2024/Datasets/RNAseq/02.NF_Core314/Merged_Results_AllBatches/filteredTPM754_star.txt", sep ="\t")


# Results
# Signatures/Immune/
load("D:/UROSCANSEQ_2024/Analysis/02.New_data/resultsStarandSalmon754filt.rda", verbose = T)
TMP_results_fullSTAR <- resultsStarandSalmon754filt[[1]]

load("D:/UROSCANSEQ_2024/Analysis/02.New_data/results665_TPM.RData")
rm(resultsStarandSalmon754filt)

# Metadata
load("D:/UROSCANSEQ_2024/Datasets/RNAseq/02.NF_Core314/Merged_Results_AllBatches/FinalPreprocessing/Metadata.Rdata")
load("D:/UROSCANSEQ_2024/Analysis/02.New_data/Metadata/Metadata665.RData")
load("D:/UROSCANSEQ_2024/Analysis/02.New_data/Metadata/Metadata536.RData")
load("D:/UROSCANSEQ_2024/Analysis/02.New_data/Metadata/metadataTMA550.RData")

TPMSTAR665log2 <- log(TPMSTAR[,Metadata665$XRNA_cohort_name]+1)

# Colors
library(viridis)

test_predict_full572 <- predict_LundTax2023_wip(data = filteredgeTMM572_star_ComBatSeq_killPC6_non_centered,
                                             include_data = TRUE,
                                             include_scores = TRUE,
                                             gene_id = "ensembl_gene_id",
                                             logTransform = FALSE)

head(test_predict_full572$scores)

test_predict_full <- predict_LundTax2023_wip(data = TPMSTAR665log2,
                                             include_data = TRUE,
                                             include_scores = TRUE,
                                             gene_id = "ensembl_gene_id",
                                             logTransform = FALSE)

head(test_predict_full$scores)

sample_order <- order(test_predict_full$scores[,"ProliferationScore"])
split <- factor(test_predict_full$predictions_7classes, levels = c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))

# Annotation subtypes
library(ComplexHeatmap)
# subtype_annotation <- HeatmapAnnotation(Lund = split,
#                                         annotation_name_side = "left",
#                                         col = list(Lund = LundTax2023Classifier::lund_colors$lund_colors))
# hm_proliferation <- Heatmap(t(test_predict_full$scores[,"ProliferationScore",drop=FALSE]),
#                             top_annotation = subtype_annotation,
#                             column_order = sample_order,
#                             name = "Proliferation Score",
#                             column_split = split,
#                             row_names_side = "left",
#                             show_column_names = FALSE)
#
# # Annotation class grade3
# grade3_annotation <- HeatmapAnnotation(GradeWHO1999 = test_predict_full$scores$GradeWHO1999,
#                                         annotation_name_side = "left",
#                                         col = list(Lund = LundTax2023Classifier::lund_colors$stage_colors))
# hm_grade3 <- Heatmap(t(test_predict_full$scores[,"GradeWHO1999_score",drop=FALSE]),
#                             column_order = sample_order,
#                      top_annotation = grade3_annotation,
#                      name = "Grade 3 Score",
#                             column_split = split,
#                             row_names_side = "left",
#                             show_column_names = FALSE)
# # Annoration class HG
# gradeHG_annotation <- HeatmapAnnotation(GradeWHO2004_2016 = test_predict_full$scores$GradeWHO2004_2016,
#                                        annotation_name_side = "left",
#                                        col = list(Lund = LundTax2023Classifier::lund_colors$stage_colors))
# grade_annotation_both <- HeatmapAnnotation(GradeWHO1999 = test_predict_full$scores$GradeWHO1999,
#                                            GradeWHO2004_2016 = test_predict_full$scores$GradeWHO2004_2016,
#                                            annotation_name_side = "left",
#                                            col = list(GradeWHO1999 = c("G1_2"="lightblue","G3"="darkblue"),
#                                                       GradeWHO2004_2016 = c("HG"="black","LG"="white")))
# # think about these colors
#
# hm_grade_hg <- Heatmap(t(test_predict_full$scores[,"GradeWHO2004_2016_score",drop=FALSE]),
#                      column_order = sample_order,
#                      top_annotation = gradeHG_annotation,
#                      name = "High Grade Score",
#                      column_split = split,
#                      row_names_side = "left",
#                      show_column_names = FALSE)
# magm_colors <- c("#000004FF","#B63679FF" ,"#FCFDBFFF")
# viridis_colors <-  c("#21908CFF", "orange")
#
# col_fun_immune <-  circlize::colorRamp2(c(-2,0,1),
#                                         c("#21908CFF","white", "#B63679FF"))
#
col_fun_proliferation <-  circlize::colorRamp2(c(quantile(test_predict_full$scores[,"ProliferationScore"], 0.05),median(test_predict_full$scores[,"ProliferationScore"]),quantile(test_predict_full$scores[,"ProliferationScore"], 0.95)),
                                               # c("#7CC4F2","#E2E7E7", "#FE3D44"),
                                               c("#21908CFF","white", "#B63679FF"))
col_fun_progression <-  circlize::colorRamp2(c(quantile(test_predict_full$scores[,"ProliferationScore"], 0.05),
                                               median(test_predict_full$scores[,"ProliferationScore"]),
                                               quantile(test_predict_full$scores[,"ProliferationScore"], 0.90)),
                                               # c("#7CC4F2","#E2E7E7", "#FE3D44"),
                                               c("#FAEBDDFF","#A11A5BFF", "#4C1D4BFF"))

col_stage <- c("Cis" = "#D87405",
               "PUNLUMP" ="pink",
               "T0" ="pink",
               "Ta" = "#72B133",
               "T1" = "#547EB6",
               "≥T2" = "#C20203",
               "Tx" = "gray"
)
annotation_proliferation_grade_prog_prostate <- HeatmapAnnotation(Lund = split,
                                                                  Stage = Metadata665[rownames(test_predict_full$scores),"Stage_simp"],
                                                             ProliferationScore = test_predict_full$scores[,"ProliferationScore"],
                                                             MolecularGradeWHO1999 = test_predict_full$scores$MolecularGradeWHO1999,
                                                             MolecularGradeWHO2016 = test_predict_full$scores$MolecularGradeWHO2016,
                                                             ProgressionScore = test_predict_full$scores$ProgressionScore,
                                                             ProgressionRisk = test_predict_full$scores$ProgressionRisk,
                                                             PossibleProstate = test_predict_full$scores$PossibleProstate,
                                                             annotation_name_side = "left",
                                                             show_legend = FALSE,
                                                             border = TRUE,
                                                             col = list(Lund = LundTax2023Classifier::lund_colors$lund_colors,
                                                                        Stage = col_stage,
                                                                        ProliferationScore = col_fun_proliferation,
                                                                        MolecularGradeWHO1999 = c("G1_2"="white","G3"="black"),
                                                                        MolecularGradeWHO2016 = c("HG"="black","LG"="white"),
                                                                        ProgressionScore = col_fun_progression,
                                                                        ProgressionRisk = c("HR" = "#A11A5BFF", "LR" = "#FAEBDDFF"),
                                                                        PossibleProstate = c("NO" = "#eae5eb", "YES" = "#ee82ee")))


# Stromal & immune infiltration ###########
# scores141up <- score141up(Data = test_predict_full$data,
#                            gene_id = "ensembl_gene_id",
#                            adjust = FALSE,
#                            logTransform = FALSE)
# scores141up_adjust <- score141up(Data = test_predict_full$data,
#                            gene_id = "ensembl_gene_id",
#                            logTransform = FALSE)
# # Add immune 141  up
# immune_scores <- cbind(Immune141_UP=scores141up_adjust$Immune141_UP,test_predict_full$scores[,6:15,drop=FALSE])
# # add stroma 141 up
# stroma_scores <- cbind(Stromal141_UP=scores141up_adjust$Stromal141_UP,test_predict_full$scores[,16:18,drop=FALSE])
immune_names <- c("Immune141_UP","B-cells", "T-cells", "T-cells CD8+", "NK-cells",
                  "Cytotoxicity Score", "Neutrophils", "Monocytic lineage", "Macrophages",
                  "M2 macrophage", "Myeloid Dendritic Cells")
stromal_names <- c("Stromal141_UP", "Endothelial cells", "Fibroblasts", "Smooth muscle")

hm_immune_scores <- Heatmap(t(scale(test_predict_full$scores[,immune_names,drop=FALSE])),
                            top_annotation = annotation_proliferation_grade_prog_prostate,
                            # col = col_fun_immune,
                            column_order = sample_order,
                            height = unit(5*ncol(test_predict_full$scores[,immune_names,drop=FALSE]),"mm"),
                            name = "Immune Scores",
                            border = TRUE,
                            column_split = split,
                            cluster_rows = FALSE,
                            show_heatmap_legend = FALSE,
                            row_names_side = "left",
                            show_column_names = FALSE)
hm_stroma_scores <- Heatmap(t(scale(test_predict_full$scores[,stromal_names,drop=FALSE])),
                            # top_annotation = annotation_proliferation_grade_prostate,
                            # col = col_fun_immune,
                            column_order = sample_order,
                            height = unit(5*ncol(test_predict_full$scores[,stromal_names,drop=FALSE]),"mm"),
                            name = "Stroma Scores",
                            border = TRUE,
                            column_split = split,
                            cluster_rows = FALSE,
                            show_heatmap_legend = FALSE,
                            row_names_side = "left",
                            show_column_names = FALSE)

all_stages <- hm_immune_scores %v% hm_stroma_scores


Ta_samples <- rownames(Metadata665[Metadata665[,"Stage_simp"] == "Ta",])
T1_samples <- rownames(Metadata665[Metadata665[,"Stage_simp"] == "T1",])
MIBC_samples <- rownames(Metadata665[Metadata665[,"Stage_simp"] == "≥T2",])

sample_orderTa <- order(test_predict_full$scores[Ta_samples,"ProliferationScore"])
splitTa <- factor(test_predict_full$predictions_7classes[Ta_samples], levels = c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))

sample_orderT1 <- order(test_predict_full$scores[T1_samples,"ProliferationScore"])
splitT1 <- factor(test_predict_full$predictions_7classes[T1_samples], levels = c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))

sample_orderMIBC <- order(test_predict_full$scores[MIBC_samples,"ProliferationScore"])
splitMIBC <- factor(test_predict_full$predictions_7classes[MIBC_samples], levels = c("UroA","UroB","UroC","GU","BaSq","Mes","ScNE"))


##### Ta plot #####
annotation_proliferation_grade_prog_prostateTa <- HeatmapAnnotation(Lund = splitTa,
                                                                  ProliferationScore = test_predict_full$scores[Ta_samples,"ProliferationScore"],
                                                                  MolecularGradeWHO1999 = test_predict_full$scores[Ta_samples,"MolecularGradeWHO1999"],
                                                                  MolecularGradeWHO2016 = test_predict_full$scores[Ta_samples,"MolecularGradeWHO2016"],
                                                                  ProgressionScore = test_predict_full$scores[Ta_samples,"ProgressionScore"],
                                                                  ProgressionRisk = test_predict_full$scores[Ta_samples,"ProgressionRisk"],
                                                                  PossibleProstate = test_predict_full$scores[Ta_samples,"PossibleProstate"],
                                                                  annotation_name_side = "left",
                                                                  show_legend = FALSE,
                                                                  border = TRUE,
                                                                  col = list(Lund = LundTax2023Classifier::lund_colors$lund_colors,
                                                                             ProliferationScore = col_fun_proliferation,
                                                                             MolecularGradeWHO1999 = c("G1_2"="white","G3"="black"),
                                                                             MolecularGradeWHO2016 = c("HG"="black","LG"="white"),
                                                                             ProgressionScore = col_fun_progression,
                                                                             ProgressionRisk = c("HR" = "#A11A5BFF", "LR" = "#FAEBDDFF"),
                                                                             PossibleProstate = c("NO" = "#eae5eb", "YES" = "#ee82ee")))

hm_immune_scoresTa <- Heatmap(t(scale(test_predict_full$scores[Ta_samples,immune_names,drop=FALSE])),
                            top_annotation = annotation_proliferation_grade_prog_prostateTa,
                            # col = col_fun_immune,
                            column_order = sample_orderTa,
                            height = unit(5*ncol(test_predict_full$scores[Ta_samples,immune_names,drop=FALSE]),"mm"),
                            name = "Immune Scores",
                            border = TRUE,
                            column_split = splitTa,
                            cluster_rows = FALSE,
                            show_heatmap_legend = FALSE,
                            row_names_side = "left",
                            show_column_names = FALSE)
hm_stroma_scoresTa <- Heatmap(t(scale(test_predict_full$scores[Ta_samples,stromal_names,drop=FALSE])),
                            # top_annotation = annotation_proliferation_grade_prostate,
                            # col = col_fun_immune,
                            column_order = sample_orderTa,
                            height = unit(5*ncol(test_predict_full$scores[Ta_samples,stromal_names,drop=FALSE]),"mm"),
                            name = "Stroma Scores",
                            border = TRUE,
                            column_split = splitTa,
                            cluster_rows = FALSE,
                            show_heatmap_legend = FALSE,
                            row_names_side = "left",
                            show_column_names = FALSE)
Ta_plot <- hm_immune_scoresTa %v% hm_stroma_scoresTa


#### T1 plot #######
annotation_proliferation_grade_prog_prostateT1 <- HeatmapAnnotation(Lund = splitT1,
                                                                    ProliferationScore = test_predict_full$scores[T1_samples,"ProliferationScore"],
                                                                    MolecularGradeWHO1999 = test_predict_full$scores[T1_samples,"MolecularGradeWHO1999"],
                                                                    MolecularGradeWHO2016 = test_predict_full$scores[T1_samples,"MolecularGradeWHO2016"],
                                                                    ProgressionScore = test_predict_full$scores[T1_samples,"ProgressionScore"],
                                                                    ProgressionRisk = test_predict_full$scores[T1_samples,"ProgressionRisk"],
                                                                    PossibleProstate = test_predict_full$scores[T1_samples,"PossibleProstate"],
                                                                    annotation_name_side = "left",
                                                                    show_legend = FALSE,
                                                                    border = TRUE,
                                                                    col = list(Lund = LundTax2023Classifier::lund_colors$lund_colors,
                                                                               ProliferationScore = col_fun_proliferation,
                                                                               MolecularGradeWHO1999 = c("G1_2"="white","G3"="black"),
                                                                               MolecularGradeWHO2016 = c("HG"="black","LG"="white"),
                                                                               ProgressionScore = col_fun_progression,
                                                                               ProgressionRisk = c("HR" = "#A11A5BFF", "LR" = "#FAEBDDFF"),
                                                                               PossibleProstate = c("NO" = "#eae5eb", "YES" = "#ee82ee")))

hm_immune_scoresT1 <- Heatmap(t(scale(test_predict_full$scores[T1_samples,immune_names,drop=FALSE])),
                              top_annotation = annotation_proliferation_grade_prog_prostateT1,
                              # col = col_fun_immune,
                              column_order = sample_orderT1,
                              height = unit(5*ncol(test_predict_full$scores[T1_samples,immune_names,drop=FALSE]),"mm"),
                              name = "Immune Scores",
                              border = TRUE,
                              column_split = splitT1,
                              cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              row_names_side = "left",
                              show_column_names = FALSE)
hm_stroma_scoresT1 <- Heatmap(t(scale(test_predict_full$scores[T1_samples,stromal_names,drop=FALSE])),
                              # top_annotation = annotation_proliferation_grade_prostate,
                              # col = col_fun_immune,
                              column_order = sample_orderT1,
                              height = unit(5*ncol(test_predict_full$scores[T1_samples,stromal_names,drop=FALSE]),"mm"),
                              name = "Stroma Scores",
                              border = TRUE,
                              column_split = splitT1,
                              cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              row_names_side = "left",
                              show_column_names = FALSE)
T1_plot <- hm_immune_scoresT1 %v% hm_stroma_scoresT1

#### >T2 plot #######
annotation_proliferation_grade_prog_prostateMIBC <- HeatmapAnnotation(Lund = splitMIBC,
                                                                    ProliferationScore = test_predict_full$scores[MIBC_samples,"ProliferationScore"],
                                                                    MolecularGradeWHO1999 = test_predict_full$scores[MIBC_samples,"MolecularGradeWHO1999"],
                                                                    MolecularGradeWHO2016 = test_predict_full$scores[MIBC_samples,"MolecularGradeWHO2016"],
                                                                    ProgressionScore = test_predict_full$scores[MIBC_samples,"ProgressionScore"],
                                                                    ProgressionRisk = test_predict_full$scores[MIBC_samples,"ProgressionRisk"],
                                                                    PossibleProstate = test_predict_full$scores[MIBC_samples,"PossibleProstate"],
                                                                    annotation_name_side = "left",
                                                                    show_legend = FALSE,
                                                                    border = TRUE,
                                                                    col = list(Lund = LundTax2023Classifier::lund_colors$lund_colors,
                                                                               ProliferationScore = col_fun_proliferation,
                                                                               MolecularGradeWHO1999 = c("G1_2"="white","G3"="black"),
                                                                               MolecularGradeWHO2016 = c("HG"="black","LG"="white"),
                                                                               ProgressionScore = col_fun_progression,
                                                                               ProgressionRisk = c("HR" = "#A11A5BFF", "LR" = "#FAEBDDFF"),
                                                                               PossibleProstate = c("NO" = "#eae5eb", "YES" = "#ee82ee")))

hm_immune_scoresMIBC <- Heatmap(t(scale(test_predict_full$scores[MIBC_samples,immune_names,drop=FALSE])),
                              top_annotation = annotation_proliferation_grade_prog_prostateMIBC,
                              # col = col_fun_immune,
                              column_order = sample_orderMIBC,
                              height = unit(5*ncol(test_predict_full$scores[MIBC_samples,immune_names,drop=FALSE]),"mm"),
                              name = "Immune Scores",
                              border = TRUE,
                              column_split = splitMIBC,
                              cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              row_names_side = "left",
                              show_column_names = FALSE)
hm_stroma_scoresMIBC <- Heatmap(t(scale(test_predict_full$scores[MIBC_samples,stromal_names,drop=FALSE])),
                              # top_annotation = annotation_proliferation_grade_prostate,
                              # col = col_fun_immune,
                              column_order = sample_orderMIBC,
                              height = unit(5*ncol(test_predict_full$scores[MIBC_samples,stromal_names,drop=FALSE]),"mm"),
                              name = "Stroma Scores",
                              border = TRUE,
                              column_split = splitMIBC,
                              cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              row_names_side = "left",
                              show_column_names = FALSE)
MIBC_plot <- hm_immune_scoresMIBC %v% hm_stroma_scoresMIBC

big_heatmap_grab <- grid.grabExpr(draw(all_stages,ht_gap = unit(1, "mm")))
Ta_heatmap_grab <- grid.grabExpr(draw(Ta_plot,ht_gap = unit(1, "mm")))
T1_heatmap_grab <- grid.grabExpr(draw(T1_plot,ht_gap = unit(1, "mm")))
MIBC_heatmap_grab <- grid.grabExpr(draw(MIBC_plot,ht_gap = unit(1, "mm")))

library(cowplot)
# pdf("Figures/test_plot_grid_half_noborder.pdf",width = 16, height =  16)
# plot_grid(big_heatmap_grab,
#           plot_grid(Ta_heatmap_grab,T1_heatmap_grab,MIBC_heatmap_grab,rel_widths = c(0.48,0.27,0.25),
#                                       nrow = 1),
#           rel_heights = c(2.2,1),
#           ncol = 1, nrow = 2)
# dev.off()

pdf("Figures/test_plot_grid_half.pdf",width = 32, height =  16)
plot_grid(big_heatmap_grab,
          plot_grid(Ta_heatmap_grab,T1_heatmap_grab,MIBC_heatmap_grab,rel_widths = c(0.48,0.27,0.25),
                    nrow = 1),
          rel_heights = c(2.2,1),
          ncol = 1, nrow = 2)
dev.off()
