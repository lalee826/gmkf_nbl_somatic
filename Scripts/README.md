## Scripts

**Module 1: Somatic variant refinement classifier**

- 01-01_crossValidation.py: Using training data, performs 5-fold cross validation and searches for best model parameters for logistic regression with elastic net regularization, SVM with a radial basis function kernel, random forest, and gradient boosted trees
- 01-02_trainClassifierAndPredict.py: Trains a neural net on DeepSVR data, finds the best parameters in 5-fold cross-validation, and makes predictions on our dataset, including both consensus variants and unique variants.
- 01-03_trainClassifierAndPredictXGBoost.py: Trains a gradient boosted tree model using the xgboost package on DeepSVR data, finds the best parameters in 5-fold cross-validation, and makes predictions on our dataset including both consensus and unique variants.

**Module 2: Creating consensus variant set and plotting**

- 02-01_integrate_variants_and_plot.R: Creates a consensus variant set from both the consensus variants and unique variants using the logic detailed in the manuscript. Plots overlap of the variant caller ensemble after integrating the variant refinement classifier, as well as the proportion of calls made by different numbers of callers in the final consensus.

**Module 3: Parsing structural variants, overlapping with genes, and creating oncoplot**

- 03_01_prepareStructuralVariantVCF.py: Transforms structural variant VCF files from GRIDSS, Lumpy, Manta, and Novobreak into events tables ready for downstream processing by clearing headers and splitting out call metrics into individual columns. Events files can be filtered further based on user-defined thresholds for any given metric.
- 03-02_makeSVEventsTableGridss.py: Parses GRIDSS structural variant events tables by joining similar breakpoints into events with multiple ends and tabulating significant metrics for further filtering based on user-defined thresholds. 
- 03-03_makeSVEventsTableLumpy.py: Parses Lumpy structural variant events tables by joining similar breakpoints into events with multiple ends and tabulating significant metrics for further filtering based on user-defined thresholds.
- 03-04_makeSVEventsTableManta.py: Parses Manta structural variant events tables by joining similar breakpoints into events with multiple ends and tabulating significant metrics for further filtering based on user-defined thresholds.
- 03-05_makeSVEventsTableNovobreak.py: Parses Novobreak structural variant events tables by joining similar breakpoints into events with multiple ends and tabulating significant metrics for further filtering based on user-defined thresholds.
- 03-06_IntersectFilteredSVsWithGenes.R: Finds structural variant mutations that potentially disrupt coding genes. This script uses helper functions found in the script sv_overlap_functions.R that perform fast intersection of structural variant breakpoints with bed files.
- 03-07_MakeOncoplotAllVariantTypes.R: This script generates a custom oncoplot with all pathogenic events, including SNVs, small insertion/deletions, structural variants, and copy number events. The logic for calling pathogenic mutations can be found in the manuscript.

**Module 4: Calculating tumor mutation burden and plotting**

- 04-01_call_tmb.R: Calculated tumor mutation burden in coding regions using bed files from hg38 canonical exons to determine coding bases.
- 04-02_calculate_mesadrn_scores.R: Calculates tumor mesenchymal and adrenergic scores from bulk RNA sequencing data
- 04-03_tmb_plot.R: Visualizes tumor mutation burden across risk status and MYCN amplification status along with tumor purity, mes/adrn classification, and tumor mutational signature cluster.

**Module 5: Mutation signature analysis, clustering, plotting, and survival analysis**

- 05-01_signature_decomposition.R: Demonstrates preparation of example MAF files for input into the deconstructSigs package for mutational signature decomposition.
- 05-02_cluster_signatures.py: This script imports mutational signature decomposition results and clusters samples by their signature profiles using either a Gaussian mixture model or Kmeans. Ideal cluster number is visualized and determined by BIC analysis.
- 05-03_survival_analysis_signatures.R: Performs and visualizes survival analysis of samples based on mutational signature clustering results.

**Module 6: RNA pathway signature identification and validation in external data**

- 06-01_gsvaAnalysis.R: Performs single sample gene set enrichment analysis using manually curated gene sets related to cancer using the GSVA package.
- 06-02_clusterGSVAandPlot.R: Performs clustering and visualization through custom heatmaps of RNA signature profiles obtained through ssGSEA analysis.
- 06-03_gsvaAnalysisTarget.R: Performs single sample gene set enrichment analysis using manually curated gene sets related to cancer using the GSVA package. The external dataset being analyzed in this script is the neuroblastoma TARGET cohort.
- 06-04_trainSVMClassifyTarget.R: External validation of gene signature clustering results is performed on the TARGET cohort by training an SVM-RBF model on the GMKF dataset. The classifier is used to predict TARGET samples that had identical analysis performed. Group assignment is made by running 50 iterations and taking the mode group. Visualization is performed through tSNE and UMAP.
- 06-05_gsvaSurvival.R: Performs and visualizaes survival analysis of RNA pathway signature groups

**Module 7: Gene set enrichment analysis of KEGG pathways against an external data set**

- 07-01_GOEnrichmentAnalysis.R: Performs and visualizes gene ontology term enrichment in a ranked list of pathogenically mutated genes in the GMKF dataset.
- 07-02_externalGSEAKEGGTerms.R: Visualizes the results of GSEA of a ranked list of pathogenically mutated genes in the GMKF dataset against a ranked list of pathogenically mutated genes in an external cohort described in the manuscript.

**Module 8: Multivariate Cox proportional hazard modeling**

-08-01_CPHRegressionHRM.R: Performs covariate regression in CPH models after aggregating recurrently altered molecular features. Visualizes the regression and plots hazard ratios and univariate p-values and significance of model fit for each feature in high-risk, MYCN-amplified tumors.
-08-02_CPHRegressionHRN.R: Performs covariate regression in CPH models after aggregating recurrently altered molecular features. Visualizes the regression and plots hazard ratios and univariate p-values and significance of model fit for each feature in high-risk, non MYCN-amplified tumors. 
