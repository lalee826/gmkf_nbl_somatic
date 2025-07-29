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

**Module 4: Calculating tumor mutation burden and plotting**

**Module 5: Mutation signature analysis, clustering, plotting, and survival analysis**

**Module 6: RNA pathway signature identification and validation in external data**

**Module 7: Gene set enrichment analysis of KEGG pathways against an external data set**
