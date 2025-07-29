## Data

The files here are results files or supplementary data files obtained from analyses performed in the manuscript that can be used in the scripts to produce the visualizations or perform further analyses.

- gene_signatures: collection of manually curated gene sets used in ssGSEA analysis
- DeepSVR_training_data_preprocessed.pkl: pickle file of DeepSVR training data
- RNAseq_rsem_coding_TPM.tsv: Bulk RNA sequencing count data in the GMKF cohort
- TARGET_harmonized_2018-03-31.csv: clinical data of the TARGET neuroblastoma cohort
- adrn_mes_genes.csv: Canonical gene lists used in determing mes/adrn scores in neuroblastoma
- allSNVList.RDS: An R data file of SNVs made ready for signature decomposition from the GMKF cohort
- classifier_results_vars_predicted.tsv: Classifier predictions used in the selecting the final consensus variant set in the manuscript
- coding_mutations_annotated.tsv: All variants called by two or more callers with annotations of pathogenecity made in the GMKF cohort
- cognbl_gmkf_clindata.csv: Clinical data of the GMKF cohort
- cosmicv3_sig_contribution.tsv: Mutational signature decomposition results of the GMKF cohort
- excluded_samples.RDS: R data file of samples that are excluded from analysis in the GMKF cohor due to suspicion of normal or cross-contamination
- gmkf_gsva_for_classification.tsv: ssGSEA results of bulk RNA sequencing data in the GMKF cohort. Results shown here are scores binned into one of 40 quantiles from -20 to 20.
- gmkf_matched_features.pkl: bam-readcount output for GMKF tumors input into the variant refinment classifier in order to be assigned predicted labels
- gmkf_tumor_allPathEvents.tsv: A comprehensive table of all pathogenic events detected in the GMKF cohort
- gsva_results_quantiles_plot_allgenes.tsv: Plotting data for ssGSEA results in the GMKF cohort
- gsva_results_quantiles_plot_allgenes_target.tsv: Plotting data for ssGSEA results in the TARGET cohort
- mes_adrn_scores.tsv: Analysis results of adrn/mes status of GMKF neuroblastoma tumors
- mutsig_cluster_results.tsv: Results table of clustering GMKF tumors by mutational signature profiles
- sample_manifest.csv: Metadata file for GMKF tumors hosted on the cloud platform Cavatica
- target_expmat.tsv: Bulk RNA expression count data for the TARGET cohort
- target_gsva_for_classification.tsv: ssGSEA results of bulk RNA sequencing data in the TARGET cohort. Results shown here are scores binned into one of 40 quantiles from -20 to 20.
- tumor_clinical_data.tsv: Clinical and molecular features table of the GMKF cohort
- tumor_purity_results.tsv: Results of tumor purity analysis performed in bulk RNA sequencind data of the GMKF cohort
- unique_coding_mutations_full.tsv: Unannotated SNVs and short insertion/deletion mutations called by only one variant caller
- variants_bamrc_analyzed.tsv: Bam-readcount output for all variants called with the variant caller ensemble
