# MoTrPAC: Analysis of publicly available exercise genomics datasets

## Human data

This repository can be used to reproduce the preprocessing and (meta-)analysis of publicly available genomic exercise studies. Currently, only human gene expression data were analyzed. Manually curated metadata of all samples can be found [here](https://storage.googleapis.com/motrpac-portal-user-davidama/GEO_sample_metadata.xlsx).

Assume we are interested in analysis of acute bout data from muscle tissues. The flow of the analysis is as follows:
(TODO: specify the scripts and their order)


The output of the preprocessing flow is an RData file per subsection of the database. These are the current available files:
1. [human_ge_cohort_preprocessed_db_acute.RData](https://storage.googleapis.com/motrpac-portal-user-davidama/human_ge_cohort_preprocessed_db_acute.RData) - acute bout cohorts
2. [human_ge_cohort_preprocessed_db_longterm.RData](https://storage.googleapis.com/motrpac-portal-user-davidama/human_ge_cohort_preprocessed_db_longterm.RData) - long-term training program cohorts

Loading each of these files will add the following objects to the R session(comments in parenthesis are about reproducing the data using the scripts):

1. cohort_data: a list of cohort datasets. Each entry in the list has a name and a group of data matrices:
  1. probe_data - probe level of microarray experiments. Can be transcript data for RNA-Seq datasets.
  2. gene_data - a gene level data that results from averaging rows in the probe_data matrix
  3. (Optional if preprocessing was successful) gene_fold_changes - for each gene and subject this matrix contains the fold change of a time-point vs. the base line.
  4. (Optional if preprocessing was successful) probe_fchanges - same as above, but computed using a slower algorithm. Was currently kept for QA reasons.
  5. (Optional if preprocessing was successful) time2ttest_stats - this matrix gives the **summary statistics** for each gene and each time point that is not the baseline: mean and variance of fold change, and the p-value of a paired t-test. 

2. cohort_metadata: a list of cohort metadata. Cohort names fit those in cohort_data. Each cohort has a set of sample identifiers (usually GEO GSM ids), and additional information such as the tissue, training type, study (GEO GSE id), etc.

3. sample_metadata: a list of sample features. Each feature has a vector whose names are the sample ids. Currently available features: sex, age (partial information), time (e.g., 1h after acute bout), subject id, training type, dataset id.
  
## Meta-analysis

Using the data in the objects above it is very easy to perform either a pooled analysis at the sample level or a meta-analysis of observed effects across the cohorts.

For example, using the metafor package we can perform a mixed-effects meta-analysis for each gene:

g_fc ~ b0 + b1*x_tr + 1|dataset 

where b0 is the intersect, b1 is the training effect (e.g., binary endurance or not), and 1|dataset is a random effect.

In theory we may want to include all possible moderators (tissue, time, etc) and perform a single analysis per gene. However, such analysis will be confounded in our case as the exercise cohorts are very heterogeneous and there will not be a reasonable coverage of the moderators. Moreover, most datasets did not include untrained controls.

We therefore chose to analyze selected time windows in each tissue: acute,muscle: 2-5h, acute,blood: 0-2h, and longterm (muscle and blood): 70-180 days of training.

The results are given [here](https://storage.googleapis.com/motrpac-portal-user-davidama/time_window_meta_analysis_results.RData). Loading this RData file will result in the following objects added to the R workspace:
1. window_based_metadnalysis_raw_results: raw results of the metanalaysis for all genes
2. selected_genes: selected genes for each analysis above (acute or longterm x muscle or blood)
3. topgo_enrichments: enriched GO terms (0.01 FDR) for the four gene sets above

Genes were selected based on the following criteria to minimize the inflation of false positives observed in past meta-analyses (Sweeny et al. NAR 2016). Note that each gene can have several measured effects (b0, b1 above)
1. FDR < 0.01 for at least one of the measured fold change effects
2. Effect sizes should be at least 0.5 for acute data analysis or 0.25 for longterm data analysis
3. Effect sizes should be different than the average fold change observed in all available controls (e.g., muscle biopsies from untrained subjects). Thresholds are the same as in step 2 above. 

Code for reproducing the analysis will be added in the future.

## Contact

For questions and suggestions please contact us at davidama@stanford.edu
