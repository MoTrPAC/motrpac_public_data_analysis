# MoTrPAC: Analysis of publicly available exercise datasets

## Human data

This repository can be used to reproduce the preprocessing and (meta-)analysis of publicly available gene expression exercise studies. Currently, only human gene expression data were analyzed. Manually curated metadata of all samples can be found [here](https://storage.googleapis.com/motrpac-portal-user-davidama/GEO_sample_metadata.xlsx).

Assume we are interested in analysis of acute bout data from muscle tissues. The scripts in data_collection can be used to generate the database. The flow takes as input the manually curated sample information and creates the gene expression profiles. 

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

Using the data in the objects above it is very easy to perform a meta-analysis of the calculated effects across the cohorts.

For example, using the metafor package we can perform a mixed-effects meta-analysis for each gene:

g_fc ~ b0 + b1*x_tr + 1|dataset 

where b0 is the intersect, b1 is the training effect (e.g., binary endurance or not), and 1|dataset is a random effect.

In theory we may want to include all possible moderators (tissue, time, etc) and perform a single analysis per gene. However, such analysis will be confounded in our case as the exercise cohorts are very heterogeneous and there will not be a reasonable coverage of the moderators. Moreover, most datasets did not include untrained controls.

We therefore simplified the moderators before the analysis. For example, we defined time windows in each tissue: acute,muscle: 2-5h, acute,blood: 0-2h, and longterm (muscle and blood): 70-180 days of training. 

The main analysis script that takes the database and generates all results and diplay items is metaanalysis/simplified_moderators_metaanalysis.R. It takes the curated database above, simplifies the moderators, and defines the input datasets for the analysis (a table for each gene). It then collects the mete-analysis results and perform the gene selection and interpretation analyses (details below).

The actual meta-analysis is performed using auxiliary scripts called wrapper_for_run_meta_analysis_on_sherlock.R and run_meta_analysis_on_sherlock.R. These were created for running the meta-analysis in parallel using multiple CPUs. Nevertheless, the code for running all possible models per gene are implemented in these scripts. Also, these scripts take as input an RData file called [meta_analysis_input.RData](https://storage.googleapis.com/motrpac-portal-user-davidama/meta_analysis_input.RData). Their output can be downloaded [here](https://storage.googleapis.com/motrpac-portal-user-davidama/meta_analysis_results.RData). It is used within the metaanalysis/simplified_moderators_metaanalysis.R script for gene selection and for creating figures.

Genes were selected based several criteria to minimize the inflation of false positives observed in past meta-analyses (some are based on Sweeny et al. NAR 2016). 
1. The selected model of a gene either include moderators or not. Models that include moderators must have an improvement in their AICc score of at least five as compared to the model without any moderators. If no such model exists, then the model without moderators is selected only if its I^2 score is lower than 50%
2. There is at least one effect size > 0.5 for acute data analysis or > 0.25 for longterm data analysis
3. FDR < 0.01 for the selected model

## Contact

For questions and suggestions please contact us at davidama@stanford.edu
