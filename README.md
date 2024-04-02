# MoTrPAC: Analysis of publicly available exercise datasets

## Human data

This repository can be used to reproduce the preprocessing and (meta-)analysis/regression of our publicly available human transcriptomics exercise studies. The manually curated metadata of all samples can be found [here](https://docs.google.com/spreadsheets/d/1QI3cPm96OYV2R9XKhpp2jC2zLQHbMojx/edit?usp=sharing&ouid=104737748241406013814&rtpof=true&sd=true).

The data here is from more than 40 studies, but the data are heterogeneous statistically and clinically. This creates several challeneges for a robust meta-analysis. We address these issues by (1) fitting a model for each gene while accounting for the explained variability by moderators including: sex, age, training type, and time post exercise (when available), (2) applying a set of filters on the results to promote replicability, as suggested in Sweeney et al. 2017 NAR, and (3) integrating the meta-analysis results with various biological networks using systems biology methods.

Below we describe the processed data that can be used for analysis. Assume we are interested in analysis of acute bout data from muscle tissues. The scripts in data_collection can be used to generate the database. The flow takes as input the manually curated sample information and creates the gene expression profiles. 

The output of the preprocessing flow is an RData file per subsection of the database. These are the current available files:
1. [human_ge_cohort_preprocessed_db_acute.RData](https://drive.google.com/file/d/12F9sBgllxpA3UhWdQ9ND-BrkIE4AeveA/view?usp=sharing) - acute bout cohorts
2. [human_ge_cohort_preprocessed_db_longterm.RData](https://drive.google.com/file/d/1gaKz4GQ9FlZJoRzUJ7Ptqob4BJ82WMwH/view?usp=sharing) - long-term training program cohorts

Loading each of these files will add the following objects to the R session(comments in parenthesis are about reproducing the data using the scripts):

1. cohort_data: a list of cohort datasets. Each entry in the list has a name and a group of data matrices:
  1. probe_data - probe level of microarray experiments. Can be transcript data for RNA-Seq datasets.
  2. gene_data - a gene level data that results from averaging rows in the probe_data matrix
  3. (Optional if preprocessing was successful) gene_fold_changes - for each gene and subject this matrix contains the fold change of a time-point vs. the base line.
  4. (Optional if preprocessing was successful) probe_fchanges - same as above, but computed using a slower algorithm. Was currently kept for QA reasons.
  5. (Optional if preprocessing was successful) time2ttest_stats - this matrix gives the **summary statistics** for each gene and each time point that is not the baseline: mean and variance of fold change, and the p-value of a paired t-test. 

2. cohort_metadata: a list of cohort metadata. Cohort names fit those in cohort_data. Each cohort has a set of sample identifiers (usually GEO GSM ids), and additional information such as the tissue, training type, study (GEO GSE id), etc.

3. sample_metadata: a list of sample features. Each feature has a vector whose names are the sample ids. Currently available features: sex, age (partial information), time (e.g., 1h after acute bout), subject id, training type, dataset id.
  
## Meta-analysis/regression

Using the data in the objects above it is very easy to perform a meta-analysis of the calculated effects across the cohorts.

For example, using the metafor package we can perform a mixed-effects meta-analysis for each gene:

g_fc ~ b0 + b1*x_tr + 1|dataset 

where b0 is the intersect, b1 is the training effect (e.g., binary endurance or not), and 1|dataset is a random effect.

In theory we may want to include all possible moderators (tissue, time, etc) and perform a single analysis per gene. However, such analysis will be confounded in our case as the exercise cohorts are very heterogeneous and there will not be a reasonable coverage of the moderators. In the paper we describe how to mitigate these issues.

There are a few scripts that are used to meta-analyze the data. *metaanalysis/simplified_mod_data_prep_and_naive_RE.R* takes the input data, simplifies the moderators, creates a table per gene and runs a simple random effects meta-analysis without using the covariates/moderators. *metaanalysis/simplified_moderators_metaanalysis.R* collects the mete-analysis results,  performs gene selection and interpretation, and generates the different display items that summarize the analysis (figures in pdf and tables).

The actual meta-analysis is performed using auxiliary scripts called *wrapper_for_run_meta_analysis_on_sherlock.R* and *run_meta_analysis_on_sherlock.R*. These were created for running the meta-analysis in parallel using multiple CPUs. Nevertheless, the code for running all possible models per gene are implemented in these scripts. Also, these scripts take as input an RData file called [meta_analysis_input.RData](https://drive.google.com/file/d/1OWt3kl0PAyx2F68qtzYytxqHc0LTqjxL/view?usp=sharing). Their output can be downloaded [here](https://drive.google.com/file/d/1rr2PXBp_DSSupmPWzmrIrd98kaX1gnHb/view?usp=sharing). It is used within the metaanalysis/simplified_moderators_metaanalysis.R script for gene selection and for creating figures.

## Contact

For questions and suggestions please contact us at davidama@stanford.edu

