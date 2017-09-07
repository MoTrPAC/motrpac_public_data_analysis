setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(metafor)
source('repos/motrpac/helper_functions.R')

# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata

# extract the data objects for the meta-analysis
acute_genes = rownames(acute_datasets[[1]]$gene_data)
for(i in 2:length(acute_datasets)){
  if(length(acute_datasets[[i]])<3){next}
  acute_genes = intersect(acute_genes,rownames(acute_datasets[[i]]$gene_data))
}
longterm_genes = rownames(longterm_datasets[[1]]$gene_data)
for(i in 2:length(longterm_datasets)){
  if(length(longterm_datasets[[i]])<3){next}
  longterm_genes = intersect(longterm_genes,rownames(longterm_datasets[[i]]$gene_data))
}

# get genes that are expressed
q_thr = 0.5
expressed_genes1 = list()
sum_counts1 = list()
for(i in 1:length(acute_datasets)){
  if (length(acute_datasets[[i]])<3){next}
  tissue = acute_metadata[[i]]$tissue
  m = acute_datasets[[i]]$gene_data[acute_genes,]
  m = apply(m,2,rank)
  m = m/nrow(m)
  curr_expressed = rowSums(m>q_thr)/ncol(m)
  if(is.null(expressed_genes1[[tissue]])){
    expressed_genes1[[tissue]] = curr_expressed
    sum_counts1[[tissue]]=1
  }
  else{
    expressed_genes1[[tissue]] = expressed_genes1[[tissue]] + curr_expressed 
    sum_counts1[[tissue]] = sum_counts1[[tissue]]+1
  }
}
expressed_genes2 = list()
sum_counts2 = list()
for(i in 1:length(longterm_datasets)){
  if (length(longterm_datasets[[i]])<3){next}
  tissue = longterm_metadata[[i]]$tissue
  m = longterm_datasets[[i]]$gene_data[longterm_genes,]
  m = apply(m,2,rank)
  m = m/nrow(m)
  curr_expressed = rowSums(m>q_thr)/ncol(m)
  if(is.null(expressed_genes2[[tissue]])){
    expressed_genes2[[tissue]] = curr_expressed
    sum_counts2[[tissue]] = 1
  }
  else{
    expressed_genes2[[tissue]] = expressed_genes2[[tissue]] + curr_expressed 
    sum_counts2[[tissue]] = sum_counts2[[tissue]]+1
  }
}
for(nn in names(expressed_genes1)){
  expressed_genes1[[nn]] = expressed_genes1[[nn]] / sum_counts1[[nn]]
}
for(nn in names(expressed_genes2)){
  expressed_genes2[[nn]] = expressed_genes2[[nn]] / sum_counts2[[nn]]
}

gs = intersect(acute_genes,longterm_genes)
plot(expressed_genes1$muscle[gs],expressed_genes2$muscle[gs])
plot(expressed_genes1$blood[gs],expressed_genes2$blood[gs])




