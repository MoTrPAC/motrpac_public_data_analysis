#' In this script we load the preprocessed datasets and reshape the data
#' so that it will be easy to do the analysis using the metafor package
#'  This is done both for the acute and longterm datasets.
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
source('/Users/David/Desktop/repos/motrpac/metaanalysis/helper_functions.R')

# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata

#' Get some basic statistics: useful to have a look at the dataset
#' 
#' @param metadata A named list. Each item represents the metadata of a cohort.
#' @return A list with some basic statistics: number of study ids (GSEs), total number
#' of samples, and additional statistics for each cohort (e.g., covered tissues, etc.)
get_basic_stats_from_metadata<-function(metadata){
  num_gses = length(unique(sapply(metadata,function(x)x$gse)))
  num_samples = length(unique(unlist(sapply(metadata,function(x)x$gsms))))
  c2n = sapply(metadata,function(x)length(x$gsms))
  c2tissue = sapply(metadata,function(x)x$tissue)
  c2tr = sapply(metadata,function(x)x$training)
  tr2samples = sapply(unique(c2tr),function(x,y,z)sum(z[y==x]),y=c2tr,z=c2n)
  tissue2samples = sapply(unique(c2tissue),function(x,y,z)sum(z[y==x]),y=c2tissue,z=c2n)
  return(list(num_gses,num_samples,c2n,c2tissue,c2tr,tr2samples,tissue2samples))
}
acute_cohorts = sapply(acute_datasets,function(x)!is.null(x$gene_fchanges))
get_basic_stats_from_metadata(acute_metadata[acute_cohorts])
longterm_cohorts = sapply(longterm_datasets,function(x)!is.null(x$gene_fchanges))
get_basic_stats_from_metadata(longterm_metadata[longterm_cohorts])

# Put the summary statistics in lists
acute_datasets_effects = lapply(acute_datasets,function(x)x$time2ttest_stats)
longterm_datasets_effects = lapply(longterm_datasets,function(x)x$time2ttest_stats)
longterm_datasets_effects = longterm_datasets_effects[sapply(longterm_datasets_effects,length)>0]
acute_datasets_effects = acute_datasets_effects[sapply(acute_datasets_effects,length)>0]

#' Put the moderator information (covariates) of a dataset in a single object
#' with: sex,age,tissue, training type, gse
#' Assumption: "other" without a description for a training type == "untrained"
#' @param metadata A named list. Each item represents the metadata of a cohort.
#' @return A named vector with the study id, tissue, training type, avg age, and proportion of males.
get_dataset_moderators<-function(metadata){
  arrs = t(sapply(metadata,function(x)c(x$gse,x$tissue,x$training,x$avg_age,x$male_prop)))
  colnames(arrs) = c("gse","tissue","training","avg_age","prop_males")
  return(arrs)
}
#' Put all data of a single gene in a single object
#' @param gene A character or an index. The gene to analyze.
#' @param dataset_effects. A list with the summary statistics.
#' @param moderators. The study covariates to add to each data frame.
#' @return A data frame with the information for the gene.
get_gene_table<-function(gene,dataset_effects,moderators){
  gene_data = lapply(dataset_effects,function(x,y)try({sapply(x,function(u,v)u[v,],v=y)}),y=gene)
  m = c()
  for(nn in names(gene_data)){
    mm = gene_data[[nn]]
    for(j in 1:ncol(mm)){
      m = rbind(m,c(nn,colnames(mm)[j],moderators[nn,],mm[1:2,j]))
    }
  }
  colnames(m)[2]="time"
  m = data.frame(m,stringsAsFactors=F)
  m = transform(m,yi=as.numeric(yi),vi=as.numeric(vi),time=as.numeric(time))
  return(m)
}

# We shall use meta-analyses for single genes at a time. 
# We therefore reshape the data to make it easier to work with
# Acute bouts
acute_genes = rownames(acute_datasets[[1]]$gene_data)
for(i in 2:length(acute_datasets)){
  if(length(acute_datasets[[i]])<3){next}
  acute_genes = intersect(acute_genes,rownames(acute_datasets[[i]]$gene_data))
}
acute_mod = get_dataset_moderators(acute_metadata)
acute_gene_tables = lapply(acute_genes,get_gene_table,dataset_effects=acute_datasets_effects,moderators=acute_mod)
names(acute_gene_tables) = acute_genes
# Longterm studies
longterm_genes = rownames(longterm_datasets[[1]]$gene_data)
for(i in 2:length(longterm_datasets)){
  if(length(longterm_datasets[[i]])<3){next}
  longterm_genes = intersect(longterm_genes,rownames(longterm_datasets[[i]]$gene_data))
}
longterm_mod = get_dataset_moderators(longterm_metadata)
longterm_gene_tables = lapply(longterm_genes,get_gene_table,dataset_effects=longterm_datasets_effects,moderators=longterm_mod)
names(longterm_gene_tables) = longterm_genes

#' Remove controls, untrained, and fat tissue from the tables.
#' We exclude these studies due to very low coverage. 
#' @param gdata A data frame with the information of a gene.
#' @return  A new filtered data frame with the information of the gene.
remove_undesired_datasets<-function(gdata){
  gdata = gdata[!grepl("control",gdata$training),]
  gdata = gdata[!grepl("other",gdata$training),]
  gdata = gdata[!grepl("yoga",gdata$training),]
  gdata = gdata[!grepl("untrained",gdata$training),]
  gdata = gdata[!grepl("fat",gdata$tissue),]
  gdata = gdata[!grepl("adipose",gdata$tissue),]
  gdata = gdata[!grepl("treatment",gdata$training),]
  return (gdata)
}
acute_gene_tables_raw = acute_gene_tables
acute_gene_tables = lapply(acute_gene_tables,remove_undesired_datasets)
longterm_gene_tables_raw = longterm_gene_tables
longterm_gene_tables = lapply(longterm_gene_tables,remove_undesired_datasets)

# Save the objects
save(acute_gene_tables_raw,acute_gene_tables,longterm_gene_tables_raw,longterm_gene_tables,
     file="PADB_dataset_level_meta_analysis_data.RData")


