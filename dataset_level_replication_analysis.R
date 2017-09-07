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

acute_datasets_effects = lapply(acute_datasets,function(x)x$time2ttest_stats)
longterm_datasets_effects = lapply(longterm_datasets,function(x)x$time2ttest_stats)
longterm_datasets_effects = longterm_datasets_effects[sapply(longterm_datasets_effects,length)>0]
acute_datasets_effects = acute_datasets_effects[sapply(acute_datasets_effects,length)>0]

# in this script meta-analyses are done by gene
# we therefore need to create a table for each gene
# get dataset info - moderators
# tissue, training type, gse
# Assumption, "other" without a description for a training type == "untrained"
get_dataset_moderators<-function(metadata){
  arrs = t(sapply(metadata,function(x)c(x$gse,x$tissue,x$training)))
  tr_desc = paste(arrs[,3],sapply(metadata,function(x)x$training_desc),sep=';')
  arrs[grepl(tr_desc,pattern = 'control',ignore.case = T) | grepl(tr_desc,pattern="untrain",ignore.case = T),3] = "control"
  arrs[grepl(tr_desc,pattern = 'other;NA',ignore.case = T) | grepl(tr_desc,pattern="untrain",ignore.case = T),3] = "control"
  colnames(arrs) = c("gse","tissue","training")
  return(arrs)
}

# extract the data objects for the meta-analysis
acute_genes = rownames(acute_datasets[[1]]$gene_data)
for(i in 2:length(acute_datasets)){
  if(length(acute_datasets[[i]])<3){next}
  acute_genes = intersect(acute_genes,rownames(acute_datasets[[i]]$gene_data))
}
acute_mod = get_dataset_moderators(acute_metadata)
longterm_genes = rownames(longterm_datasets[[1]]$gene_data)
for(i in 2:length(longterm_datasets)){
  if(length(longterm_datasets[[i]])<3){next}
  longterm_genes = intersect(longterm_genes,rownames(longterm_datasets[[i]]$gene_data))
}
longterm_mod = get_dataset_moderators(longterm_metadata)

# transform the datasets into matrices that can be plotted
acute_pvals_matrix = c()
for(d in names(acute_datasets_effects)){
  l = acute_datasets_effects[[d]]
  times = names(l)
  for (tt in times){
    currname = paste(c(d,tt),collapse=';')
    acute_pvals_matrix = cbind(acute_pvals_matrix,l[[tt]][acute_genes,"p"])
    colnames(acute_pvals_matrix)[ncol(acute_pvals_matrix)] = currname
  }
}
longterm_pvals_matrix = c()
for(d in names(longterm_datasets_effects)){
  l = longterm_datasets_effects[[d]]
  times = names(l)
  for (tt in times){
    currname = paste(c(d,tt),collapse=';')
    longterm_pvals_matrix = cbind(longterm_pvals_matrix,l[[tt]][longterm_genes,"p"])
    colnames(longterm_pvals_matrix)[ncol(longterm_pvals_matrix)] = currname
  }
}
save(longterm_pvals_matrix,acute_pvals_matrix,file="PADB_dataset_level_replicability_analysis_data.RData")

acute_col2id = sapply(colnames(acute_pvals_matrix),function(x)strsplit(x,split=';')[[1]][1])
acute_col2data = sapply(acute_metadata[acute_col2id],function(x)paste(x$tissue,x$training,x$training_desc,sep=";"))
names(acute_col2data) = names(acute_col2id)

longterm_col2id = sapply(colnames(longterm_pvals_matrix),function(x)strsplit(x,split=';')[[1]][1])
longterm_col2data = sapply(longterm_metadata[longterm_col2id],function(x)paste(x$tissue,x$training,x$training_desc,sep=";"))
names(longterm_col2data) = names(longterm_col2id)

# Separate by tissue
rep_datasets = list()
rep_datasets[["acute_muscle"]] = acute_pvals_matrix[,
    grepl("muscle",acute_col2data) & 
      !grepl(pattern = "other|untra|control",acute_col2data)
]
rep_datasets[["acute_blood"]] = acute_pvals_matrix[,
    grepl("blood",acute_col2data) & 
    !grepl(pattern = "other|untra|control",acute_col2data)
]
rep_datasets[["longterm_muscle"]] = longterm_pvals_matrix[,
    grepl("muscle",longterm_col2data) & 
    !grepl(pattern = "other|untra|control",longterm_col2data)
    ]
rep_datasets[["longterm_blood"]] = longterm_pvals_matrix[,
    grepl("blood",longterm_col2data) & 
    !grepl(pattern = "other|untra|control",longterm_col2data)
    ]
sapply(rep_datasets,dim)
rep_datasets = lapply(rep_datasets,function(x){x[is.na(x)]=0.5;return(x)})
library(kernlab)
source('~/Desktop/screen/supplementary_data/submission_code/SCREEN_code_for_submission.R')
source('~/Desktop/screen/supplementary_data/submission_code/twogroups_methods_for_submission.R')
run_leadingeigen_clustering = function(x,cor_thr=0.2,toplot=F){
  x = x >= cor_thr
  mode(x) = 'numeric';diag(x)=0
  g = graph.adjacency(x,mode='undirected',weighted=T)
  if(toplot){plot(igraph::simplify(g))}
  return (cluster_infomap(g)$membership)
}
screen_res = sapply(rep_datasets,function(x)SCREEN(x,ks=2:ncol(x),nH=10000))
sapply(screen_res,function(x)colSums(x<=0.2))
acute_genes[screen_res$acute_muscle[,"SCREEN 9"]<=0.2]

rep_num_genes_per_percent = list()
rep_lfdr = 0.2
rep_gene_sets = list(); decay_plot_data = list()
par(mfrow=c(2,2))
for(nn in names(screen_res)){
  currgenes = acute_genes
  if(grepl("longterm",nn)){
    currgenes = longterm_genes
  }
  mat = screen_res[[nn]]
  num_genes = colSums(mat<=rep_lfdr)
  percents = (1:ncol(mat))/ncol(mat)
  plot(x=percents,y=num_genes,type='b')
  ind = which(percents >= 0.5)[1]
  selected_genes = currgenes[mat[,ind]<=rep_lfdr]
  if(length(selected_genes)>0){
    rep_gene_sets[[nn]] = selected_genes
  }
  decay_plot_data[[nn]] = list()
  decay_plot_data[[nn]][["percents"]] = percents
  decay_plot_data[[nn]][["num_genes"]] = num_genes
}

save(acute_genes,longterm_genes,screen_res,
     decay_plot_data,rep_gene_sets,
     file="PADB_dataset_level_replicability_analysis_results.RData")

