setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(metafor)
source('repos/motrpac/helper_functions.R')

# Get the datasets and their metadata
acute_datasets = get(load("PADB_univariate_results_and_preprocessed_data_acute.RData"))
acute_metadata = get(load("PADB_sample_metadata_acute.RData"))
longterm_datasets = get(load("PADB_univariate_results_and_preprocessed_data_longterm.RData"))
longterm_metadata = get(load("PADB_sample_metadata_longterm.RData"))

# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_pval <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  return(t.test(x1,x2,paired = T)$p.value)
}

get_ttest_pval_per_dataset<-function(mat,metadata){
  dataset_times = metadata$time[colnames(mat)]
  if(any(is.na(dataset_times))){
    mat = mat[,!is.na(dataset_times)]
    dataset_times = metadata$time[colnames(mat)]
  }
  min_time = min(dataset_times)
  other_times = setdiff(unique(dataset_times),min_time)
  times2pvals = c()
  for(other_time in other_times){
    curr_mat = mat[,dataset_times==min_time | dataset_times==other_time]
    curr_subjects = metadata$subject[colnames(curr_mat)]
    subjects_to_keep = names(which(table(curr_subjects)==2))
    curr_mat = curr_mat[,is.element(curr_subjects,set = subjects_to_keep)]
    ord = order(metadata$subject[colnames(curr_mat)],metadata$time[colnames(curr_mat)])
    curr_mat = curr_mat[,ord]
    curr_times = metadata$time[colnames(curr_mat)]
    paired_test_data = apply(curr_mat,1,get_paired_ttest_pval,sample2time=curr_times,t1=min_time,t2=other_time)
    times2pvals = cbind(times2pvals,paired_test_data)
    colnames(times2pvals)[ncol(times2pvals)] = as.character(other_time)
  }
  return(times2pvals)
}

acute_datasets_pvals = lapply(acute_datasets,function(x,y)get_ttest_pval_per_dataset(x$gene_data,y),y=acute_metadata)
longterm_datasets_pvals = lapply(longterm_datasets,function(x,y)get_ttest_pval_per_dataset(x$gene_data,y),y=longterm_metadata)
longterm_datasets_pvals = longterm_datasets_pvals[sapply(longterm_datasets_pvals,length)>0]
acute_datasets_pvals = acute_datasets_pvals[sapply(acute_datasets_pvals,length)>0]
sapply(acute_datasets_pvals,function(x)sapply(x,function(y)sum(is.na(y)|is.nan(y))))
sapply(longterm_datasets_pvals,function(x)sapply(x,function(y)sum(is.na(y)|is.nan(y))))

# in this script meta-analyses are done by gene
# we therefore need to create a table for each gene
# get dataset info - moderators
# tissue, training type, gse
get_dataset_moderators<-function(datasets,metadata){
  arrs = t(sapply(datasets,function(x)strsplit(x,split=';')[[1]][c(1,2,4)]))
  arrs[,2] = simplify_tissue_info(arrs[,2])
  arrs[grepl(datasets,pattern = 'control',ignore.case = T) | grepl(datasets,pattern="untrain",ignore.case = T),3] = "control"
  colnames(arrs) = c("gse","tissue","training")
  return(arrs)
}

# extract the data objects for the meta-analysis
acute_genes = rownames(acute_datasets[[1]]$gene_data)
for(i in 2:length(acute_datasets)){
  acute_genes = intersect(acute_genes,rownames(acute_datasets[[i]]$gene_data))
}
acute_mod = get_dataset_moderators(names(acute_datasets_pvals),acute_metadata)
longterm_genes = rownames(longterm_datasets[[1]]$gene_data)
for(i in 2:length(longterm_datasets)){
  longterm_genes = intersect(longterm_genes,rownames(longterm_datasets[[i]]$gene_data))
}
longterm_mod = get_dataset_moderators(names(longterm_datasets_pvals),longterm_metadata)

# # remove controls, untrained, and fat from the tables
# remove_undesired_datasets<-function(gdata){
#   gdata = gdata[!grepl("control",gdata$training),]
#   gdata = gdata[!grepl("other",gdata$training),]
#   gdata = gdata[!grepl("yoga",gdata$training),]
#   gdata = gdata[!grepl("untr",gdata$training),]
#   gdata = gdata[!grepl("fat",gdata$tissue),]
#   gdata = gdata[!grepl("adipose",gdata$tissue),]
#   return (gdata)
# }

# transform the datasets into matrices that can be plotted
acute_pvals_matrix = c()
for(d in names(acute_datasets_pvals)){
  l = acute_datasets_pvals[[d]]
  times = colnames(l)
  darr = strsplit(split=';',d)[[1]]
  darr[2] = simplify_tissue_info(darr[2])
  for (tt in times){
    currname = paste(c(tt,darr[c(2,4,1)]),collapse=';')
    acute_pvals_matrix = cbind(acute_pvals_matrix,l[acute_genes,as.character(tt)])
    colnames(acute_pvals_matrix)[ncol(acute_pvals_matrix)] = currname
  }
}
longterm_pvals_matrix = c()
for(d in names(longterm_datasets_pvals)){
  l = longterm_datasets_pvals[[d]]
  times = colnames(l)
  darr = strsplit(split=';',d)[[1]]
  darr[2] = simplify_tissue_info(darr[2])
  for (tt in times){
    currname = paste(c(tt,darr[c(2,4,1)]),collapse=';')
    longterm_pvals_matrix = cbind(longterm_pvals_matrix,l[longterm_genes,as.character(tt)])
    colnames(longterm_pvals_matrix)[ncol(longterm_pvals_matrix)] = currname
  }
}

save(longterm_pvals_matrix,acute_pvals_matrix,acute_datasets_pvals,longterm_datasets_pvals,
     file="PADB_dataset_level_replicability_analysis_data.RData")

# Separate by tissue
rep_datasets = list()
rep_datasets[["acute_muscle"]] = acute_pvals_matrix[,
    grepl("muscle",colnames(acute_pvals_matrix)) & 
      !grepl(pattern = "other|untra|control",colnames(acute_pvals_matrix))
]
rep_datasets[["acute_blood"]] = acute_pvals_matrix[,
    grepl("blood",colnames(acute_pvals_matrix)) & 
    !grepl(pattern = "other|untra|control",colnames(acute_pvals_matrix))
]
rep_datasets[["longterm_muscle"]] = longterm_pvals_matrix[,
    grepl("muscle",colnames(longterm_pvals_matrix)) & 
    !grepl(pattern = "other|untra|control",colnames(longterm_pvals_matrix))
    ]
rep_datasets[["longterm_blood"]] = longterm_pvals_matrix[,
    grepl("blood",colnames(longterm_pvals_matrix)) & 
    !grepl(pattern = "other|untra|control",colnames(longterm_pvals_matrix))
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

