
# This script interprets the results of the acute and longterm analyses.
# The summarized analyses here are the mixed effects, replicability, and dataset level analyses.
# The script is structured by generating one display item or analysis at a time from scratch.

setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(corrplot);library(metafor);library(gplots);library(Vennerable)
source('/Users/David/Desktop/repos/motrpac/rscripts/helper_functions.R')
source('/Users/David/Desktop/repos/motrpac/rscripts/helper_functions_meta_analaysis.R')
library(org.Hs.eg.db)
entrez2symbol = as.list(org.Hs.egSYMBOL)

# # Vizualize overlaps
# source("https://bioconductor.org/biocLite.R")
# biocLite("RBGL")
# install.packages("reshape")
# library(RBGL)
# install.packages("Vennerable", repos="http://R-Forge.R-project.org",dependencies = T)

######################################################################
# functions
######################################################################
library(gplots)
get_tstats_heatmap<-function(xx,remove_rows=T,entrez2symbol=NULL,max_t = 10, min_t=1,genes_as_rows=F,...){
  xx[is.na(xx)|is.nan(xx)]=0
  xx[xx>max_t]=max_t;xx[xx< -(max_t)] = -max_t
  xx[xx<min_t & xx > -min_t]=0
  if(remove_rows){xx = xx[!apply(xx==0,1,all),]}
  if(!is.null(entrez2symbol)){rownames(xx) = unlist(entrez2symbol[rownames(xx)])}
  if(genes_as_rows){
    heatmap.2(xx,scale="none",trace="none",col=colorRampPalette(c("green","white","red"))(256),...)
  }else{
    heatmap.2(t(xx),scale="none",trace="none",col=colorRampPalette(c("green","white","red"))(256),...)
  }
}
plot_with_err_bars<-function(xnames,avg,sdev,add=F,...){
  if(add){
    lines(avg,pch=19,...)
  }
  else{
    plot(avg,xaxt = "n",pch=19, type='l',...)
    axis(1, at=1:length(xnames), labels=xnames)
  }
  # hack: we draw arrows but with very special "arrowheads"
  arrows(1:length(xnames), avg-sdev, 1:length(xnames), avg+sdev, length=0.05, angle=90, code=3)
}
run_corrmat_clustering_of_a_gene_set<-function(genes,m,exclude_cols_regex = "fat|treatment|yoga",
                                               min_homogn=0.5,cor_method="spearman",min_clust_size=5){
  m = m[genes,]
  m = m[,!grepl(exclude_cols_regex,colnames(m))]
  mc = cor(t(m),method = "spearman")
  m_clust = cluster_genes_by_homogeneity(mc,min_homogn,kmeans)
  selected_clusters = names(which(table(m_clust)>=min_clust_size))
  m_clust = m_clust[is.element(m_clust,set=selected_clusters)]
  return(list(clustering=m_clust,homogeneities=cluster_homogeneities(mc,m_clust),mc=mc))
}
######################################################################
######################################################################
######################################################################

############ Display item: cohort tables ############
# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
load("PADB_sample_metadata_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
acute_subj = sample_metadata$subject
acute_time = sample_metadata$time
data_gsms = unique(unlist(sapply(acute_metadata,function(x)x$gsms)))
write.table(t(table(sample_metadata$tissue[data_gsms],sample_metadata$training[data_gsms])),sep="\t",quote=F)
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
load("PADB_sample_metadata_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
longterm_subj = sample_metadata$subject
longterm_time = sample_metadata$time
data_gsms = unique(unlist(sapply(longterm_metadata,function(x)x$gsms)))
write.table(t(table(sample_metadata$tissue[data_gsms],sample_metadata$training[data_gsms])),sep="\t",quote=F)
get_metadata_for_table<-function(x,subjs,times){
  gsms = x$gsms
  N = length(gsms)
  S = subjs[gsms]
  Ns = length(unique(S))
  Tset = unique(times[gsms])
  if(min(Tset)<0){Tset[Tset<0] = "Pre"}
  else{Tset[Tset<=0] = "Pre"}
  v = c(x$gse,x$tissue,x$training,N,Ns,paste(Tset,collapse=','),x$additional_info)
  names(v) = c("Dataset","Tissue","Training","Nsample","Nsubject","Time points","AdditionalInfo")
  return(v)
}
acute_cohort_table = t(sapply(acute_metadata,get_metadata_for_table,subjs=acute_subj,times=acute_time))
longterm_cohort_table = t(sapply(longterm_metadata,get_metadata_for_table,subjs=longterm_subj,times=longterm_time))
write.table(acute_cohort_table,quote=F,sep="\t")
write.table(longterm_cohort_table,quote=F,sep="\t")

# Statistics for display items
# gd = acute_gene_tables_raw[[1]]
# tb = table(gd$time,gd$tissue)
# gd$time = simplify_time_acute(gd$time)
# tb = table(gd$time,gd$tissue)
# gd = longterm_gene_tables_raw[[1]]
# tb = table(gd$time,gd$tissue)

#####################################################

# Load or compute the gene patterns
# Gene patterns for the analysis
# # LOAD data structures needed for the lines below: ADD
# weighted_avg_matrices=list()
# weighted_avg_matrices[["acute"]] = t(sapply(acute_gene_tables_raw,get_gene_weighted_avg_pattern))
# weighted_avg_matrices[["longterm"]] = t(sapply(longterm_gene_tables_raw,get_gene_weighted_avg_pattern))
# weighted_avg_matrices = lapply(weighted_avg_matrices,reorder_weighted_avg_matrix)
# save(weighted_avg_matrices,file="effect_weighted_avg_matrices_by_training_and_tissues.RData")
# Or load
load("effect_weighted_avg_matrices_by_training_and_tissues.RData")

# Load all three analysis results
# Meta-analysis: specific time points
load("tp_meta_analysis_results.RData")
tp_meta_analysis_gene_sets = selected_genes_all_tests
tp_meta_analysis_gene_names = sapply(tp_meta_analysis_gene_sets,function(x,y)unlist(y[x]),y=entrez2symbol)
sapply(tp_meta_analysis_gene_names,function(x)x[grepl("^MY|^COL",x)])
# # Meta-analysis: all data with binary moderators or simple random effects models
# Deprecated for now: may be highly confounded
# load("PADB_metafor_simple_random_effect_results.RData")
# Replication analysis
load("PADB_dataset_level_replicability_analysis_results.RData")
rep_gene_sets = rep_gene_sets_0.4

window_based_metadnalysis_raw_results = tp_meta_analysis_results
selected_genes = selected_genes_all_tests_names
topgo_enrichments = selected_genes_all_tests_topgo
save(window_based_metadnalysis_raw_results,selected_genes,topgo_enrichments,file="time_window_meta_analysis_results.RData")

# Some overlaps
names(rep_gene_sets)
names(tp_meta_analysis_gene_sets)
ll = c(rep_gene_sets[c(1,3,5,7,8)],tp_meta_analysis_gene_sets)
names(ll)[1:5] = paste("REP",names(ll)[1:5],sep="_")
names(ll)[6:9] = paste("META",names(ll)[6:9],sep="_")
V = Venn(ll[c(2,3,6)])
plot(V,doWeights=F)
V = Venn(ll[c(4,5,8)])
plot(V,doWeights=F)
V = Venn(ll[c(1,7)])
plot(V,doWeights=F)
expression_gene_sets = ll
expression_gene_sets = lapply(expression_gene_sets,function(x,y)sort(unlist(y[x])),y=entrez2symbol)
save(expression_gene_sets,file="Exercise_data_analysis_gene_sets.RData")

# cluster the gene sets by their patterns
selected_gene_clustering = list()
for(i in 1:length(tp_meta_analysis_gene_sets)){
  nn = names(tp_meta_analysis_gene_sets)[i]
  m = weighted_avg_matrices$acute
  if(grepl("longterm",nn)){m = weighted_avg_matrices$longterm}
  m1 = m[,grepl("muscle",colnames(m))]
  if(grepl("blood",nn)){m1 = m[,grepl("blood",colnames(m))]}
  m1[is.na(m1)] = 0
  curr_genes = tp_meta_analysis_gene_sets[[nn]]
  curr_rep_genes = unique(unlist(rep_gene_sets[grepl(nn,names(rep_gene_sets))]))
  curr_genes = union(curr_genes,curr_rep_genes)
  curr_genes = intersect(curr_genes,rownames(m1)) # REQUIRED since the rep analysis contains more genes
  #curr_genes = intersect(curr_rep_genes,rownames(m1)) # test
  selected_gene_clustering[[nn]] = run_corrmat_clustering_of_a_gene_set(curr_genes,m1)
}
lapply(selected_gene_clustering,function(x)table(x[[1]]))
selected_gene_clusters=list()
for(i in 1:length(tp_meta_analysis_gene_sets)){
  nn = names(tp_meta_analysis_gene_sets)[i]
  m_clust = selected_gene_clustering[[nn]][[1]]
  for(j in sort(unique(m_clust))){
    if (sum(m_clust==j)<10){next}
    selected_gene_clusters[[paste(nn,j,sep="_")]] = names(which(m_clust==j))
  }
}
sapply(selected_gene_clusters,length)
selected_gene_clusters_names = lapply(selected_gene_clusters,function(x,y)unlist(y[x]),y=entrez2symbol)
sapply(selected_gene_clusters_names,intersect,y=known_genes)

# Enrichments
all_genes = unique(unlist(sapply(weighted_avg_matrices,rownames)))
selected_gene_clusters_topgo_res = run_topgo_enrichment_fisher(selected_gene_clusters,all_genes)
extract_top_go_results(selected_gene_clusters_topgo_res)
get_most_sig_enrichments_by_groups(extract_top_go_results(selected_gene_clusters_topgo_res),2)
table(extract_top_go_results(selected_gene_clusters_topgo_res)[,1])

# cluster plots
names(selected_gene_clusters)
par(mfrow=c(3,2))
for(i in c(1:3)){
  m_1 = weighted_avg_matrices$acute[selected_gene_clusters[[i]],]
  m_1[is.na(m_1)|is.nan(m_1)]=0
  m_1 = m_1[,!grepl("treatment",colnames(m_1))]
  plot_gene_pattern(apply(m_1,2,mean),errs = apply(m_1,2,sd),tosmooth = T,mfrow=NULL,y_lim_add = 0.5,y_lim_min = 1)
}

# cluster plots
names(selected_gene_clusters)
par(mfrow=c(3,2))
for(i in c(6:8)){
  m_1 = weighted_avg_matrices$longterm[selected_gene_clusters[[i]],]
  m_1[is.na(m_1)|is.nan(m_1)]=0
  m_1 = m_1[,!grepl("treatment",colnames(m_1))]
  plot_gene_pattern(apply(m_1,2,mean),errs = apply(m_1,2,sd),tosmooth = T,mfrow=NULL,y_lim_add = 0.5,y_lim_min = 0.5)
}
for(i in c(9:11)){
  m_1 = weighted_avg_matrices$longterm[selected_gene_clusters[[i]],]
  m_1[is.na(m_1)|is.nan(m_1)]=0
  m_1 = m_1[,!grepl("treatment",colnames(m_1))]
  plot_gene_pattern(apply(m_1,2,mean),errs = apply(m_1,2,sd),tosmooth = T,mfrow=NULL,y_lim_add = 0.5,y_lim_min = 0.5)
}

# cluster plots
names(selected_gene_clusters)
par(mfrow=c(2,2))
for(i in c(4:5)){
  m_1 = weighted_avg_matrices$acute[selected_gene_clusters[[i]],]
  m_1[is.na(m_1)|is.nan(m_1)]=0
  m_1 = m_1[,!grepl("treatment",colnames(m_1))]
  plot_gene_pattern(apply(m_1,2,mean),errs = apply(m_1,2,sd),tosmooth = T,mfrow=NULL,y_lim_add = 0.5,y_lim_min = 1)
}

# Selected intersections
names(tp_meta_analysis_gene_sets)
t(t(unname(unlist(entrez2symbol[intersect(tp_meta_analysis_gene_sets[[1]],tp_meta_analysis_gene_sets[[3]])]))))

load("PADB_dataset_level_meta_analysis_data.RData")
par(mfrow=c(1,1))
get_subset_forest_plot(acute_gene_tables[["70"]],tissue = "muscle")
par(mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$acute["70",],tosmooth = T,mfrow=NULL)
plot_gene_pattern(weighted_avg_matrices$longterm["70",],tosmooth = T,mfrow=NULL)

gene = "4627"
gdata = acute_gene_tables[[gene]]
gdata = gdata[gdata$time<10,]
get_subset_forest_plot(gdata ,tissue = "muscle")
rma.mv(yi,vi,random=~V1|gse,subset = tissue=="muscle",data=gdata,mods = ~training)

