########## !!!!!!!!!! #########
# TODO:
# A major change was performed
# Commented parts need to be revised
###############################

# This script interprets the results of the acute and longterm analyses.
# The summarized analyses here are the mixed effects, replicability, and dataset level analyses.
# The script is structured by generating one display item or analysis at a time from scratch.

setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery');library(corrplot);library(metafor)
source('repos/motrpac/helper_functions.R')
source('repos/motrpac/helper_functions_meta_analaysis.R.R')
library(org.Hs.eg.db)
entrez2symbol = as.list(org.Hs.egSYMBOL)
######################################################################
# heatmaps
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
# ######################################################################
# # Get dataset sizes
# load("PADB_univariate_results_and_preprocessed_data_acute.RData")
# acute_study_sizes = sapply(dataset2preprocessed_data,function(x)ncol(x$gene_data))
# load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
# longterm_study_sizes = sapply(dataset2preprocessed_data,function(x)ncol(x$gene_data))
# dataset_sizes = c(acute_study_sizes,longterm_study_sizes)
# 
# # Dataset sizes vs. t-tests: show effect sizes vs. dataset sizes and precision
# # Load the ts
# load("Longterm_replicability_analysis.RData")
# longterm_ts = rep_data$dataset_stats
# longterm_ps = rep_data$dataset_pvals
# load("Acute_replicability_analysis.RData")
# acute_ts = rep_data$dataset_stats
# acute_ps = rep_data$dataset_pvals
# shared_genes = intersect(rownames(longterm_ts),rownames(acute_ts))
# all_ts = cbind(longterm_ts[shared_genes,],acute_ts[shared_genes,])
# all_ts[is.na(all_ts)|is.nan(all_ts)]=0
# all_ps = cbind(longterm_ps[shared_genes,],acute_ps[shared_genes,])
# all_ps[is.na(all_ps)|is.nan(all_ps)]=0.5
# ts_info_table = sapply(colnames(all_ts),function(x)strsplit(x,split=';')[[1]])
# ts_info_table = rbind(c(rep("longterm",ncol(longterm_ts)),rep("acute",ncol(acute_ts))),ts_info_table)
# 
# # correlation of effects vs study size
# dataset_sizes_info = sapply(names(dataset_sizes),function(x)strsplit(x,split=';')[[1]][1:4])
# dataset_sizes_info = rbind(dataset_sizes_info,dataset_sizes)
# dataset_sizes_info[2,] = sapply(dataset_sizes_info[2,],simplify_tissue_info)
# gse_tissue_to_size = as.numeric(dataset_sizes_info[5,])
# names(gse_tissue_to_size) = apply(dataset_sizes_info[c(1,2),],2,paste,collapse=",")
# all_ts_col_to_sizes = as.numeric(gse_tissue_to_size[apply(ts_info_table[c(5,3),],2,paste,collapse=",")])
# ts_sizes_corrs = cor(abs(t(all_ts)),all_ts_col_to_sizes,method="spearman")[,1]
# hist(ts_sizes_corrs,main="Correlation between t-statistic and sample size",breaks=100)
# abline(v=0,col="red",lwd=5,lty=2)

# Load effect sizes
load("PADB_dataset_level_meta_analysis_data.RData")
# ts1 = acute_ts
# ts2 = acute_effects_matrix/acute_sds_matrix
# x1 = ts1[,"0;blood;endurance;GSE46075"]
# x2 = ts2[,"0;blood;endurance;GSE46075"]
# plot(-x1,x2);abline(0,1)
# corrs1 = c()
# for(j in 1:nrow(acute_effects_matrix)){
#   corrs1[j] = cor(longterm_effects_matrix[j,],longterm_sds_matrix[j,],method="spearman")
# }
# hist(corrs1,breaks=100)

# Load the gene sets from the different analyses
load("metafor_gene_sets.RData")
load("PADB_dataset_level_replicability_analysis_results.RData")
metafor_gene_sets = metafor_gene_sets[[1]]
all_selected_pe_genes = unique(unlist(metafor_gene_sets))
all_selected_pe_genes = union(all_selected_pe_genes,unlist(rep_gene_sets))
all_selected_pe_genes = sort(all_selected_pe_genes)
gene_data_summary_table = c()
for(nn in names(rep_gene_sets)){
  v = rep(F,length(all_selected_pe_genes));names(v)=all_selected_pe_genes
  v[rep_gene_sets[[nn]]] = T
  gene_data_summary_table = cbind(gene_data_summary_table,v)
  colnames(gene_data_summary_table)[ncol(gene_data_summary_table)] = nn
}
colnames(gene_data_summary_table) = paste("rep_",colnames(gene_data_summary_table),sep="")
for(nn in names(metafor_gene_sets)){
  v = rep(F,length(all_selected_pe_genes));names(v)=all_selected_pe_genes
  v[metafor_gene_sets[[nn]]] = T
  gene_data_summary_table = cbind(gene_data_summary_table,v)
  colnames(gene_data_summary_table)[ncol(gene_data_summary_table)] = nn
}
gene_data_summary_table = gene_data_summary_table[,-ncol(gene_data_summary_table)]
#write.table(gene_data_summary_table,file="replicability_and_metaanalysis_genes.txt",sep="\t",quote=F)
sort(unlist(entrez2symbol[rownames(gene_data_summary_table)]))

# Get the meta-analysis p-values
selected_genes_meta_analysis_pvals = c()
for(g in rownames(gene_data_summary_table)){
  gname = entrez2symbol[g]
  g_acute_pvals = rep(NA,ncol(acute_ps));names(g_acute_pvals) = colnames(acute_ps)
  if(g %in% rownames(acute_ps)){g_acute_pvals = acute_ps[g,]}
  names(g_acute_pvals) = paste("acute",names(g_acute_pvals),sep=",")
  g_longterm_pvals = longterm_ps[g,]
  names(g_longterm_pvals) = paste("longterm",names(g_longterm_pvals),sep=",")
  selected_genes_meta_analysis_pvals = rbind(selected_genes_meta_analysis_pvals,c(g_acute_pvals,g_longterm_pvals))
}
selected_genes_meta_analysis_pvals[is.na(selected_genes_meta_analysis_pvals)] = 1
rownames(selected_genes_meta_analysis_pvals) = rownames(gene_data_summary_table)

selected_genes_meta_analysis_mod_pvals = list()
for(g in rownames(gene_data_summary_table)){
  acute_res = sapply(acute_meta_analysis_results,function(x)x[[g]])
  acute_res = acute_res[grepl(names(acute_res),pattern="mixed")]
  acute_mod_ps = unlist(sapply(acute_res,get_modifiers_pvals))
  names(acute_mod_ps) = gsub(names(acute_mod_ps),pattern = "(mixed_effects_)|(pval_)",replace="")
  longterm_res = sapply(longterm_meta_analysis_results,function(x)x[[g]])
  longterm_res = longterm_res[grepl(names(longterm_res),pattern="mixed")]
  longterm_mod_ps = unlist(sapply(longterm_res,get_modifiers_pvals))
  names(longterm_mod_ps) = gsub(names(longterm_mod_ps),pattern = "(mixed_effects_)|(pval_)",replace="")
  selected_genes_meta_analysis_mod_pvals[[g]] = list(acute = acute_mod_ps, longterm = longterm_mod_ps)
}

# Look at specific examples
load("PADB_metafor_meta_analysis_results.RData")
sort(rowSums(gene_data_summary_table))

# Common examples that we expect to see
gene = "10891" # PGC1 gene
gene = "6678"

gname = entrez2symbol[gene]
sig_ps = which(gene_data_summary_table[as.character(gene),])
acute_res = sapply(acute_meta_analysis_results,function(x)x[[gene]])
longterm_res = sapply(longterm_meta_analysis_results,function(x)x[[gene]])
acute_gene_table = acute_gene_tables[[gene]]
get_subset_forest_plot(acute_gene_table,"muscle")
get_subset_forest_plot(acute_gene_table,"blood")
longterm_gene_table = longterm_gene_tables[[gene]]
get_subset_forest_plot(longterm_gene_table,"muscle")
get_subset_forest_plot(longterm_gene_table,"blood")

get_subset_forest_plot<-function(gdata,tissue="all",training="all",sortby = "time"){
  if(!tissue=="all"){
    gdata = gdata[gdata$tissue==tissue,]
  }
  if(!training=="all"){
    gdata = gdata[gdata$training==training,]
  }
  ord = order(gdata[,sortby])
  gdata = gdata[ord,]
  res = rma.mv(yi,vi,random = ~ 1|gse,data=gdata)
  forest(res,slab=paste(gdata$gse,gdata$time,sep=','))
}

# Show heatmaps
get_tstats_heatmap(longterm_effects_matrix[metafor_gene_sets$`longterm,random_effects_muscle`,],
                   entrez2symbol=entrez2symbol,mar=c(12,15),min_t = 0)
get_tstats_heatmap(acute_effects_matrix[metafor_gene_sets$`acute,random_effects_muscle`,grepl("muscle",colnames(acute_effects_matrix))],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)
get_tstats_heatmap(acute_effects_matrix[metafor_gene_sets$`acute,mixed_effects_muscle`,grepl("muscle",colnames(acute_effects_matrix))],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)
get_tstats_heatmap(acute_effects_matrix[metafor_gene_sets$`acute,muscle,modifiers`,grepl("muscle",colnames(acute_effects_matrix))],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)

top_acue_avg = names(sort(abs(rowMeans(acute_effects_matrix)),decreasing = T)[1:20])
get_tstats_heatmap(acute_effects_matrix[top_acue_avg,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)
top_longterm_avg = names(sort(abs(rowMeans(longterm_effects_matrix)),decreasing = T)[1:20])
get_tstats_heatmap(longterm_effects_matrix[top_longterm_avg,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)

# Classification overall results using the gene intensity data
# load("Longterm_data_analysis_expression_intensities_data_classification_tests.RData")
# longterm_classification_perf = main_configurations_lso_perf_scores
# load("Acute_data_analysis_expression_intensities_data_classification_tests.RData")
# acute_classification_perf = main_configurations_lso_perf_scores
# acute_rocs = sapply(acute_classification_perf,function(x)mean(x$aucs[,3]))
# longterm_rocs = sapply(longterm_classification_perf,function(x)mean(x$aucs[,3]))
# acute_auprs = sapply(acute_classification_perf,function(x)mean(x$aucs[,2]))
# longterm_auprs = sapply(longterm_classification_perf,function(x)mean(x$aucs[,2]))

# Compare mixed effects analyses
load("Longterm_data_analysis_mixed_effect_analysis.RData")
longterm_bic = mixed_effect_pvals_apprx_bic_diff
load("Acute_data_analysis_mixed_effect_analysis.RData")
acute_bic = mixed_effect_pvals_apprx_bic_diff
shared_genes = intersect(names(longterm_bic),names(acute_bic))
plot(acute_bic[shared_genes],longterm_bic[shared_genes])
# Gene enrichment
mixed_effects_gene_sets = list()
mixed_effects_gene_sets[["longterm"]] = names(which(longterm_bic>0))
mixed_effects_gene_sets[["acute"]] = names(which(acute_bic>0))
mixed_effects_gene_sets[["both"]] = intersect(mixed_effects_gene_sets[[1]],mixed_effects_gene_sets[[2]])
sapply(mixed_effects_gene_sets,length)
all_genes = union(names(longterm_bic),names(acute_bic))
# mixed_effects_gene_sets_topgo = run_topgo_enrichment_fisher(mixed_effects_gene_sets,all_genes)
# extract_top_go_results(mixed_effects_gene_sets_topgo,0.1)
# # GSEA
# pathways = reactomePathways(names(longterm_bic))
# system.time({longterm_gsea_res = fgsea_wrapper(pathways,longterm_bic,nperm=10000,minSize = 5)})
# qvals = longterm_gsea_res$padj; ES = longterm_gsea_res$ES
# longterm_selected_pathways = longterm_gsea_res$pathway[qvals<0.01 & ES > 0]
# pathways = reactomePathways(names(acute_bic))
# system.time({acute_gsea_res = fgsea_wrapper(pathways,acute_bic,nperm=10000,minSize = 5)})
# qvals = acute_gsea_res$padj; ES = acute_gsea_res$ES
# acute_selected_pathways = acute_gsea_res$pathway[qvals<0.01 & ES > 0]
# Analysis of replicated results
rep_gene_sets = list()
load("Longterm_replicability_analysis.RData")
colSums(screen_res<=0.2)
longterm_ts = rep_data$dataset_stats
longrerm_list_screen = screen_res
rep_gene_sets[["longterm"]] = rownames(rep_data$dataset_pvals)[screen_res[,"SCREEN 5"]<=0.2]
intersect(rep_gene_sets[["longterm"]],mixed_effects_gene_sets$longterm)
intersect(rep_gene_sets[["longterm"]],mixed_effects_gene_sets$both)
load("Acute_replicability_analysis.RData")
acute_ts = rep_data$dataset_stats
acute_list_screen = screen_res
rep_gene_sets[["acute"]] = rownames(rep_data$dataset_pvals)[screen_res[,"SCREEN 12"]<=0.2]
intersect(rep_gene_sets[["acute"]],mixed_effects_gene_sets$acute)
intersect(rep_gene_sets[["acute"]],mixed_effects_gene_sets$both)
rep_gene_sets[["both"]] = intersect(rep_gene_sets$longterm,rep_gene_sets$acute)
sapply(rep_gene_sets,length)
# # rep gene sets
# rep_gene_sets_topgo = run_topgo_enrichment_fisher(rep_gene_sets,all_genes)
# extract_top_go_results(rep_gene_sets_topgo)
# all_diff_genes = unique(c(unlist(mixed_effects_gene_sets),unlist(rep_gene_sets)))
# ks.test(ts_sizes_corrs,ts_sizes_corrs[all_diff_genes])
get_tstats_heatmap(longterm_effects_matrix[mixed_effects_gene_sets$longterm,],entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)
get_tstats_heatmap(all_ps[mixed_effects_gene_sets$both,],entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0)

xx = all_ts[intersect(all_diff_genes,rownames(all_ts)),]
colnames(xx) = apply(ts_info_table[c(3,5),],2,paste,collapse=",")
get_tstats_heatmap(xx,entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0,max_t=5)

classes = apply(ts_info_table[c(3,4),],2,paste,collapse=",")
num_cls = sum(table(classes)>1)
par(mfrow=c(2,num_cls/2))
num_clusters=3
for(cl in classes){
  cols = classes == cl
  if(sum(cols)<=1){next}
  curxx = xx[,cols]
  clusters = kmeans(curxx,num_clusters)
  clusters = clusters$cluster
  cluster_col = heat.colors(num_clusters)
  ord = order(as.numeric(ts_info_table[2,cols]))
  colnames(curxx) = ts_info_table[2,cols]
  times = as.numeric(colnames(curxx))[ord]
  for(j in unique(clusters)){
    clust_inds = clusters==j
    plot_with_err_bars(times,colMeans(curxx[clust_inds,]),apply(curxx[clust_inds,],2,sd),
                       add=(j!=1),ylim=c(-12,12),main=cl,col=cluster_col[j])
  }
}

pathways = reactomePathways(all_genes)
all_pathway_genes = unique(unlist(pathways))
non_pathway_genes =sapply(rep_gene_sets,setdiff,y=all_pathway_genes)
non_pathway_genes = sapply(non_pathway_genes,function(x,y)unlist(y[x]),y=entrez2symbol)

# Simple unsupervised analysis: biclustering
source('http://acgt.cs.tau.ac.il/twigs/TWIGS_gibbs_sampling_v1.R')
library(isa2)
run_isa<-function(I,min_row,min_col,direction = c('updown','updown'),
                  thr.row = c(1.5,2,2.5),thr.col = c(1.5,2,2.5),no.seeds = 1000,...){
  isa.result=isa(I,direction = direction,thr.row = thr.row,thr.col = thr.col,no.seeds = no.seeds)
  rows_mat = isa.result$rows
  cols_mat = isa.result$columns
  to_keep = apply(rows_mat>0,2,sum) >= min_row & apply(cols_mat>0,2,sum) >= min_col
  rows_mat = rows_mat[,to_keep]
  cols_mat = cols_mat[,to_keep]
  rows = list(); cols = list()
  if (ncol(rows_mat)==0){return (list(rows=rows,cols=cols))}
  for (j in 1:ncol(rows_mat)){
    rows[[j]] = rownames(I)[rows_mat[,j]>0]
    cols[[j]] = colnames(I)[cols_mat[,j]>0]
  }
  return (list(rows=rows,cols=cols))
}
all_ts[is.na(all_ts)|is.nan(all_ts)] = 0
isa_res = run_isa(all_ts,5,5,thr.row = c(1,2,3),thr.col = c(1,2,3),no.seeds=250)
isa_res = reduce_overlaps_debi(isa_res$rows,isa_res$cols,0.8,50)
length(isa_res[[1]])
sapply(isa_res$cols,length)
sapply(isa_res$rows,length)
isa_res$cols
topgo_bic_enrichments = run_topgo_enrichment_fisher(isa_res$rows,rownames(all_ts),go_dags=c("BP"))
ps = as.numeric(topgo_bic_enrichments$classicFisher)
ps[is.na(ps)] = 1e-30
qs = p.adjust(ps)
table(qs<0.1)
extracted_res = topgo_bic_enrichments[qs<0.1,]

# get some stats of the bics
bic_mats = sapply(1:length(isa_res[[1]]),function(x,y,z,m)m[y[[x]],z[[x]]],
                       y=isa_res$rows,z=isa_res$cols,m=all_ts)
sapply(bic_mats,function(x)colMeans(abs(x)))
sapply(sapply(bic_mats,rowMeans),min)
get_tstats_heatmap(bic_mats[[2]][1:100,],min_t = 0, max_t = 10)

# bimax
all_ps = cbind(longterm_ps[shared_genes,],acute_ps[shared_genes,])
colnames(all_ps)=colnames(all_ts)
num_ts = apply(all_ps,2,p.adjust,method='fdr') < 0.1
mode(num_ts) = 'numeric'
num_ts[is.na(num_ts)] = 0
table(num_ts)
library(biclust)
bimax_sol = run_bimax(num_ts,1000,10,20)
bimax_sol = reduce_overlaps_debi(bimax_sol$rows,bimax_sol$cols,0.6,100)
length(bimax_sol[[1]])
sapply(bimax_sol$cols,length)
sapply(bimax_sol$rows,length)
bic_mats = sapply(1:length(bimax_sol[[2]]),function(x,y,z,m)m[y[[x]],z[[x]]],
                  y=bimax_sol$rows,z=bimax_sol$cols,m=all_ts)
get_tstats_heatmap(bic_mats[[1]],min_t = 2, max_t = 10)

# plaid models
plaid_res = biclust(all_ts, method=BCPlaid(), cluster="b", fit.model = y ~ m + b,
        background = TRUE, background.layer = NA, background.df = 1, row.release = 0.7,
        col.release = 0.7, shuffle = 3, back.fit = 0, max.layers = 20, iter.startup = 5,
        iter.layer = 10, verbose = TRUE)
plaid_bics = bicluster(all_ts,plaid_res)
plaid_bics = plaid_bics[sapply(plaid_bics,length)>0]
sapply(plaid_bics,dim)
get_tstats_heatmap(plaid_bics[[3]],min_t = 2, max_t = 10)
