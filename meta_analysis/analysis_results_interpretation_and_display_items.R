# This script interprets the results of the acute and longterm analyses
# The script is structured by generating one display item or analysis at a time from scratch

setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery');library(corrplot)
source('helper_functions.R')
library(org.Hs.eg.db)
entrez2symbol = as.list( org.Hs.egSYMBOL)

# Classification overall results using the gene intensity data
load("Longterm_data_analysis_expression_intensities_data_classification_tests.RData")
longterm_classification_perf = main_configurations_lso_perf_scores
load("Acute_data_analysis_expression_intensities_data_classification_tests.RData")
acute_classification_perf = main_configurations_lso_perf_scores
acute_rocs = sapply(acute_classification_perf,function(x)mean(x$aucs[,3]))
longterm_rocs = sapply(longterm_classification_perf,function(x)mean(x$aucs[,3]))
acute_auprs = sapply(acute_classification_perf,function(x)mean(x$aucs[,2]))
longterm_auprs = sapply(longterm_classification_perf,function(x)mean(x$aucs[,2]))

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
mixed_effects_gene_sets_topgo = run_topgo_enrichment_fisher(mixed_effects_gene_sets,all_genes)
extract_top_go_results(mixed_effects_gene_sets_topgo,0.1)
# GSEA
pathways = reactomePathways(names(longterm_bic))
system.time({longterm_gsea_res = fgsea_wrapper(pathways,longterm_bic,nperm=10000,minSize = 5)})
qvals = longterm_gsea_res$padj; ES = longterm_gsea_res$ES
longterm_selected_pathways = longterm_gsea_res$pathway[qvals<0.01 & ES > 0]
pathways = reactomePathways(names(acute_bic))
system.time({acute_gsea_res = fgsea_wrapper(pathways,acute_bic,nperm=10000,minSize = 5)})
qvals = acute_gsea_res$padj; ES = acute_gsea_res$ES
acute_selected_pathways = acute_gsea_res$pathway[qvals<0.01 & ES > 0]

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
# rep gene sets
rep_gene_sets_topgo = run_topgo_enrichment_fisher(rep_gene_sets,all_genes)
extract_top_go_results(rep_gene_sets_topgo)

# particularily interesting gene sets
unlist(entrez2symbol[intersect(longterm_rep_gene_list,mixed_effects_gene_sets$longterm)])
unlist(entrez2symbol[mixed_effects_gene_sets$both])
unlist(entrez2symbol[rep_gene_sets$both])

# heatmaps
library(gplots)
get_tstats_heatmap<-function(xx,remove_rows=T,entrez2symbol=NULL,max_t = 6, min_t=2,genes_as_rows=F,...){
  xx[is.na(xx)|is.nan(xx)]=0
  xx[xx>max_t]=max_t;xx[xx< -(max_t)] = -max_t
  xx[xx<min_t & xx > -min_t]=0
  if(remove_rows){xx = xx[!apply(xx==0,1,all),]}
  if(!is.null(entrez2symbol)){rownames(xx) = unlist(entrez2symbol[rownames(xx)])}
  if(genes_as_rows){
    heatmap.2(xx,scale="none",trace="none",col=colorRampPalette(c("red","white","blue"))(256),...)
  }else{
    heatmap.2(t(xx),scale="none",trace="none",col=colorRampPalette(c("red","white","blue"))(256),...)
  }
}
get_tstats_heatmap(longterm_ts[mixed_effects_gene_sets$both,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 1)
get_tstats_heatmap(longterm_ts[mixed_effects_gene_sets$both,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 2)
get_tstats_heatmap(longterm_ts[mixed_effects_gene_sets$both,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 2,genes_as_rows=T)
get_tstats_heatmap(longterm_ts[mixed_effects_gene_sets$longterm,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 2,genes_as_rows=F)
get_tstats_heatmap(longterm_ts[rep_gene_sets$longterm,],
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 2,genes_as_rows=F)

pathways = reactomePathways(all_genes)
all_pathway_genes = unique(unlist(pathways))
non_pathway_genes =sapply(rep_gene_sets,setdiff,y=all_pathway_genes)
non_pathway_genes = sapply(non_pathway_genes,function(x,y)unlist(y[x]),y=entrez2symbol)

xx1 = acute_ts[mixed_effects_gene_sets$both,]
xx2 = longterm_ts[mixed_effects_gene_sets$both,]
get_tstats_heatmap(cbind(xx1,xx2),
                   entrez2symbol=entrez2symbol,mar=c(8,15),min_t = 0.2)

