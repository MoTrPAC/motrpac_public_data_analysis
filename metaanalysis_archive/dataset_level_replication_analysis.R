setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(metafor)
source('repos/motrpac/helper_functions.R')
library(org.Hs.eg.db)
entrez2symbol = as.list(org.Hs.egSYMBOL)

# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata

acute_datasets_tests = lapply(acute_datasets,function(x)x$time2ttest_stats)
longterm_datasets_tests = lapply(longterm_datasets,function(x)x$time2ttest_stats)
longterm_datasets_tests = longterm_datasets_tests[sapply(longterm_datasets_tests,length)>0]
acute_datasets_tests = acute_datasets_tests[sapply(acute_datasets_tests,length)>0]

get_genes_intersect<-function(l,internal_list_ind = 1){
  x = rownames(l[[1]][[internal_list_ind]])
  for(j in 2:length(l)){
    x = intersect(x,rownames(l[[j]][[internal_list_ind]]))
  }
  return(x)
}

all_cohorts = union(names(acute_datasets_tests),names(longterm_datasets_tests))
is_acute = grepl("GE_A_",all_cohorts);names(is_acute)=all_cohorts
tissues = c("blood","muscle")
trainings = c("endurance","resistance","all")
types = c("acute","longterm")
rep_analysis_datasets = list()
for(type in types){
  datasets = longterm_datasets_tests;metadata = longterm_metadata;cohorts = all_cohorts[!is_acute]
  if(type=="acute"){
    datasets = acute_datasets_tests;metadata = acute_metadata;cohorts = all_cohorts[is_acute]
  }
  curr_tissues = sapply(metadata[cohorts],function(x)x$tissue)
  curr_trainings = sapply(metadata[cohorts],function(x)x$training)
  table(curr_tissues,curr_trainings)
  for(tissue in tissues){
    for(training in trainings){
      curr_cohorts = cohorts[curr_tissues == tissue&(curr_trainings==training | curr_trainings=="both")]
      if(training=="all"){
        curr_cohorts = cohorts[curr_tissues == tissue&(curr_trainings=="endurance" | curr_trainings=="resistance" | curr_trainings=="both")]
      }
      else{
        #next
      }
      if(length(curr_cohorts)<3){next}
      curr_genes = get_genes_intersect(datasets[curr_cohorts])
      m = c()
      curr_times = c()
      for(cc in curr_cohorts){
        dd = datasets[[cc]]
        for (tt in names(dd)){
          curr_times = c(curr_times,tt)
          m = cbind(m,dd[[tt]][curr_genes,"p"])
          colnames(m)[ncol(m)] = paste(cc,";tp_",tt,sep="")
        }
      }
      currname = paste(type,tissue,training,sep=",")
      if(tissue=="muscle" && type=="acute"){
        m1 = m[,curr_times>=2 & curr_times < 6]
        rep_analysis_datasets[[paste(currname,"2-5h",sep=',')]] = m1
      }
      if(tissue=="blood" && type=="acute"){
        m1 = m[,curr_times < 2]
        rep_analysis_datasets[[paste(currname,"0-1h",sep=',')]] = m1
      }
      rep_analysis_datasets[[currname]] = m
    }
  }
}
sapply(rep_analysis_datasets,dim)
sapply(rep_analysis_datasets,function(x)sum(is.na(x)))
rep_analysis_datasets = lapply(rep_analysis_datasets,function(x){x[is.na(x)]=0.5;return(x)})
sapply(rep_analysis_datasets,dim)
save(rep_analysis_datasets,file="PADB_dataset_level_replicability_analysis_data.RData")

library(kernlab);library(corrplot)
source('~/Desktop/screen/supplementary_data/submission_code/SCREEN_code_for_submission.R')
source('~/Desktop/screen/supplementary_data/submission_code/twogroups_methods_for_submission.R')
extract_study_pairwise_correlations<-function(pvals,lfdr_method='znormix',use_power=T,threegroups=T,B=50,...){
  print ("Analyzing each study")
  mar_est = get_study_marginal_estimation(pvals,lfdr_method=lfdr_method,use_power=use_power,threegroups=threegroups,...)
  lfdr_objs = mar_est$lfdr_objs
  f_1_mat = mar_est$f_1_mat
  f_0_mat = mar_est$f_0_mat
  print ("Done")
  corrs = get_study_pair_corr_matrix(f_1_mat,f_0_mat,B=B,convergenceEps=1e-4)
  colnames(corrs) = colnames(pvals)
  rownames(corrs) = colnames(corrs)
  return(corrs)
}
#corrplot(extract_study_pairwise_correlations(rep_analysis_datasets[[3]],B=10),order="hclust")
screen_res = lapply(rep_analysis_datasets,function(x)SCREEN(x,ks=2:ncol(x),nH=10000))
corr_mats = lapply(rep_analysis_datasets,extract_study_pairwise_correlations)
for(i in 1:length(rep_analysis_datasets)){
  rownames(screen_res[[i]]) = rownames(rep_analysis_datasets[[i]])
}

# Currently we ran screen_res on non "all" datasets
save(rep_analysis_datasets,screen_res,corr_mats,file="PADB_dataset_level_replicability_analysis_data.RData")
sapply(screen_res,function(x)colSums(x<=0.2))
sapply(screen_res,function(x)x["10891",])

rep_num_genes_per_percent = list()
rep_lfdr = 0.2
rep_gene_sets = list(); decay_plot_data = list()
par(mfrow=c(3,4))
for(nn in names(screen_res)){
  currgenes = rownames(screen_res[[nn]])
  mat = screen_res[[nn]]
  num_genes = colSums(mat<=rep_lfdr)
  percents = (2:(1+ncol(mat)))/(1+ncol(mat))
  plot(x=percents,y=num_genes,type='b',main=nn,ylim = c(0,500))
  ind = which(percents >= 0.5)[1]
  selected_genes = currgenes[mat[,ind]<=rep_lfdr]
  if(length(selected_genes)>0){rep_gene_sets[[nn]] = selected_genes}
  decay_plot_data[[nn]] = list()
  decay_plot_data[[nn]][["percents"]] = percents
  decay_plot_data[[nn]][["num_genes"]] = num_genes
}

# Select genes, look at known genes and enrichments
sapply(rep_gene_sets,length)
rep_gene_sets_names = lapply(rep_gene_sets,function(x,y)sort(unlist(y[x])),y=entrez2symbol)
known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
                "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")
sapply(rep_gene_sets_names,intersect,y=known_genes)

all_genes = unique(unlist(sapply(screen_res,rownames)))
rep_gene_sets_topgo = run_topgo_enrichment_fisher(rep_gene_sets,all_genes)
table(extract_top_go_results(rep_gene_sets_topgo)[,1])
get_most_sig_enrichments_by_groups(extract_top_go_results(rep_gene_sets_topgo,0.1),4)

# A similar analysis using a higher lfdr analysis
rep_lfdr = 0.4
rep_gene_sets_0.4 = list()
for(nn in names(screen_res)){
  currgenes = rownames(screen_res[[nn]])
  mat = screen_res[[nn]]
  num_genes = colSums(mat<=rep_lfdr)
  percents = (2:(1+ncol(mat)))/(1+ncol(mat))
  ind = which(percents >= 0.5)[1]
  selected_genes = currgenes[mat[,ind]<=rep_lfdr]
  if(length(selected_genes)>0){rep_gene_sets_0.4[[nn]] = selected_genes}
}
sapply(rep_gene_sets_0.4,length)
rep_gene_sets_0.4_names = lapply(rep_gene_sets_0.4,function(x,y)sort(unlist(y[x])),y=entrez2symbol)
sapply(rep_gene_sets_0.4_names,intersect,y=known_genes)
rep_gene_sets_0.4_topgo = run_topgo_enrichment_fisher(rep_gene_sets_0.4,all_genes)
table(extract_top_go_results(rep_gene_sets_0.4_topgo)[,1])
get_most_sig_enrichments_by_groups(extract_top_go_results(rep_gene_sets_0.4_topgo,0.1),4)

library(Vennerable)
names(rep_gene_sets_0.4_names)
V = Venn(rep_gene_sets_0.4_names[c(3,5,7,8)])
plot(V,doWeights=F)

save(screen_res,decay_plot_data,rep_gene_sets,rep_gene_sets_topgo,rep_gene_sets_names,
     rep_gene_sets_0.4,rep_gene_sets_0.4_names,rep_gene_sets_0.4_topgo,
     file="PADB_dataset_level_replicability_analysis_results.RData")






