# In this script we load the datasets and perform meta-analysis using the metafor package
# This is done both for the acute and longterm datasets
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
sapply(acute_datasets_effects,length)
sapply(acute_datasets_effects,names)
sapply(acute_datasets_effects,function(x)sapply(x,function(y)sum(is.na(y)|is.nan(y))))
sapply(longterm_datasets_effects,function(x)sapply(x,function(y)sum(is.na(y)|is.nan(y))))

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

# extract the data objects for the meta-analysis
acute_genes = rownames(acute_datasets[[1]]$gene_data)
for(i in 2:length(acute_datasets)){
  if(length(acute_datasets[[i]])<3){next}
  acute_genes = intersect(acute_genes,rownames(acute_datasets[[i]]$gene_data))
}
acute_mod = get_dataset_moderators(acute_metadata)
acute_gene_tables = lapply(acute_genes,get_gene_table,dataset_effects=acute_datasets_effects,moderators=acute_mod)
names(acute_gene_tables) = acute_genes
longterm_genes = rownames(longterm_datasets[[1]]$gene_data)
for(i in 2:length(longterm_datasets)){
  if(length(longterm_datasets[[i]])<3){next}
  longterm_genes = intersect(longterm_genes,rownames(longterm_datasets[[i]]$gene_data))
}
longterm_mod = get_dataset_moderators(longterm_metadata)
longterm_gene_tables = lapply(longterm_genes,get_gene_table,dataset_effects=longterm_datasets_effects,moderators=longterm_mod)
names(longterm_gene_tables) = longterm_genes

# remove controls, untrained, and fat from the tables
remove_undesired_datasets<-function(gdata){
  gdata = gdata[!grepl("control",gdata$training),]
  gdata = gdata[!grepl("other",gdata$training),]
  gdata = gdata[!grepl("yoga",gdata$training),]
  gdata = gdata[!grepl("untr",gdata$training),]
  gdata = gdata[!grepl("fat",gdata$tissue),]
  gdata = gdata[!grepl("adipose",gdata$tissue),]
  return (gdata)
}
acute_gene_tables_raw = acute_gene_tables
acute_gene_tables = lapply(acute_gene_tables,remove_undesired_datasets)
longterm_gene_tables_raw = longterm_gene_tables
longterm_gene_tables = lapply(longterm_gene_tables,remove_undesired_datasets)

# # transform the datasets into matrices that can be plotted
# acute_effects_matrix = c();acute_sds_matrix = c()
# for(d in names(acute_datasets_effects)){
#   l = acute_datasets_effects[[d]]
#   times = names(l)
#   darr = strsplit(split=';',d)[[1]]
#   darr[2] = simplify_tissue_info(darr[2])
#   for (tt in times){
#     currname = paste(c(tt,darr[c(2,4,1)]),collapse=';')
#     acute_effects_matrix = cbind(acute_effects_matrix,l[[tt]][acute_genes,1])
#     acute_sds_matrix = cbind(acute_sds_matrix,l[[tt]][acute_genes,2])
#     colnames(acute_sds_matrix)[ncol(acute_sds_matrix)] = currname
#   }
# }
# colnames(acute_effects_matrix) = colnames(acute_sds_matrix)
# longterm_effects_matrix = c();longterm_sds_matrix = c()
# for(d in names(longterm_datasets_effects)){
#   l = longterm_datasets_effects[[d]]
#   times = names(l)
#   darr = strsplit(split=';',d)[[1]]
#   darr[2] = simplify_tissue_info(darr[2])
#   for (tt in times){
#     currname = paste(c(tt,darr[c(2,4,1)]),collapse=';')
#     longterm_effects_matrix = cbind(longterm_effects_matrix,l[[tt]][longterm_genes,1])
#     longterm_sds_matrix = cbind(longterm_sds_matrix,l[[tt]][longterm_genes,2])
#     colnames(longterm_sds_matrix)[ncol(longterm_sds_matrix)] = currname
#   }
# }
# colnames(longterm_effects_matrix) = colnames(longterm_sds_matrix)

save(acute_gene_tables_raw,acute_gene_tables,longterm_gene_tables_raw,longterm_gene_tables,
     file="PADB_dataset_level_meta_analysis_data.RData")

# ... e.g., subset = (tissue=="blood")
get_gene_analysis_pvals<-function(gdata,use_mods=T,func=rma,...){
  gdata = gdata[gdata$vi>0,]
  v = NULL
  try({
    if(use_mods){
      res1 = func(yi=yi,vi=vi, mods = ~ training + time ,data=gdata, control=list(maxiter=10000,stepadj=0.5),...) 
    }
    else{
      res1 = func(yi=yi,vi=vi,data=gdata, control=list(maxiter=10000,stepadj=0.5),...) 
    }
    ps1 = res1$pval;names(ps1) = rownames(res1$b)
    ps1 = c(ps1,res1$QMp);names(ps1)[length(ps1)] = "AllMods"
    stats1 = res1$b[,1];ses1 = res1$se
    names(ps1) = paste("pval_",names(ps1),sep="")
    names(stats1) = paste("stat_",names(stats1),sep="")
    v = c(ps1,stats1)
  })
  return(v)
}
get_gene_analysis_pvals_with_gse_correction<-function(gdata,use_mods=T,func=rma.mv,...){
  gdata = gdata[gdata$vi>0,]
  v = NULL
  try({
    if(use_mods){
      res1 = func(yi,vi,mods = ~ training + time ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...) 
    }
    else{
      res1 = func(yi,vi,data=gdata, random = ~ 1|gse,control=list(maxiter=10000),...)
    }
    ps1 = res1$pval;names(ps1) = rownames(res1$b)
    ps1 = c(ps1,res1$QMp);names(ps1)[length(ps1)] = "AllMods"
    stats1 = res1$b[,1];ses1 = res1$se
    names(ps1) = paste("pval_",names(ps1),sep="")
    names(stats1) = paste("stat_",names(stats1),sep="")
    v = c(ps1,stats1)
  })
  return(v)
}

# save results as lists: sometimes the solvers fail
acute_meta_analysis_results = list()
acute_meta_analysis_results[["random_effects_blood"]] = lapply(acute_gene_tables,
                function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="blood",],use_mods=F))
acute_meta_analysis_results[["random_effects_muscle"]] = lapply(acute_gene_tables,
                function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="muscle",],use_mods=F))
acute_meta_analysis_results[["mixed_effects_blood"]] = lapply(acute_gene_tables,
                function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="blood",],use_mods=T))
acute_meta_analysis_results[["mixed_effects_muscle"]] = lapply(acute_gene_tables,
                function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="muscle",],use_mods=T))

longterm_meta_analysis_results = list()
longterm_meta_analysis_results[["random_effects_blood"]] = lapply(longterm_gene_tables,
        function(x)get_gene_analysis_pvals(x[x$tissue=="blood",],use_mods=F))
longterm_meta_analysis_results[["random_effects_muscle"]] = lapply(longterm_gene_tables,
        function(x)get_gene_analysis_pvals(x[x$tissue=="muscle",],use_mods=F))
longterm_meta_analysis_results[["mixed_effects_blood"]] = lapply(longterm_gene_tables,
        function(x)get_gene_analysis_pvals(x[x$tissue=="blood",],use_mods=T))
longterm_meta_analysis_results[["mixed_effects_muscle"]] = lapply(longterm_gene_tables,
        function(x)get_gene_analysis_pvals(x[x$tissue=="muscle",],use_mods=T))

longterm_meta_analysis_results[["random_effects_muscle_with_gse"]] = lapply(longterm_gene_tables, 
      function(x)get_gene_analysis_pvals_with_gse_correction(
      x[x$tissue=="muscle",],use_mods=F))

save(acute_meta_analysis_results,longterm_meta_analysis_results,file="PADB_metafor_meta_analysis_results.RData")

load("PADB_metafor_meta_analysis_results.RData")
get_ps<-function(x,ind=1){
  nn = names(x)
  xl = lapply(x,function(x)x[ind])
  xl[sapply(xl,length)==0] = NA
  v = unlist(xl)
  names(v) = nn
  return (v)
}
library(locfdr)
get_lfdrs<-function(pp){
  pp[is.na(pp)]=0.5
  if(sd(pp)==0){return(matrix(1,nrow=length(pp),ncol=3))}
  pp[pp==0] = min(pp[pp>0])
  pp[pp==1] = max(pp[pp<1])
  zz = qnorm(pp)
  lfdr = locfdr(zz)
  return(cbind(pp,lfdr$fdr,zz))
}
acute_ps_intersect = sapply(acute_meta_analysis_results,get_ps)
acute_ps_time = sapply(acute_meta_analysis_results[3:4],get_ps,ind="pval_time")
acute_ps_tr = sapply(acute_meta_analysis_results[3:4],get_ps,ind="pval_trainingresistance")
acute_ps = cbind(acute_ps_intersect,acute_ps_time,acute_ps_tr)
colnames(acute_ps)[5:6] = c("blood, time","muscle,time")
colnames(acute_ps)[7:8] = c("blood, training","muscle,training")
longterm_ps_intersect = sapply(longterm_meta_analysis_results[-5],get_ps)
longterm_ps_time = sapply(longterm_meta_analysis_results[3:4],get_ps,ind="pval_time")
longterm_ps_tr = sapply(longterm_meta_analysis_results[3:4],get_ps,ind="pval_trainingresistance")
longterm_ps = cbind(longterm_ps_intersect,longterm_ps_time,longterm_ps_tr)
colnames(longterm_ps)[5:6] = c("blood, time","muscle,time")
colnames(longterm_ps)[7:8] = c("blood, training","muscle,training")
acute_qs = apply(acute_ps,2,p.adjust,method='fdr')
longterm_qs = apply(longterm_ps,2,p.adjust,method='fdr')
acute_lfdrs = apply(acute_ps,2,function(x)get_lfdrs(x)[,2])
longterm_lfdrs = apply(longterm_ps,2,function(x)get_lfdrs(x)[,2])
longterm_zzs = apply(longterm_ps,2,function(x)get_lfdrs(x)[,3])
acute_zzs = apply(acute_ps,2,function(x)get_lfdrs(x)[,3])

# par(mfrow=c(2,2))
# hist(acute_ps[,1],main="Acute, blood")
# hist(acute_ps[,2],main="Acute, muscle")
# hist(longterm_ps[,1],main="Longterm, blood")
# hist(longterm_ps[,2],main="Longterm, muscle")

metafor_gene_sets = list()
metafor_gene_sets[["0.1fdr"]] = c(
  apply(acute_qs,2,function(x,y)y[x<=0.1 & !is.na(x)],y=rownames(acute_ps)),
  apply(longterm_qs,2,function(x,y)y[x<=0.1 & !is.na(x)],y=rownames(longterm_ps))
)
names(metafor_gene_sets[["0.1fdr"]]) = paste("longterm",names(metafor_gene_sets[["0.1fdr"]]),sep=",")
names(metafor_gene_sets[["0.1fdr"]])[1:ncol(acute_ps)] = gsub(pattern="longterm,",
                  names(metafor_gene_sets[["0.1fdr"]])[1:ncol(acute_ps)],replace="acute,")

metafor_gene_sets[["0.1lfdr"]] = list()
for(j in 1:ncol(acute_ps)){
  nn = paste("acute",colnames(acute_ps)[j],sep=',')
  metafor_gene_sets[["0.1lfdr"]][[nn]] = rownames(acute_ps)[acute_lfdrs[,j]<=0.1 & acute_zzs[,j]<(-1)]
}
for(j in 1:ncol(longterm_ps)){
  nn = paste("longterm",colnames(longterm_ps)[j],sep=',')
  metafor_gene_sets[["0.1lfdr"]][[nn]] = rownames(longterm_ps)[longterm_lfdrs[,j]<=0.1 & longterm_zzs[,j]<(-1)]
}

metafor_gene_sets = lapply(metafor_gene_sets,function(x)x[sapply(x,length)>0])
metafor_sets_enrichments = run_topgo_enrichment_fisher(metafor_gene_sets[[2]],rownames(acute_ps))
extract_top_go_results(metafor_sets_enrichments)
save(metafor_gene_sets,metafor_sets_enrichments,acute_ps,longterm_ps,file="metafor_gene_sets.RData")

# Tests and comments from the paper of metafor (2010)
gdata = acute_gene_tables_raw[["10891"]] # PGC1 in acute response
gdata = longterm_gene_tables[["1282"]] # COL1 gene
gdata = longterm_gene_tables[["4168"]] # A negative example
gdata = acute_gene_tables[["5166"]] # survived rep but not meta
# gdata = acute_gene_tables[["11326"]] # gene with significant modifiers
# another test: make a dataset with large effects
# gdata$yi = rnorm(nrow(gdata),mean=5)
plot(gdata[,"yi"],gdata[,"vi"],pch=as.numeric(as.factor(gdata$tissue)))
# knha - a correction that accounts for the uncertainty in the random effect
res0 = rma.mv(yi,vi,random = ~ 1|gse,data=gdata,subset = (tissue=="muscle"))
# explanation of the result above:
#   mu - the average effect is 0.067, the CI contains zero
res = rma.mv(yi,vi,mods = ~  training + time, random= ~1|gse ,data=gdata[gdata$tissue=="muscle",])
# CIs of the anova stats:
confint(res0)
# Forest plot - very informative
forest(res)
# difference in tau before and after using moderators - teaches us about the
# percentage of explained variance due to the moderators
# when the test for residuals (QE) is signficant - we may be missing additional
# moderators
# We can use predict to get expected effects for new moderators:
# Currently does not work because factors should be transformed into dummy vars.
# Also pages 18-19 show nice figures and analysis of the predictions.
predict(res,newmods = as.matrix(data.frame("",time=5,"",training="endurance",tissue="muscle")), addx = TRUE)
predict(res,newmods = as.matrix(gdata))
# look at the fitted values
predict(res)
# Separate by tissue
res2 = rma(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, subset = (tissue=="muscle"))
res1 = rma(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, subset = (tissue=="blood"))
# Residual analysis for detecting outlier datasets
barplot(as.numeric(rstudent(res)$z));abline(-2,0);abline(2,0)
# Residual analysis is informative but not enough
# case deletion diagnostics are informatice as well
plot(influence(res))
# funnel plots
funnel(res0)
# radial plots: useful for consistency analysis
# can be used only for models without moderators
radial(res0)
# qq plots for the standardized residuals
qqnorm(res0,main="random")
qqnorm(res,main="mixed, both tissues")
qqnorm(res1,main="blood")
qqnorm(res2,main="muscle")
# tests for publication bias
regtest(res0,predictor="vi",model="lm")
regtest(res2,predictor="vi",model="lm")
# anova tests
anova(res0,res)

# # Compare to the rmeta package
# install.packages('rmeta')
# library(rmeta)
# blood_gdata = gdata[gdata$tissue=="blood",]
# rmeta_res0 = meta.summaries(d=blood_gdata$yi,se = blood_gdata$vi)
# summary(rmeta_res0)[[3]]
# plot(rmeta_res0)
# funnelplot(rmeta_res0)
# gdata$yi
# 
# # Tests on simulated data
# n = 5;effect = 3; hetero=0.5; effect2=6
# vi = rep(0.5,n)
# yi = rnorm(n,sd=vi) + effect + rnorm(n,sd=hetero)
# rma(yi,vi)$pval
# dummy = 1:n
# rma.mv(yi,vi,random=~ 1|dummy)$pval
# # vs.
# rma(c(yi,rnorm(n,sd=vi) + rnorm(n,sd=hetero)),c(vi,vi))$pval
# # vs.
# ref = factor(c(rep(1,n),rep(2,n)))
# rma.mv(c(yi,rnorm(n,sd=vi) + rnorm(n,sd=hetero)),c(vi,vi),random = ~1|ref)$pval
# 
# # random effects
# beta1 = 1 ; beta2 = 0.5
# vi = c(vi,vi)
# yi = c(yi+beta1,yi+beta2)
# ref = factor(c(rep(1,n),rep(2,n)))
# rma(yi,vi)
# rma.mv(yi,vi,rand = ~1|ref)










