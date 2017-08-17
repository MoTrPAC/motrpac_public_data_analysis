# In this script we load the datasets and perform meta-analysis using the metafor package
# This is done both for the acute and longterm datasets
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(metafor)
source('repos/motrpac/helper_functions.R')

# Get the datasets and their metadata
acute_datasets = get(load("PADB_univariate_results_and_preprocessed_data_acute.RData"))
acute_metadata = get(load("PADB_sample_metadata_acute.RData"))
longterm_datasets = get(load("PADB_univariate_results_and_preprocessed_data_longterm.RData"))
longterm_metadata = get(load("PADB_sample_metadata_longterm.RData"))

# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_yi_vi <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  d = x2-x1
  sdd = sd(d)/sqrt(length(x1))
  return(c(yi=mean(d),vi=sdd))
}

get_ttest_yi_vi_per_dataset<-function(mat,metadata){
  dataset_times = metadata$time[colnames(mat)]
  if(any(is.na(dataset_times))){
    mat = mat[,!is.na(dataset_times)]
    dataset_times = metadata$time[colnames(mat)]
  }
  min_time = min(dataset_times)
  other_times = setdiff(unique(dataset_times),min_time)
  times2effects = list()
  for(other_time in other_times){
    curr_mat = mat[,dataset_times==min_time | dataset_times==other_time]
    curr_subjects = metadata$subject[colnames(curr_mat)]
    subjects_to_keep = names(which(table(curr_subjects)==2))
    curr_mat = curr_mat[,is.element(curr_subjects,set = subjects_to_keep)]
    ord = order(metadata$subject[colnames(curr_mat)],metadata$time[colnames(curr_mat)])
    curr_mat = curr_mat[,ord]
    curr_times = metadata$time[colnames(curr_mat)]
    paired_test_data = t(apply(curr_mat,1,get_paired_ttest_yi_vi,sample2time=curr_times,t1=min_time,t2=other_time))
    times2effects[[as.character(other_time)]] = paired_test_data
  }
  return(times2effects)
}

acute_datasets_effects = lapply(acute_datasets,function(x,y)get_ttest_yi_vi_per_dataset(x$gene_data,y),y=acute_metadata)
longterm_datasets_effects = lapply(longterm_datasets,function(x,y)get_ttest_yi_vi_per_dataset(x$gene_data,y),y=longterm_metadata)
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
get_dataset_moderators<-function(datasets,metadata){
  arrs = t(sapply(datasets,function(x)strsplit(x,split=';')[[1]][c(1,2,4)]))
  arrs[,2] = simplify_tissue_info(arrs[,2])
  arrs[grepl(datasets,pattern = 'control',ignore.case = T) | grepl(datasets,pattern="untrain",ignore.case = T),3] = "control"
  colnames(arrs) = c("gse","tissue","training")
  return(arrs)
}
get_gene_table<-function(gene,dataset_effects,moderators){
  gene_data = lapply(dataset_effects,function(x,y)sapply(x,function(u,v)u[v,],v=y),y=gene)
  m = c()
  for(nn in names(gene_data)){
    mm = gene_data[[nn]]
    for(j in 1:ncol(mm)){
      m = rbind(m,c(nn,colnames(mm)[j],moderators[nn,],mm[,j]))
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
  acute_genes = intersect(acute_genes,rownames(acute_datasets[[i]]$gene_data))
}
acute_mod = get_dataset_moderators(names(acute_datasets_effects),acute_metadata)
acute_gene_tables = lapply(acute_genes,get_gene_table,dataset_effects=acute_datasets_effects,moderators=acute_mod)
names(acute_gene_tables) = acute_genes
longterm_genes = rownames(longterm_datasets[[1]]$gene_data)
for(i in 2:length(longterm_datasets)){
  if(length(acute_gene_tables)==0){next}
  longterm_genes = intersect(longterm_genes,rownames(longterm_datasets[[i]]$gene_data))
}
longterm_mod = get_dataset_moderators(names(longterm_datasets_effects),longterm_metadata)
longterm_gene_tables = lapply(longterm_genes,get_gene_table,dataset_effects=longterm_datasets_effects,moderators=longterm_mod)
names(longterm_gene_tables) = longterm_genes

get_gene_analysis_pvals<-function(gdata,use_mods=T,func=rma,...){
  gdata = gdata[gdata$vi>0,]
  v = NULL
  try({
    if(use_mods){
      res1 = func(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, 
                 subset = (tissue=="blood"),control=list(maxiter=10000,stepadj=0.5),...) 
      res2 = func(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, 
                 subset = (tissue=="muscle"),control=list(maxiter=10000,stepadj=0.5),...)
    }
    else{
      res1 = func(yi=yi,vi=vi,data=gdata,knha=T, subset = (tissue=="blood"),control=list(maxiter=10000,stepadj=0.5),...) 
      res2 = func(yi=yi,vi=vi,data=gdata,knha=T, subset = (tissue=="muscle"),control=list(maxiter=10000,stepadj=0.5),...)
    }
    ps1 = res1$pval;names(ps1) = rownames(res1$b)
    ps1 = c(ps1,res1$QMp);names(ps1)[length(ps1)] = "AllMods"
    stats1 = res1$b[,1];ses1 = res1$se
    ps2 = res2$pval;names(ps2) = rownames(res2$b)
    stats2 = res2$b[,1];ses2 = res2$se
    ps2 = c(ps2,res2$QMp);names(ps2)[length(ps2)] = "AllMods"
    names(ps1) = paste("blood_p_",names(ps1),sep="")
    names(ps2) = paste("muscle_p_",names(ps2),sep="")
    names(stats1) = paste("blood_s_",names(stats1),sep="")
    names(stats2) = paste("muscle_s_",names(stats2),sep="")
    v = c(ps1,ps2,stats1,stats2)
  })
  return(v)
}
acute_meta_analysis_res = t(sapply(acute_gene_tables,get_gene_analysis_pvals,use_mods=F))
longterm_meta_analysis_res = t(sapply(longterm_gene_tables,get_gene_analysis_pvals,use_mods=F))
# longterm had several fails, the above line returns a list
names(longterm_meta_analysis_res) = names(longterm_gene_tables)
longterm_meta_analysis_res = longterm_meta_analysis_res[!sapply(longterm_meta_analysis_res,is.null)]
tmp_mat = t(sapply(longterm_meta_analysis_res,function(x) x ))
longterm_meta_analysis_res = tmp_mat
save(acute_meta_analysis_res,longterm_meta_analysis_res,file="PADB_metafor_meta_analysis_results.RData")

par(mfrow=c(2,2))
hist(acute_meta_analysis_res[,"blood_p_intrcpt"])
hist(longterm_meta_analysis_res[,"blood_p_intrcpt"])
hist(acute_meta_analysis_res[,"muscle_p_intrcpt"])
hist(longterm_meta_analysis_res[,"muscle_p_intrcpt"])

metafor_gene_sets = list()
metafor_gene_sets[["acute_muscle"]] = names(which(p.adjust(acute_meta_analysis_res[,"muscle_p_intrcpt"],method='fdr')<0.1))
metafor_gene_sets[["acute_blood"]] = names(which(p.adjust(acute_meta_analysis_res[,"blood_p_intrcpt"],method='fdr')<0.1))
metafor_gene_sets[["longterm_muscle"]] = names(which(p.adjust(longterm_meta_analysis_res[,"muscle_p_intrcpt"],method='fdr')<0.1))
metafor_gene_sets[["longterm_blood"]] = names(which(p.adjust(longterm_meta_analysis_res[,"blood_p_intrcpt"],method='fdr')<0.1))
metafor_gene_sets_enrichments = run_topgo_enrichment_fisher(metafor_gene_sets,union(acute_genes,longterm_genes))

# tests and comments from the paper of metafor (2010)
gdata = get_gene_table("8863",acute_datasets_effects,acute_mod)
gdata = get_gene_table("65121",longterm_datasets_effects,longterm_mod)
# another test: make a dataset with large effects
# gdata$yi = rnorm(nrow(gdata),mean=5)
plot(gdata[,"yi"],gdata[,"vi"],pch=as.numeric(as.factor(gdata$tissue)))
# knha - a correction that accounts for the uncertainty in the random effect
res0 = rma(yi,vi,data=gdata,knha=T)
res0 = rma.mv(yi,vi,random = ~ 1|gse,data=gdata)
res0 = rma(yi,vi,data=gdata,knha=T,subset = (tissue=="muscle"))
res0 = rma.mv(yi,vi,random = ~ 1|gse,data=gdata,subset = (tissue=="muscle"))
# explanation of the result above:
#   mu - the average effect is 0.067, the CI contains zero
res = rma(yi=yi,vi=vi,mods = ~ tissue + training + time ,data=gdata,knha=T)
# CIs of the anova stats:
confint(res)
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
funnel(res)
# radial plots: useful for consistency analysis
# can be used only for models without moderators
radial(res0)
# qq plots for the standardized residuals
par(mfrow=c(2,2))
qqnorm(res0,main="random")
qqnorm(res,main="mixed, both tissues")
qqnorm(res1,main="blood")
qqnorm(res2,main="muscle")
# tests for publication bias
regtest(res0,predictor="vi",model="lm")
regtest(res2,predictor="vi",model="lm")
# anova tests
anova(res0,res)




