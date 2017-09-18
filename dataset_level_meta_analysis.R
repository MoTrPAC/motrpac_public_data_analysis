# In this script we load the datasets and perform meta-analysis using the metafor package
# This is done both for the acute and longterm datasets
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(metafor)
library(org.Hs.eg.db)
source('repos/motrpac/helper_functions.R')

# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata

# Some basic statistics
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
get_basic_stats_from_metadata(acute_metadata[!acute_cohorts])
longterm_cohorts = sapply(longterm_datasets,function(x)!is.null(x$gene_fchanges))
get_basic_stats_from_metadata(longterm_metadata[longterm_cohorts])
get_basic_stats_from_metadata(longterm_metadata[!longterm_cohorts])

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

save(acute_gene_tables_raw,acute_gene_tables,longterm_gene_tables_raw,longterm_gene_tables,
     file="PADB_dataset_level_meta_analysis_data.RData")

########################################################################
#######################  Meta-analysis functions  ######################
########################################################################
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
      res1 = get_rma_obj_with_mods(gdata,func)
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

# If training has controls: add a new binary vector for controls
get_rma_obj_with_mods<-function(gdata,func=rma.mv){
  has_times = length(unique(gdata$time))>1
  has_tr = length(unique(gdata$training))>1
  if(has_tr && has_times){
    res1 = func(yi,vi,mods = ~ training + time ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...) 
  }
  if(!has_tr && !has_times){
    res1 = func(yi,vi,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...) 
  }
 if(has_times && ! has_tr){
    res1 = func(yi,vi,mods = ~ time ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...)
  }
 if(has_tr && ! has_times){
    res1 = func(yi,vi,mods = ~ training ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...)
  }
  return(res1)
}

perform_control_test<-function(gdata,func=rma.mv,...){
  try({
  has_times = length(unique(gdata$time))>1
  gdata = gdata[gdata$training!="other",]
  is_ctrl = factor(gdata$training=="control")
  gdata$is_ctrl=is_ctrl
  if(has_times){
    res1 = func(yi,vi,mods = ~ is_ctrl + training + time ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...) 
  }
  if(!has_times){
    res1 = func(yi,vi,mods = ~ is_ctrl + training ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000),...)
  }
  ind = grepl("is_ctrl",rownames(res1$beta))
  return(res1$pval[ind])
  })
  return(NA)
}

get_ps<-function(x,ind=1){
  nn = names(x)
  xl = lapply(x,function(x)x[ind])
  xl[sapply(xl,length)==0] = NA
  v = unlist(xl)
  names(v) = nn
  return (v)
}
library(locfdr)
get_lfdrs<-function(pp,...){
  pp[is.na(pp)]=0.5
  if(sd(pp)==0){return(matrix(1,nrow=length(pp),ncol=3))}
  pp[pp==0] = min(pp[pp>0])
  pp[pp==1] = max(pp[pp<1])
  zz = qnorm(pp)
  lfdr = locfdr(zz,...)
  return(cbind(pp,lfdr$fdr,zz))
}
library(lme4)
# For egger and correlations: take the first time point from each gse
# Use mixed effects with all to measure bias
perform_bias_analysis<-function(gdata,c2n,tissue="muscle"){
  gdata = gdata[gdata$vi>0,]
  gdata = gdata[gdata$tissue==tissue,]
  d = data.frame(y=gdata$yi,a=gdata$vi,b=as.numeric(c2n[gdata$V1]),z=as.character(gdata$gse))
  # lme1 = lmer(formula = y ~ a+b+(1|z),d,REML=F)
  # lme0 = lmer(formula = y ~ (1|z),d,REML=F)
  # av_p = anova(lme0,lme1)[2,8]
  lm1 = lm(y~a+z,d)
  av_p = summary(lm1)[[4]][2,4]
  gdata1 = apply(gdata,2,function(x,y)tapply(x,y,function(x)x[1]),y=gdata$V1)
  gdata = data.frame(gdata1)
  gdata$yi = as.numeric(as.character(gdata$yi))
  gdata$vi = as.numeric(as.character(gdata$vi))
  eg_test = NA
  try({
    res = rma(yi=yi,vi=vi,data=gdata, control=list(maxiter=10000,stepadj=0.5)) 
    eg_test = regtest(res)$pval
  })
  size_effect_test = cor.test(gdata$yi,c2n[gdata$V1],method="spearman")
  size_effect_cor = size_effect_test$estimate
  size_effect_p = size_effect_test$p.value
  size_precision_cor = cor(gdata$vi,c2n[gdata$V1],method="spearman")
  return(c(eggerp = eg_test,size_effect_cor=size_effect_cor,
           size_effect_p=size_effect_p,size_precision_cor=size_precision_cor,
           me_with_vi_and_size=av_p))
}
#perform_bias_analysis(acute_gene_tables[["10891"]],acute_c2n)
#perform_bias_analysis(longterm_gene_tables[["10891"]],longterm_c2n)

########################################################################
########################################################################
########################################################################

load("PADB_dataset_level_meta_analysis_data.RData")

################################
# Bug fix: Sept 11 2017: the objects have standard errors in the vi
# we use the square for correction
correct_vi<-function(gdata){gdata$vi = gdata$vi^2;return(gdata)}
acute_gene_tables = lapply(acute_gene_tables,correct_vi)
acute_gene_tables_raw = lapply(acute_gene_tables_raw,correct_vi)
longterm_gene_tables = lapply(longterm_gene_tables,correct_vi)
longterm_gene_tables_raw = lapply(longterm_gene_tables_raw,correct_vi)
#################################

# Estimation of publication bias
publication_bias_res = list()

acute_c2n = sapply(acute_metadata,function(x)length(x$gsms))
c2num_ts = table(acute_gene_tables_raw[[1]]$V1)
acute_c2n[names(c2num_ts)] = acute_c2n[names(c2num_ts)]/c2num_ts
publication_bias_res[["acute,muscle"]] = sapply(acute_gene_tables_raw,perform_bias_analysis,c2n=acute_c2n)
publication_bias_res[["acute,blood"]] = sapply(acute_gene_tables_raw,perform_bias_analysis,c2n=acute_c2n,tissue="blood")

longterm_c2n = sapply(longterm_metadata,function(x)length(x$gsms))
c2num_ts = table(longterm_gene_tables_raw[[1]]$V1)
longterm_c2n[names(c2num_ts)] = longterm_c2n[names(c2num_ts)]/c2num_ts
publication_bias_res[["longterm,muscle"]] = sapply(longterm_gene_tables_raw,perform_bias_analysis,c2n=longterm_c2n)
publication_bias_res[["longterm,blood"]] = sapply(longterm_gene_tables_raw,perform_bias_analysis,c2n=longterm_c2n,tissue="blood")

par(mfrow=c(2,2))
hist(publication_bias_res[["acute,muscle"]][1,],main="Acute, Egger test", xlab="p-value")
hist(publication_bias_res[["longterm,muscle"]][1,],main="Longterm, Egger test", xlab="p-value")
hist(publication_bias_res[["acute,blood"]][1,],main="Acute, Egger test", xlab="p-value")
hist(publication_bias_res[["longterm,blood"]][1,],main="Longterm, Egger test", xlab="p-value")
# sapply(publication_bias_res,function(x)table(x[1,]<0.01))
# 
# gene = "10891"
# gdata = acute_gene_tables[[gene]]
# res = rma(yi,vi,data=gdata)
# funnel(res)
# regtest(res)
# perform_bias_analysis(gdata,acute_c2n)
# perform_bias_analysis(gdata,acute_c2n,"blood")
# 
# all_genes = unique(unlist(metafor_gene_sets[[2]]))
# gss = lapply(publication_bias_res,function(x,y)intersect(y,colnames(x)),y=all_genes)
# hist(publication_bias_res[["acute"]][5,gss[[1]]])
# hist(publication_bias_res[["longterm"]][5,gss[[2]]])

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

save(acute_meta_analysis_results,longterm_meta_analysis_results,publication_bias_res,file="PADB_metafor_meta_analysis_results.RData")

################################
# New analysis: Sept 12 2017: merge time points, used as ordered factors
# we merge time points such that there are at least 3 datasets per TP
# Additional analyses:
# Analyze controls vs. non controls
# Check model with controls
# Shave off conrol and/or publication bias genes
simplify_time_in_gdata<-function(gdata,func=simplify_time_acute){
  gdata$time = func(gdata$time)
  gdata$time = ordered(gdata$time)
  return(gdata)
}
acute_gene_tables_simpletime = lapply(acute_gene_tables,simplify_time_in_gdata)
acute_gene_tables_raw_simpletime = lapply(acute_gene_tables_raw,simplify_time_in_gdata)
longterm_gene_tables_simpletime = lapply(longterm_gene_tables,simplify_time_in_gdata,func=simplify_time_longterm)
longterm_gene_tables_raw_simpletime = lapply(longterm_gene_tables_raw,simplify_time_in_gdata,func=simplify_time_longterm)

meta_analysis_results = list()
meta_analysis_results[["acute,muscle"]] = lapply(acute_gene_tables_simpletime,
    function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="muscle",],use_mods=T))
meta_analysis_results[["acute,blood"]] = lapply(acute_gene_tables_simpletime,
    function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="blood",],use_mods=T))
meta_analysis_results[["longterm,muscle"]] = lapply(longterm_gene_tables_simpletime,
    function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="muscle",],use_mods=T))
meta_analysis_results[["longterm,blood"]] = lapply(longterm_gene_tables_simpletime,
    function(x)get_gene_analysis_pvals_with_gse_correction(x[x$tissue=="blood",],use_mods=T))
meta_analysis_results[["acute,controls,muscle"]] = sapply(acute_gene_tables_raw_simpletime,
    function(x)perform_control_test(x[x$tissue=="muscle",]))
meta_analysis_results[["acute,controls,blood"]] = sapply(acute_gene_tables_raw_simpletime,
    function(x)perform_control_test(x[x$tissue=="blood",]))
meta_analysis_results[["longterm,controls,muscle"]] = sapply(longterm_gene_tables_raw_simpletime,
    function(x)perform_control_test(x[x$tissue=="muscle",]))
# # longterm blood does not have controls...
# meta_analysis_results[["longterm,controls,blood"]] = sapply(longterm_gene_tables_raw_simpletime[1:3],
#     function(x)perform_control_test(x[x$tissue=="blood",]))
sapply(meta_analysis_results[5:7],hist)

get_ps_for_simpletime_analysis<-function(x){
  p1 = get_ps(x)
  p2 = get_ps(x,ind="pval_AllMods")
  return(cbind(p1,p2))
}
meta_analysis_pvals = lapply(meta_analysis_results[1:4],get_ps_for_simpletime_analysis)
x1 = meta_analysis_pvals[[1]][,2]
x3 = meta_analysis_pvals[[1]][,1]
x2 = meta_analysis_results$`acute,controls,muscle`
plot(-log(x1),-log(x2))
sort(names(which((x3<0.001 | x1<0.0001)&x2<0.0001)))
table((x3<0.001 | x1<0.0001),x2<0.0001)
x1["10891"]
x2["10891"]
save(meta_analysis_results,file="New_metafor_analysis_sept_12.RData")

all_ps = c(unlist(meta_analysis_interc_pvals),unlist(meta_analysis_mod_pvals))
all_ps = all_ps[!is.na(all_ps)]
pthr = max(all_ps[p.adjust(all_ps,method='BY')<0.01])

sapply(meta_analysis_mod_pvals,function(x)table(x<pthr))
sapply(meta_analysis_interc_pvals,function(x)table(x<pthr))
x = get_lfdrs(meta_analysis_mod_pvals[[4]])
which(x[,2]<0.01)

gene = "634" # TNNT2 - expected but has poor signal
entrez2symbol[gene]
gdata = acute_gene_tables[[gene]] 
get_subset_forest_plot(gdata,"muscle")
get_subset_forest_plot(gdata,"blood")
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = F,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL)
gene_ps = sapply(acute_datasets_effects,function(x)sapply(x,function(y)y[gene,"p"]))
gene_ps[["GE_A_8"]]

gdata = longterm_gene_tables_simpletime[[gene]]
#gdata = acute_gene_tables[[gene]]
gdata = gdata[gdata$tissue=="blood",]
res1 = func(yi,vi,mods = ~ training + time ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000))
res2 = func(yi,vi,data = gdata, random = ~ 1|gse, control=list(maxiter=10000))
summary(res2)[[2]]
publication_bias_res[[1]][,gene]
anova(res1,res3)
summary(lm(yi~time+training+factor(gse),data=gdata))

#################################


##########################################################
##########################################################
##########################################################
##########################################################
# Gene selection and display items
get_matrix_p_adjust<-function(x,q=0.1,...){
  v = c(x)
  v = v[!is.na(x)]
  vq = p.adjust(v,...)
  thr = max(v[vq<=q])
  return(thr)
}

load("PADB_metafor_meta_analysis_results.RData")
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
acute_qs_thr = get_matrix_p_adjust(acute_ps,0.1,method='fdr')
longterm_qs_thr = get_matrix_p_adjust(longterm_ps,0.1,method='fdr')
acute_lfdrs = apply(acute_ps,2,function(x)get_lfdrs(x)[,2])
longterm_lfdrs = apply(longterm_ps,2,function(x)get_lfdrs(x)[,2])
longterm_zzs = apply(longterm_ps,2,function(x)get_lfdrs(x)[,3])
acute_zzs = apply(acute_ps,2,function(x)get_lfdrs(x)[,3])

metafor_gene_sets = list()
metafor_gene_sets[["0.1fdr"]] = c(
  apply(acute_ps,2,function(x,y)y[x<=acute_qs_thr & !is.na(x)],y=rownames(acute_ps)),
  apply(longterm_ps,2,function(x,y)y[x<=longterm_qs_thr & !is.na(x)],y=rownames(longterm_ps))
)
names(metafor_gene_sets[["0.1fdr"]]) = paste("longterm",names(metafor_gene_sets[["0.1fdr"]]),sep=",")
names(metafor_gene_sets[["0.1fdr"]])[1:ncol(acute_ps)] = gsub(pattern="longterm,",
                  names(metafor_gene_sets[["0.1fdr"]])[1:ncol(acute_ps)],replace="acute,")
sapply(metafor_gene_sets[[1]],length)

sapply(metafor_gene_sets[[1]],length)
pb_genes = publication_bias_res$`acute,blood`[1,]<0.01
intersect(metafor_gene_sets[[1]]$`acute,blood, time`,names(which(pb_genes)))

metafor_gene_sets[["0.1lfdr"]] = list()
for(j in 1:ncol(acute_ps)){
  nn = paste("acute",colnames(acute_ps)[j],sep=',')
  metafor_gene_sets[["0.1lfdr"]][[nn]] = rownames(acute_ps)[acute_lfdrs[,j]<=0.1 & acute_zzs[,j]<(-1)]
}
for(j in 1:ncol(longterm_ps)){
  nn = paste("longterm",colnames(longterm_ps)[j],sep=',')
  metafor_gene_sets[["0.1lfdr"]][[nn]] = rownames(longterm_ps)[longterm_lfdrs[,j]<=0.1 & longterm_zzs[,j]<(-1)]
}
sapply(metafor_gene_sets[[2]],length)

metafor_gene_sets = lapply(metafor_gene_sets,function(x)x[sapply(x,length)>0])
metafor_sets_enrichments = run_topgo_enrichment_fisher(metafor_gene_sets[[2]],rownames(acute_ps))
enrich_res = extract_top_go_results(metafor_sets_enrichments)
get_most_sig_enrichments_by_groups(enrich_res)
save(metafor_gene_sets,metafor_sets_enrichments,acute_ps,longterm_ps,file="metafor_gene_sets.RData")
######################################################################################################
######################################################################################################

# Main display items
load("PADB_dataset_level_meta_analysis_data.RData")
load("tissue_expression_scores.RData")
load("metafor_gene_sets.RData")
gs = tissue_expression_scores$shared_genes
muscle_expression_scores = tissue_expression_scores$longterm$muscle[gs] + 
  tissue_expression_scores$acute$muscle[gs]
muscle_expression_scores = muscle_expression_scores/2
# show the effect of "unxepressed genes"
ps1 = acute_ps[gs,2]
ps2 = acute_ps[gs,4]
par(mfrow=c(1,3))
hist(c(runif(10000),runif(300)/100),main="Expected behaviour",xlab="p-value")
hist(ps1,main="RE",xlab="p-value")
hist(ps2,main="ME",xlab="p-value")
cor.test(-log(ps1+1e-20),muscle_expression_scores)$p.value
par(mfrow=c(1,1))
get_lfdrs(ps1,plot=1)
fdrs = get_lfdrs(ps1,plot=1,nulltype=1)
table(fdrs[,2]<0.1 & ps1<0.01)

# number of selected genes
m = c(length(metafor_gene_sets[[1]]$`acute,random_effects_muscle`),0)
m = rbind(m,c(length(metafor_gene_sets[[2]]$`acute,random_effects_muscle`),
              length(metafor_gene_sets[[2]]$`acute,mixed_effects_muscle`)))
rownames(m) = c("BH","locfdr")
colnames(m) = c("RE","ME")
barplot(m,beside=T,legend=T,args.legend = c(x="top"),ylab="Number of selected genes")

# Gene numbers and specific examples
gene_nums = t(t(sapply(metafor_gene_sets[[2]],length)))
enrich_res = extract_top_go_results(metafor_sets_enrichments)
get_most_sig_enrichments_by_groups(enrich_res)


######## "Advanced" analyses: functions ############
# Get average effect per time point x training
get_subset_forest_plot<-function(gdata,tissue="all",training="all",sortby = "time",labelsby=c("gse","time"),...){
  if(!tissue=="all"){
    gdata = gdata[gdata$tissue==tissue,]
  }
  if(!training=="all"){
    gdata = gdata[gdata$training==training,]
  }
  ord = order(gdata[,sortby])
  gdata = gdata[ord,]
  res = rma.mv(yi,vi,random = ~ 1|gse,data=gdata)
  res0 = rma(yi,vi,data=gdata)
  rtest = regtest(res0)
  print(rtest)
  forest(res,slab=apply(gdata[,labelsby],1,paste,collapse=','),...)
}

get_gene_weighted_avg_pattern <-function(gdata){
  gdata$wt = 1/(sqrt(gdata$vi))
  res = by(gdata, paste(gdata$tissue,gdata$training,gdata$time,sep=','), function(x) weighted.mean(x$yi, x$wt),simplify = T)
  v = as.numeric(res);names(v) = names(res)
  return(v)
}
reorder_weighted_avg_matrix<-function(x){
  cn = colnames(x)
  cn = sapply(cn,function(x)strsplit(x,split=',')[[1]])
  ord = order(cn[1,],cn[2,],as.numeric(cn[3,]))
  x = x[,ord]
  return(x)
}
# geneset - a list of gene symbols
print_drem_matrices<-function(x,dir_name,gene_conv = entrez2symbol,min_time_points=3,geneset=NULL){
  cn = colnames(x)
  cn = sapply(cn,function(x)strsplit(x,split=',')[[1]])
  tb = table(cn[1,],cn[2,])
  for(tissue in rownames(tb)){
    xx = x[,cn[1,]==tissue]
    times = as.numeric(cn[3,cn[1,]==tissue])
    trs = cn[2,cn[1,]==tissue]
    for(tr in unique(trs)){
      if(sum(trs==tr)<min_time_points){next}
      curr_xx = xx[,trs==tr]
      curr_times = times[trs==tr]
      colnames(curr_xx) = paste(curr_times,"h",sep='')
      rownames(curr_xx) = gene_conv[rownames(curr_xx)]
      if(!is.null(geneset)){curr_xx = curr_xx[intersect(rownames(curr_xx),geneset),]}
      fname = paste(dir_name,tissue,"_",tr,'.txt',sep='')
      curr_xx = cbind(rownames(curr_xx),curr_xx)
      colnames(curr_xx)[1] = "UID"
      write.table(curr_xx,file=fname,quote=F,sep="\t",col.names = T,row.names = F)
    }
  }
}
plot_gene_pattern<-function(x,errs = NULL,tosmooth=T,min_time_points=3,
                            tr2col=c(control="red",endurance="blue",resistance="green",both="purple"),
                            legend.x="topright",legend.cex=0.8,legend.ncol=2,
                            mfrow = c(1,2),main_prefix = "acute",y_lim_add = 0.5,y_lim_min = 1.5,
                            ...){
  cn = names(x)
  cn = sapply(cn,function(x)strsplit(x,split=',')[[1]])
  tb = table(cn[1,],cn[2,])
  if(!is.null(mfrow)){par(mfrow=mfrow)}
  for(tissue in rownames(tb)){
    xx = x[cn[1,]==tissue]
    if(!is.null(errs)){yy = errs[cn[1,]==tissue]}
    if(length(xx)<min_time_points){next}
    times = as.numeric(cn[3,cn[1,]==tissue])
    trs = cn[2,cn[1,]==tissue]
    ylim = c(min(xx)-y_lim_add,max(xx)+y_lim_add)
    ylim[1] = min(-y_lim_min,ylim[1])
    ylim[2] = max(y_lim_min,ylim[2])
    plot(xx,x=times,col="white",las=2,main=paste(main_prefix,tissue),xlim=c(-1,max(times)),ylim=ylim,ylab="log2 Fold Change",xlab="time")
    for(tr in unique(trs)){
      if(tosmooth){xx1 = c(0,smooth(xx[trs==tr]))}
      else{xx1 = c(0,xx[trs==tr])}
      if(length(xx1)<min_time_points){next}
      tt1 = c(-1,times[trs==tr])
      lines(xx1,x=tt1,type='l',lwd=2,col=tr2col[tr],pch=20)
      if(!is.null(errs)){
        yy1 = c(0.01,yy[trs==tr])
        arrows(tt1, xx1-yy1, tt1, xx1+yy1, length=0.05, angle=90, code=3,col=tr2col[tr])
      }
    }
    legend(x=legend.x,legend = unique(trs),fill=tr2col[unique(trs)],cex=legend.cex,ncol=legend.ncol)  
  }
}
library(cluster)
pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
# x: rows are genes
perform_gene_clustering<-function(x,fun=pam1,standardize_genes=T,num_pcs=5){
  xx = x
  if(standardize_genes){xx = t(apply(x,1,function(x)(x-mean(x))/(sd(x))))}
  if(num_pcs>1){xx = prcomp(xx)$x[,1:num_pcs]}
  gs.pam.xx <- clusGap(xx, FUN = fun, K.max = 10, B = 20)
  num_c = maxSE(gs.pam.xx[[1]][,3],gs.pam.xx[[1]][,4])
  cc = fun(xx,num_c)$cluster
  return(cc)
}
cluster_homogeneities<-function(x,cl,...){
  hs = c()
  for(cc in unique(cl)){
    cc = as.character(cc)
    if(sum(cl==cc)==1){
      hs[cc] = 1
    }
    else{
      rhos = cor(t(x[cl==cc,]),...)
      hs[cc] = mean(rhos[lower.tri(rhos)])
    }
  }
  names(hs) = as.character(unique(cl))
  return(hs)
}
standardize_rows = function(x){
  xx = t(apply(x,1,function(x)(x-mean(x))/(sd(x))))
  return(x)
}
options(expressions=10000)
cluster_genes_by_homogeneity <- function(x,thr=0.5,func=kmeans,...){
  if(is.null(dim(x))){
    return (1)
  }
  if(nrow(x)<2){
    return(1)
  }
  v = rep(1,nrow(x));names(v) = rownames(x)
  hs = cluster_homogeneities(x,v)
  if(hs>=thr){
    return(v)
  }
  k=2
  while(all(hs<thr) && k < nrow(x)){
    v = as.numeric(func(x,k)[[1]])
    names(v) = rownames(x)
    hs = cluster_homogeneities(x,v)
    k = k+1
  }
  if(k==nrow(x)){
    v = 1:nrow(x);names(x) = rownames(x)
    return(v)
  }
  for(cl in names(hs)){
    if(hs[cl]>=thr){next}
    newx = x[v==cl,]
    newx_v = as.numeric(cluster_genes_by_homogeneity(newx,thr,func))
    newx_v = as.numeric(newx_v) + max(as.numeric(v))
    names(newx_v) = rownames(newx)
    v[names(newx_v)] = newx_v
    # print(table(v))
    # print("###")
    # print(table(newx_v))
    # print("###")
    # print(cluster_homogeneities(newx,newx_v))
  }
  v = as.numeric(as.factor(v))
  names(v) = rownames(x)
  return(v)
}

# mm = standardize_rows(m)
# vv = cluster_genes_by_homogeneity(mm,0.5,pam1)
# sort(table(vv))
# cluster_homogeneities(mm,vv)

#########################################################
library(org.Hs.eg.db)
entrez2symbol = as.list(org.Hs.egSYMBOL)

weighted_avg_matrices=list()
weighted_avg_matrices[["acute"]] = t(sapply(acute_gene_tables_raw,get_gene_weighted_avg_pattern))
weighted_avg_matrices[["longterm"]] = t(sapply(longterm_gene_tables_raw,get_gene_weighted_avg_pattern))
# reorder the matrices
weighted_avg_matrices = lapply(weighted_avg_matrices,reorder_weighted_avg_matrix)
dir.create("processed_avg_effects_matrices")
dir.create("processed_avg_effects_matrices/acute")
dir.create("processed_avg_effects_matrices/longterm")
print_drem_matrices(weighted_avg_matrices$acute,"processed_avg_effects_matrices/acute/")
print_drem_matrices(weighted_avg_matrices$longterm,"processed_avg_effects_matrices/longterm/")
metafor_gene_sets_names_0.2 = lapply(metafor_gene_sets[[3]],function(x,y)sort(unlist(y[x])),y=entrez2symbol)
metafor_gene_sets_names_0.2_all_genes = sort(unique(unlist(metafor_gene_sets_names_0.2)))

# Specific examples
metafor_gene_sets_names = lapply(metafor_gene_sets[[1]],function(x,y)sort(unlist(y[x])),y=entrez2symbol)
pb_genes = lapply(publication_bias_res,function(x)names(which(x[1,]<0.01)))
pb_genes = lapply(pb_genes,function(x,y)sort(unlist(y[x])),y=entrez2symbol)
all_genes = sort(unique(unlist(metafor_gene_sets_names)))
known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
                "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")
which(sapply(known_genes,function(x,y)any(grepl(x,y)),y=all_genes))
intersect(all_genes,known_genes)
intersect(known_genes,unlist(entrez2symbol[rownames(acute_ps)]))
sapply(metafor_gene_sets_names,intersect,y=known_genes)
sapply(pb_genes,intersect,y=known_genes)

gene = "10891"
gdata = acute_gene_tables[[gene]] # PGC1 in acute response
get_subset_forest_plot(gdata,"muscle",main="PGC1, acute, muscle")
get_subset_forest_plot(gdata,"blood",main="PGC1,acute,blood")
acute_ps[gene,]
acute_gene_tables[[gene]]
acute_gene_tables_raw[[gene]]
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = T,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL)

gene = "1282"
entrez2symbol[gene]
gdata = acute_gene_tables[[gene]] # COL4A1 in acute response
get_subset_forest_plot(gdata,"muscle",main="COL4A1, acute, muscle")
gdata = longterm_gene_tables[[gene]] # COL4A1 in longterm response
get_subset_forest_plot(gdata,"muscle",main="COL4A1, longterm, muscle")
acute_ps[gene,]
longterm_ps[gene,]
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = T,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL)

gene = "4318" # MMP9
entrez2symbol[gene]
gdata = acute_gene_tables[[gene]]
get_gene_analysis_pvals_with_gse_correction(gdata[gdata$tissue=="blood",])
get_subset_forest_plot(gdata,"blood",main="MMP9, acute, blood")
acute_ps[gene,]
longterm_ps[gene,]
gdata$time = ordered(gdata$time)
get_gene_analysis_pvals_with_gse_correction(gdata[gdata$tissue=="blood",])
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = T,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL)

gene = "7139" # TNNT2 - expected but has poor signal
entrez2symbol[gene]
gdata = acute_gene_tables[[gene]] 
get_subset_forest_plot(gdata,"muscle")
get_subset_forest_plot(gdata,"blood")
gdata = longterm_gene_tables[[gene]] 
get_subset_forest_plot(gdata,"muscle")
get_subset_forest_plot(gdata,"blood")
acute_ps[gene,]
longterm_ps[gene,]
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = T,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL)

gene = "7070" # THY1
gdata = longterm_gene_tables[[gene]] 
get_subset_forest_plot(gdata,"muscle")
get_subset_forest_plot(gdata,"blood")
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = F,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL,tosmooth = F)

gene = "70" # ACTC1
gdata = longterm_gene_tables[[gene]] 
get_subset_forest_plot(gdata,"muscle")
gdata = acute_gene_tables[[gene]] 
get_subset_forest_plot(gdata,"muscle",main = "ACTC1, acute, muscle")
plot_gene_pattern(weighted_avg_matrices$acute[gene,],tosmooth = T,mfrow=c(2,2))
plot_gene_pattern(weighted_avg_matrices$longterm[gene,],main_prefix = "long-term",mfrow=NULL,tosmooth = T)

gdata = longterm_gene_tables_simpletime[[gene]]
gdata = acute_gene_tables[[gene]]
gdata = gdata[gdata$tissue=="blood",]
gdata = gdata[gdata$tissue=="muscle",]
res1 = func(yi,vi,mods = ~ training + time ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000))
res2 = func(yi,vi,mods = ~ training + ordered(time) ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000))
res3 = func(yi,vi,1/sqrt(vi),data=gdata, random = ~ 1|gse, control=list(maxiter=10000))
res4 = func(yi,vi,1/sqrt(vi),mods = ~ training ,data=gdata, random = ~ 1|gse, control=list(maxiter=10000))
summary(res2)[[2]]
publication_bias_res[[1]][,gene]
anova(res1,res3)
summary(lm(yi~time+training+factor(gse),data=gdata))

table(gdata$training,gdata$time)
table(gdata$training,gdata$gse)

# # simulations, metafor
# ns = 10
# yi = rt(ns,df=10) 
# vi = rchisq(ns,df=10)/10
# yi[1] = yi[1] + 10
# yi[3] = yi[3] + 10
# gse = 1:10
# m = rep(0,10)
# m[1:3] = 1
# rma.mv(yi,vi,random=~1|gse,mods = ~m)
# rma.mv(yi,vi,random=~1|m,mods = ~gse)
# rma.mv(yi,vi,random=~1|m)
# regtest(rma(yi,vi))
# funnel(rma.mv(yi,vi))
# forest(rma(yi,vi))

# Clustering and plots
# Step 1: take a selected list of genes, in a specific tissue 
genes = metafor_gene_sets$`0.2lfdr`$`acute,mixed_effects_muscle`
genes = unique(unlist(metafor_gene_sets[[2]]))

genes = intersect(names(acute_gene_tables),genes)
tissue = "muscle"
tissues = acute_gene_tables_raw[[genes[1]]]$tissue
gene_patterns = sapply(acute_gene_tables_raw[genes],function(x)x$yi)[tissues==tissue,]
colnames(gene_patterns) = entrez2symbol[genes]

# Step 2: order by time and training
times = acute_gene_tables_raw[[genes[1]]]$time[tissues==tissue]
trs = acute_gene_tables_raw[[genes[1]]]$training[tissues==tissue]
ord = order(times,trs)
times=times[ord];trs=trs[ord];gene_patterns=t(gene_patterns[ord,])
gene_clusters = perform_gene_clustering(gene_patterns,num_pcs=-1,standardize_genes = T)
cluster_homogeneities(gene_patterns,gene_clusters,method='spearman')
table(gene_clusters)

tr2col = c("red","blue","green");names(tr2col) = unique(trs)
tr2lty = c(2,1,1);names(tr2lty) = unique(trs)
num_c = length(unique(gene_clusters))
num_c=8
par(mfrow=c(2,num_c/2))
for(i in 1:num_c){
  cl_data = gene_patterns[gene_clusters==i,]
  plot(colMeans(cl_data),x=times,col="white",las=2)
  for(tr in unique(trs)){
    curr_data = cl_data[,trs==tr]
    print(intersect(rownames(curr_data),known_genes))
    curr_times = times[trs==tr]
    curr_profile = get_avg_merged_pattern(curr_data,curr_times)
    lines(curr_profile,x=unique(curr_times),ylim=c(-1,1),main=tr,type='b',
          lwd=2,col=tr2col[tr],lty=tr2lty[tr],pch=20)
    abline(h=0)
  }
  legend(x="top",legend = unique(trs),fill=tr2col,cex=1)  
}

# Step 3: cluster


xx = ordered(paste(times,trs,sep=","))
plot(gene_patterns[,1],type='b',ylim=c(min(gene_patterns),max(gene_patterns)),x=xx)
for(j in 1:ncol(gene_patterns)){
  lines(gene_patterns[,j],x=xx,type='l')
}

# x: rows are genes
get_avg_merged_pattern<-function(x,times,smooth=F){
  xx = t(apply(x,1,function(x,y)tapply(x,factor(times),mean),y=times))
  xx_means = colMeans(xx)
  if(smooth){xx_means = smooth.spline(xx_means)}
  return(xx_means)
}

library(gplots)
heatmap.2(cor(t(gene_patterns),method = "spearman"),trace="none",
          scale="none",col=colorRampPalette(c("red","white","blue"))(256))


# # Tests and comments from the paper of metafor (2010)
# gdata = acute_gene_tables_raw[["10891"]] # PGC1 in acute response
# gdata = longterm_gene_tables[["1282"]] # COL1 gene
# gdata = longterm_gene_tables[["4168"]] # A negative example
# gdata = acute_gene_tables[["5166"]] # survived rep but not meta
# # gdata = acute_gene_tables[["11326"]] # gene with significant modifiers
# # another test: make a dataset with large effects
# # gdata$yi = rnorm(nrow(gdata),mean=5)
# plot(gdata[,"yi"],gdata[,"vi"],pch=as.numeric(as.factor(gdata$tissue)))
# # knha - a correction that accounts for the uncertainty in the random effect
# res0 = rma.mv(yi,vi,random = ~ 1|gse,data=gdata,subset = (tissue=="muscle"))
# # explanation of the result above:
# #   mu - the average effect is 0.067, the CI contains zero
# res = rma.mv(yi,vi,mods = ~  training + time, random= ~1|gse ,data=gdata[gdata$tissue=="muscle",])
# # CIs of the anova stats:
# confint(res0)
# # Forest plot - very informative
# forest(res)
# # difference in tau before and after using moderators - teaches us about the
# # percentage of explained variance due to the moderators
# # when the test for residuals (QE) is signficant - we may be missing additional
# # moderators
# # We can use predict to get expected effects for new moderators:
# # Currently does not work because factors should be transformed into dummy vars.
# # Also pages 18-19 show nice figures and analysis of the predictions.
# predict(res,newmods = as.matrix(data.frame("",time=5,"",training="endurance",tissue="muscle")), addx = TRUE)
# predict(res,newmods = as.matrix(gdata))
# # look at the fitted values
# predict(res)
# # Separate by tissue
# res2 = rma(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, subset = (tissue=="muscle"))
# res1 = rma(yi=yi,vi=vi,mods = ~ training + time ,data=gdata,knha=T, subset = (tissue=="blood"))
# # Residual analysis for detecting outlier datasets
# barplot(as.numeric(rstudent(res)$z));abline(-2,0);abline(2,0)
# # Residual analysis is informative but not enough
# # case deletion diagnostics are informatice as well
# plot(influence(res))
# # funnel plots
# funnel(res0)
# # radial plots: useful for consistency analysis
# # can be used only for models without moderators
# radial(res0)
# # qq plots for the standardized residuals
# qqnorm(res0,main="random")
# qqnorm(res,main="mixed, both tissues")
# qqnorm(res1,main="blood")
# qqnorm(res2,main="muscle")
# # tests for publication bias
# regtest(res0,predictor="vi",model="lm")
# regtest(res2,predictor="vi",model="lm")
# # anova tests
# anova(res0,res)
# 
# # # Compare to the rmeta package
# # install.packages('rmeta')
# # library(rmeta)
# # blood_gdata = gdata[gdata$tissue=="blood",]
# # rmeta_res0 = meta.summaries(d=blood_gdata$yi,se = blood_gdata$vi)
# # summary(rmeta_res0)[[3]]
# # plot(rmeta_res0)
# # funnelplot(rmeta_res0)
# # gdata$yi
# # 
# # # Tests on simulated data
# # n = 5;effect = 3; hetero=0.5; effect2=6
# # vi = rep(0.5,n)
# # yi = rnorm(n,sd=vi) + effect + rnorm(n,sd=hetero)
# # rma(yi,vi)$pval
# # dummy = 1:n
# # rma.mv(yi,vi,random=~ 1|dummy)$pval
# # # vs.
# # rma(c(yi,rnorm(n,sd=vi) + rnorm(n,sd=hetero)),c(vi,vi))$pval
# # # vs.
# # ref = factor(c(rep(1,n),rep(2,n)))
# # rma.mv(c(yi,rnorm(n,sd=vi) + rnorm(n,sd=hetero)),c(vi,vi),random = ~1|ref)$pval
# # 
# # # random effects
# # beta1 = 1 ; beta2 = 0.5
# # vi = c(vi,vi)
# # yi = c(yi+beta1,yi+beta2)
# # ref = factor(c(rep(1,n),rep(2,n)))
# # rma(yi,vi)
# # rma.mv(yi,vi,rand = ~1|ref)










