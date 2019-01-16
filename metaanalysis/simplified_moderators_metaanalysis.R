# Our algorithm for analysis of a single gene
# Input: datasets for time point t and gene g
# 1. Run meta-analysis for endurence and resistance (and both)
# 2. Run meta-analysis for controls
# 3. Get controls intercept estimation beta_c
# 4. Get current gene's overall diff expression: beta_e = max(beta_0,beta_0 + beta_endurence)
# 5. SCORE 1: abs(beta_e-beta_c)
# 6. SCORE 2: Egger's test of the exercise meta-analysis
# 7. Return SCORE 1, SCORE 2, and the betas (and their significance) from the exercise meta-analysis

setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(org.Hs.eg.db);library(metafor)
source('/Users/David/Desktop/repos/motrpac_public_data_analysis/metaanalysis/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

############################################################################
############################################################################
############################################################################
# Functions for data filtering

#' Change the time field from numeric to ordered factor.
#' @param gdata A data frame with a column "time"
#' @return A new data frame.
simplify_time_in_gdata<-function(gdata,func=binarize_feature,...){
  gdata$time = func(gdata$time,...)
  gdata$time = ordered(gdata$time)
  return(gdata)
}
binarize_feature<-function(tt,thr=5){return(as.numeric(tt <= thr))}
simplify_age_gdata<-function(gdata,func=binarize_feature,...){
  gdata$avg_age = func(as.numeric(gdata$avg_age),...)
  gdata$avg_age = ordered(gdata$avg_age)
  return(gdata)
}

# Clean the datasets: from muscle, remove 1 cohort with time point zero
# From long term exclude all studies with time > 150 days
clean_acute_table <-function(gdata,tissue="muscle"){
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  if(tissue=="muscle"){
    gdata = gdata[as.numeric(gdata$time) > 1,]
  }
  gdata = simplify_time_in_gdata(gdata,simplify_time_acute)
  gdata = simplify_age_gdata(gdata,thr=40)
  return(gdata)
}
clean_longterm_table <-function(gdata,tissue="muscle"){
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  gdata = gdata[as.numeric(gdata$time) < 150 ,]
  gdata = simplify_age_gdata(gdata,thr=40)
  return(gdata)
}
get_untrained_table <-function(gdata,tissue="muscle"){
  gdata = gdata[gdata$tissue==tissue,]
  gdata$vi = pmax(gdata$vi,1e-5)
  rows = is.element(gdata$training,set=c("yoga","control","untrained"))
  return(gdata[rows,])
}

############################################################################
############################################################################
############################################################################
# Functions for the meta-analysis
# ... specifies random effects formula and structure
model_selection_meta_analysis<-function(gdata,...){
  mod_names = c("training","time","avg_age","prop_males")
  gdata_basic = gdata[,c("yi","vi","V1","gse")]
  gdata_mods = gdata[,mod_names]
  to_rem = apply(gdata_mods,2,function(x)length(unique(x,na.rm=T))<2)
  gdata_mods = gdata_mods[,!to_rem]
  mod_names = colnames(gdata_mods)
  models = list(); aics = c()
  N = 2^length(mod_names)
  for(j in 1:N){
    curr_group = as.integer(intToBits(j-1))[1:length(mod_names)]
    curr_mods = mod_names[curr_group==1]
    curr_frm=NULL;curr_model = NULL;curr_name="base_model"
    if(length(curr_mods)>0){
      curr_frm = as.formula(paste("~",paste(curr_mods,collapse="+")))
      curr_name = paste(curr_mods,collapse=";")
    }
    curr_model = meta_analysis_wrapper(gdata,mods=curr_frm,...)
    if(is.element("rma",class(curr_model))){
      models[[curr_name]] = curr_model
      aics[curr_name] = fitstats(curr_model)[5,1]
    }
  }
  return(list(models=models,aics=aics))
}
select_model_return_p<-function(res,aic_thr=2){
  aics_thr = (res$aics[1] - res$aics) > aic_thr
  best_ind = which(res$aics == min(res$aics,na.rm=T) & aics_thr)[1]
  if(is.na(best_ind) || length(best_ind)==0){best_ind=1}
  best_m = res$models[[best_ind]]
  best_m_name = names(res$models)[best_ind]
  sapply(res$models,function(x)x$QMp)
  pval = best_m$QMp
  return(list(pval=pval,name=best_m_name,model=best_m,aic = min(res$aics,na.rm=T)))
}

meta_analysis_wrapper<-function(gdata,func = rma.mv,...){
  for(rel.tol in c(1e-8,1e-7,1e-6)){
    cc=list(iter.max=10000,rel.tol=rel.tol)
    res = NULL
    try({
      res = func(yi,vi,data=gdata,...)
    })
    if(!is.null(res)){return(res)}
  }
  return(NA)
}
num_params<-function(model){
  x = fitstats(model)[,1]
  k = (2*x[1]+x[3])/2
  return(unname(k))
}

add_prefix_to_names<-function(pref,l,sep=":"){
  names(l) = paste(pref,names(l),sep=sep)
  return(l)
}

acute_gdata_metaanalysis<-function(gdata,permtest=F){
  res1 = model_selection_meta_analysis(gdata,random=list(~ V1|gse),struct="AR")
  res0 = model_selection_meta_analysis(gdata,func=rma.uni)
  l = list(
    models = c(add_prefix_to_names("simple",res0$models),
               add_prefix_to_names("time_ar",res1$models)),
    aics = c(add_prefix_to_names("simple",res0$aics),
             add_prefix_to_names("time_ar",res1$aics))
  )
  sel = select_model_return_p(l)
  if(permtest && is.element("rma.uni",set=class(sel$model))){
    sel[["permp"]] = permutest(sel$model)
  }
  return(list(all_models=l,selected=sel))
}

longterm_gdata_metaanalysis<-function(gdata,permtest=F){
  res1 = model_selection_meta_analysis(gdata,random=list(~ 1|gse))
  res0 = model_selection_meta_analysis(gdata,func=rma.uni)
  l = list(
    models = c(add_prefix_to_names("simple",res0$models),
               add_prefix_to_names("time_ar",res1$models)),
    aics = c(add_prefix_to_names("simple",res0$aics),
             add_prefix_to_names("time_ar",res1$aics))
  )
  sel = select_model_return_p(l)
  if(permtest && is.element("rma.uni",set=class(sel$model))){
    sel[["permp"]] = permutest(sel$model)
  }
  return(list(all_models=l,selected=sel))
}

simple_stouffer_meta_analysis<-function(gdata){
  gdata = gdata[!is.na(as.numeric(gdata$p)),]
  p = metap::sumz(as.numeric(gdata$p),weights = as.numeric(gdata$df))$p[1,1]
  return(p)
}

run_simple_re<-function(gdata){return(meta_analysis_wrapper(gdata,func=rma.uni))}
try_get_field<-function(res,fname){
  if(is.element(fname,set=names(res))){
    return(res[[fname]])
  }
  return(NA)
}
pvalue_qqplot<-function(ps,...){
  x = runif(5000)
  qqplot(y=-log(ps,10),x=-log(x,10),xlab="Theoretical",ylab="Observed",...)
  abline(0,1,lty=2,lwd=2)
}

############################################################################
############################################################################
############################################################################

# Prepare the datasets for the different analyses below
# Load the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
load("PADB_dataset_level_meta_analysis_data.RData")

# Get the cleaned, filtered datasets
datasets = list()
datasets[["acute,muscle"]] = lapply(acute_gene_tables,clean_acute_table,tissue="muscle")
datasets[["acute,blood"]] = lapply(acute_gene_tables,clean_acute_table,tissue="blood")
datasets[["longterm,muscle"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="muscle")
datasets[["longterm,blood"]] = lapply(longterm_gene_tables,clean_longterm_table,tissue="blood")

# Reshape data for replication analysis
rep_datasets = lapply(datasets,function(x)t(sapply(x,function(y)as.numeric(y$p))))
for(nn in names(rep_datasets)){
  colnames(rep_datasets[[nn]]) = 
    paste(rownames(datasets[[nn]][[1]]),datasets[[nn]][[1]]$V1,
          datasets[[nn]][[1]]$time,sep=";")
}
sapply(rep_datasets,dim)
sapply(rep_datasets,rownames)
sapply(rep_datasets,colnames)
# Get the untrained controls datasets
untrained_datasets = list()
untrained_datasets[["acute,muscle"]] = lapply(acute_gene_tables_raw,get_untrained_table,tissue="muscle")
untrained_datasets[["acute,blood"]] = lapply(acute_gene_tables_raw,get_untrained_table,tissue="blood")
untrained_datasets[["longterm,muscle"]] = lapply(longterm_gene_tables_raw,get_untrained_table,tissue="muscle")
untrained_datasets[["longterm,blood"]] = lapply(longterm_gene_tables_raw,get_untrained_table,tissue="blood")
sapply(untrained_datasets,function(x)dim(x[[1]]))
############################################################################
############################################################################
############################################################################
# Simple random effects analysis - to understand the data
simple_REs = lapply(datasets, function(x)lapply(x,run_simple_re))
simple_RE_pvals = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="QMp"))
simple_RE_I2s = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="I2"))
simple_RE_tau2s = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="tau2"))
simple_RE_beta = lapply(simple_REs,function(x)sapply(x,try_get_field,fname="beta"))
st_pvals = lapply(datasets, function(x)lapply(x,simple_stouffer_meta_analysis))
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  pvalue_qqplot(simple_RE_pvals[[nn]],main=nn,pch=20,cex=0.5)
}
for(nn in names(simple_RE_pvals)){
  r = hist(simple_RE_I2s[[nn]],plot=F)
  per = sum(simple_RE_I2s[[nn]]>70,na.rm = T)/length(simple_RE_I2s[[nn]])
  per = format(per*100,digits=3)
  cols=rep("white",length(r$mids))
  cols[r$mids > 70] = "blue"
  hist(simple_RE_I2s[[nn]],main=paste(nn,"(",per,"%)",sep=""),xlab = "I^2(%)",col=cols)
}
for(nn in names(simple_RE_pvals)){
  pvalue_qqplot(unlist(st_pvals[[nn]]),main=nn,pch=20,cex=0.5)
}
for(nn in names(simple_RE_pvals)){
  x1 = simple_RE_I2s[[nn]]
  x2 = log(simple_RE_pvals[[nn]])
  inds = !is.na(x1) & !is.na(x2)
  print(nn)
  print(cor(x1[inds],x2[inds],method="spearman"))
}
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  x1 = simple_RE_I2s[[nn]]
  x2 = simple_RE_tau2s[[nn]]
  inds = !is.na(x1) & !is.na(x2)
  print(nn)
  print(cor(x1[inds],x2[inds],method="spearman"))
  plot(x=x1,y=x2,xlab="I2(%)",ylab="Tau^2",pch=20,cex=0.2,main=nn)
}
par(mfrow=c(2,2))
for(nn in names(simple_RE_pvals)){
  x1 = simple_RE_I2s[[nn]]
  x2 = simple_RE_beta[[nn]]
  inds = !is.na(x1) & !is.na(x2)
  print(nn)
  print(cor(x1[inds],x2[inds],method="spearman"))
  plot(x=x1,y=x2,xlab="I2(%)",ylab="beta",pch=20,cex=0.2,main=nn)
}

# Example for high I2 and low tau
# TODO: revise and add some stats to the report
nn  = 2
selected_is = names(which(simple_RE_I2s[[nn]] > 60 &
    simple_RE_tau2s[[nn]] < 0.1 & abs(simple_RE_beta[[nn]]) > 0.25))
selected_i = sample(selected_is)[1]
forest(simple_REs[[nn]][[selected_i]])
median(as.numeric(datasets[[nn]][[selected_i]]$p))
screen_res_znormix[[nn]][selected_i,]
colSums(screen_res_znormix[[nn]]<0.25)
colSums(screen_res_locfdr[[nn]]<0.25)

save.image(file="workspace_before_rep_analysis.RData")

############################################################################
############################################################################
############################################################################
# Replication analysis
# Naive analysis
naive_rep_analysis<-function(pvals,gses,thr=0.05,nrep=2){
  v = pvals <= thr
  if(sum(v,na.rm=T)==0){return(F)}
  tt = table(v,gses)["TRUE",]
  return(sum(tt>0) >= nrep)
}

naive_rep_analysis_results = list()
for(nn in names(rep_datasets)){
  gses = datasets[[nn]][[1]]$gse
  currn = ceiling(length(unique(gses))/2)
  print(currn)
  naive_rep_analysis_results[[nn]] = apply(rep_datasets[[nn]],1,
      naive_rep_analysis,gses=gses,thr=0.05,nrep=currn)
}
sapply(naive_rep_analysis_results,table)
sapply(naive_rep_analysis_results,function(x)table(x)/length(x))
sapply(naive_rep_analysis_results,function(x,y)x[y],y="10891")
sapply(naive_rep_analysis_results,function(x,y)x[y],y="7139")
sapply(naive_rep_analysis_results,function(x,y)x[y],y="1277")
sapply(naive_rep_analysis_results,function(x,y)x[y],y="70")

# # Load repfdr results and compare
# scr_path = "/Users/David/Desktop/MoTrPAC/PA_database/screen_res/"
# pvals_files = list.files(scr_path)
# pvals_files = pvals_files[grepl("pvals.txt",pvals_files)]
# screen_results = list()
# for (ff in pvals_files){
#   currname = gsub("_pvals.txt","",ff)
#   screen_results[[currname]] = list()
#   for(m in c("bum","znormix")){
#     outfile = paste(scr_path,currname,"_",m,".txt",sep="")
#     screen_results[[currname]][[m]] = read.delim(outfile,row.names = 1,header=T)
#   }
# }
# names(screen_results) = gsub("_",",",names(screen_results))
# save(screen_results,file=paste(scr_path,"screen_results.RData",sep=""))
# for(nn in names(screen_results)){
#   gses = datasets[[nn]][[1]]$gse
#   currn = ceiling(length(unique(gses))/2)
#   curr_genes1 = names(which(naive_rep_analysis_results[[nn]]))
#   x1 = screen_results[[nn]]$bum[curr_genes1,currn]
#   x2 = screen_results[[nn]]$znormix[curr_genes1,currn]
#   max_fdr_gene = curr_genes1[x2==max(x2)][1]
#   print(screen_results[[nn]]$bum[max_fdr_gene,])
#   print(screen_results[[nn]]$znormix[max_fdr_gene,])
#   print(rep_datasets[[nn]][max_fdr_gene,])
#   plot(x1,x2)
# }
# sapply(screen_results,function(x)colSums(x$bum < 0.1))

# # Some QA
# nn = "longterm,muscle"
# bum_res = screen_results[[nn]]$bum[,7]
# names(bum_res) = rownames(screen_results[[nn]]$bum)
# plot(bum_res,rowMeans(rep_datasets[[nn]]))
# hist(bum_res[naive_rep_analysis_results[[nn]]])
# xx = rep_datasets[[nn]][names(bum_res[bum_res>0.19 & bum_res < 0.2]),]
# xx = rep_datasets[[nn]][names(which(naive_rep_analysis_results[[nn]])),]
# library(gplots)
# heatmap.2(xx)
# dim(xx)

# # Select genes, look at known genes and enrichments
# rep_gene_sets_names = lapply(rep_gene_sets,function(x,y)sort(unlist(y[x])),y=entrez2symbol)
# known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
#                 "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")
# sapply(rep_gene_sets_names,intersect,y=known_genes)
# sapply(naive_rep_analysis_results,function(x,y)x[y],y=known_genes)

# Save the data needed for the metanaanalysis and model selection below
save(naive_rep_analysis_results,datasets,untrained_datasets,file="meta_analysis_input.RData")


############################################################################
############################################################################
############################################################################
# Meta-regression and model selection for selected genes
library(parallel)

# all_meta_analysis_res <- list()
# for(nn in names(naive_rep_analysis_results)){
#   curr_dataset = datasets[[nn]][naive_rep_analysis_results[[nn]]]
#   if(grepl("acute",nn)){
#     analysis1 = mclapply(curr_dataset,acute_gdata_metaanalysis,mc.cores = 4)
#   }
#   else{
#     analysis1 = mclapply(curr_dataset[1:3],longterm_gdata_metaanalysis,mc.cores = 4)
#   }
#   analysis2 = unlist(mclapply(curr_dataset,simple_stouffer_meta_analysis,mc.cores=4))
#   all_meta_analysis_res[[nn]] = list(model_selection = analysis1,simple_stouffer = analysis2)
#   forest(analysis1[[3]]$selected$model)
# }
# save(all_meta_analysis_res,naive_rep_analysis_results,file="meta_analysis_results.RData")

# Load the results (run on sherlock) instead of running the code above
load("meta_analysis_results.RData")

for(nn in names(all_meta_analysis_res)){
  analysis1 = all_meta_analysis_res[[nn]][[1]]
  ps1 = sapply(analysis1,function(x)x$selected$pval)
}

ps1 = sapply(analysis1,function(x)x$pval);abline(0,1)
cor(-log(ps1),-log(analysis2),method="spearman")
qqplot(-log(ps1),-log(analysis2))
hist(ps1)
qqplot(-log(ps1),x=-log(runif(1000)));abline(0,1)
table(sapply(analysis1,function(x)x$name))
min_i = which(p.adjust(analysis2)<0.01)
sorted_diffs = sort(ps1-analysis2,decreasing=F)[1:10]
gdata = dataset[["54407"]]
res = acute_gdata_metaanalysis(gdata)
res$pval
analysis2[["54407"]]
forest(res$model)
metap::allmetap(as.numeric(gdata$p),method="all")
res$name
class(res$model)
permutest(res$model)
hist(analysis2)
qqplot(-log(analysis2),x=-log(runif(1000)));abline(0,1)
qqplot(-log(pmax(ps1,analysis2)),x=-log(runif(1000)));abline(0,1)
table(p.adjust((pmax(ps1,analysis2)),method='fdr')<0.1)
which(p.adjust((pmax(ps1,analysis2)),method='fdr')<0.1)

# # Examine stats - we use this to make decisions 
# gd = acute_gene_tables_raw[[1]]
# tb = table(gd$time,gd$tissue)
# gd$time = simplify_time_acute(gd$time)
# tb = table(gd$time,gd$tissue)
# gd = longterm_gene_tables_raw[[1]]
# tb = table(gd$time,gd$tissue)

# # A series of helper functions for implementing our algorithm
# controls_meta_analysis<-function(gdata){
#   gdata = gdata[!grepl("treatment",gdata$training),]
#   gdata_c = gdata[grepl("untrained",gdata$training),]
#   beta_c_est = NA; beta_c_lb = NA; beta_c_ub = NA
#   if(nrow(gdata_c)>0){
#     obj_c = rma(yi,vi,weights = 1/vi,data=gdata_c)
#     beta_c_est = obj_c$beta[[1]]
#     beta_c_lb = obj_c$ci.lb[1]
#     beta_c_ub = obj_c$ci.ub[1]
#   }
#   return(c(beta_c_est,beta_c_lb,beta_c_ub))
# }

# # Some simple tests
# gdata = datasets[[1]][["9202"]]
# gdata = datasets[[3]][["55222"]]
# colnames(obj)
# # Very simple models without specifying repeated measures
# res0_1 = rma.uni(yi,vi,data = gdata,weights=1/vi)
# AIC(res0_1)
# res0_2 = rma.uni(yi,vi,data=gdata)
# AIC(res0_2)
# res0_3 = rma.mv(yi,vi,data=gdata,random=~1|gse)
# AIC(res0_3)
# res0_4 = rma.mv(yi,vi,data=gdata,random=~V1|gse,struct="AR")
# AIC(res0_4)
# permutest(res0_2)
# res0_2$QMp
# forest(res0_2)
# res_mv1 = rma.mv(yi,vi,data = gdata,mods = ~1+training+time,
#                  random = ~ V1|gse,struct="AR")
# res_mv2 = rma.mv(yi,vi,data = gdata,
#                  mods = ~1+training+time,
#                  random = list(~ V1|gse),struct="AR")
# AIC(res_mv1)
# AIC(res_mv2)
# forest(res_mv1)
# # useful functions: funnel, forest, radial, confint, summary, rstandard, rstudent, residuals
# # Add noise structure to random effects
# # The nested structure is required to have different correlated
# # random effects for time points within a cohort.
# plot(residuals(res_mv2),col=as.factor(gdata$gse))
# forest(res_mv2)
# funnel(res_mv1)
# funnel(res_mv2)



