setwd("/home/users/davidama/motrpac_metaanalysis/meta_analysis/")
library(parallel)
library(metafor,lib.loc="~/R/packages")

args = commandArgs(trailingOnly=TRUE)
num_cores = 4
if(length(args)>0){
  num_cores = as.numeric(args[1])
}

############################################################################
############################################################################
############################################################################
# TODO: merge with the functions in simplified_moderators_metaanalysis.R
# TODO: we currently have code dup
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
    curr_model = suppressWarnings(meta_analysis_wrapper(gdata,mods=curr_frm,...))
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
    },silent=TRUE)
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
  # sel = select_model_return_p(l)
  # if(permtest && is.element("rma.uni",set=class(sel$model))){
  #   sel[["permp"]] = permutest(sel$model)
  # }
  # keep the simple models and the top 2
  l = keep_main_results_for_model_list(l)
  gc()
  return(l)
}

longterm_gdata_metaanalysis<-function(gdata,permtest=F,simple_output=T){
  res1 = model_selection_meta_analysis(gdata,random=list(~ 1|gse))
  res0 = model_selection_meta_analysis(gdata,func=rma.uni)
  l = list(
    models = c(add_prefix_to_names("simple",res0$models),
               add_prefix_to_names("time_ar",res1$models)),
    aics = c(add_prefix_to_names("simple",res0$aics),
             add_prefix_to_names("time_ar",res1$aics))
  )
  # sel = select_model_return_p(l)
  # if(permtest && is.element("rma.uni",set=class(sel$model))){
  #   sel[["permp"]] = permutest(sel$model)
  # }
  # keep the simple models and the top 2
  l = keep_main_results_for_model_list(l)
  gc()
  return(l)
}

keep_main_results_for_model_list<-function(l,num=2){
  aics = sort(l$aics)
  selected_models = union(names(aics)[1:num],names(aics)[grepl("base_model",names(aics))])
  l = list(models = l$models[selected_models],aics=l$aics[selected_models])
  model_names = names(l$aics)
  new_l = list()
  for(nn in model_names){
    curr_m = l$models[[nn]]
    coeffs = cbind(curr_m$beta,curr_m$zval,curr_m$pval,curr_m$ci.lb,curr_m$ci.ub)
    colnames(curr_coeffs) = c("beta","zval","pval","lb","ub")
    new_l[[nn]] = list(
      aic_c = l$aics[nn],
      coeffs = coeffs,
      mod_p = try_get_field(curr_m,"QMp"),
      het_p = try_get_field(curr_m,"QEp"),
      sigma2 = try_get_field(curr_m,"sigma2"),
      tau2 = try_get_field(curr_m,"tau2"),
      I2 = try_get_field(curr_m,"I2")
    )
  }
  return(new_l)
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

# Load the input
load("meta_analysis_input.RData")
# Run the analysis
all_meta_analysis_res <- list()
for(nn in names(naive_rep_analysis_results)){
  curr_dataset = datasets[[nn]][naive_rep_analysis_results[[nn]]]
  if(grepl("acute",nn)){
    analysis1 = mclapply(curr_dataset,acute_gdata_metaanalysis,mc.cores = num_cores)
  }
  else{
    analysis1 = mclapply(curr_dataset,longterm_gdata_metaanalysis,mc.cores = num_cores)
  }
  analysis2 = unlist(mclapply(curr_dataset,simple_stouffer_meta_analysis,mc.cores=num_cores))
  all_meta_analysis_res[[nn]] = list(model_selection = analysis1,simple_stouffer = analysis2)
  forest(analysis1[[3]]$selected$model)
}
save(all_meta_analysis_res,naive_rep_analysis_results,file="meta_analysis_results.RData")




