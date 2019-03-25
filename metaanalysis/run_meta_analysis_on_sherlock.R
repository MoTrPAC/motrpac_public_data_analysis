
library(parallel)
library(metafor,lib.loc="~/R/packages")
library(nloptr,lib.loc="~/R/packages")

print("Usage:<working dir with the input RData>
      <num cores><indices of datasets in meta_reg_datasets, comma seperated single string>
      <gene start index><gene end index>")

args = commandArgs(trailingOnly=TRUE)
print("Input args are:")
print(args)

if(length(args)!=5){
  print("Number of command line arguments should be 3, check the command and rerun")
  q("no")
}

num_cores = as.numeric(args[2])
inds = as.numeric(strsplit(args[3],split=",")[[1]])
setwd(args[1])
start = as.numeric(args[4])
end = as.numeric(args[5])

############################################################################
############################################################################
############################################################################
# TODO: merge with the functions in simplified_moderators_metaanalysis.R
# TODO: we currently have code dup
# Functions for the meta-analysis
# ... specifies random effects formula and structure
model_selection_meta_analysis<-function(gdata,
    mod_names = c("training","time","avg_age","prop_males"),...){
  to_rem = sapply(mod_names,function(x,y)length(unique(y[,x],na.rm=T))<2,y=gdata)
  mod_names = mod_names[!to_rem] 
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
  return(list(models=models,aics=aics,mod_names=mod_names))
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
  res=NULL
  try({
    res = func(yi,vi,data=gdata,...)
  },silent=TRUE)
  if(!is.null(res)){return(res)}
  
  for(rel.tol in c(1e-8,1e-7,1e-6)){
    cc=list(iter.max=10000,rel.tol=rel.tol)
    res = NULL
    try({
      res = func(yi,vi,data=gdata,control=cc,...)
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

# acute_gdata_metaanalysis<-function(gdata,permtest=F,
#     mod_names = c("training","time","avg_age","prop_males")){
#   res1 = model_selection_meta_analysis(gdata,
#        random=list(~ V1|gse),struct="CS",mod_names=mod_names)
#   res0 = model_selection_meta_analysis(gdata,func=rma.uni,mod_names=mod_names)
#   l = list(
#     models = c(add_prefix_to_names("simple",res0$models),
#                add_prefix_to_names("base2",res1$models)),
#     aics = c(add_prefix_to_names("simple",res0$aics),
#              add_prefix_to_names("base2",res1$aics))
#   )
#   # keep the simple models and the top 2
#   l = keep_main_results_for_model_list(l)
#   gc()
#   return(l)
# }

gdata_metaanalysis<-function(gdata,permtest=F,
     mod_names = c("training","time","avg_age","prop_males")){
  res1 = model_selection_meta_analysis(gdata,random=list(~ V1|gse),struct="CS",mod_names=mod_names)
  res0 = model_selection_meta_analysis(gdata,func=rma.uni,mod_names=mod_names)
  if(length(res1$models)==0 && length(res0$models)==0){return(NULL)}
  if(length(res1$models)>0){
    l = list(
      models = c(add_prefix_to_names("simple",res0$models),
                 add_prefix_to_names("base2",res1$models)),
      aics = c(add_prefix_to_names("simple",res0$aics),
               add_prefix_to_names("base2",res1$aics))
    )
  }
  else{
    l = list(
      models = add_prefix_to_names("simple",res0$models),
      aics = add_prefix_to_names("simple",res0$aics)
    )
  }
  # keep the simple models and the top 2
  l = keep_main_results_for_model_list(l)
  gc()
  return(l)
}

# longterm_gdata_metaanalysis<-function(gdata,permtest=F,simple_output=T,
#     mod_names = c("training","time","avg_age","prop_males")){
#   res1 = model_selection_meta_analysis(gdata,random=list(~ 1|gse),mod_names=mod_names)
#   res0 = model_selection_meta_analysis(gdata,func=rma.uni,mod_names=mod_names)
#   l = list(
#     models = c(add_prefix_to_names("simple",res0$models),
#                add_prefix_to_names("base2",res1$models)),
#     aics = c(add_prefix_to_names("simple",res0$aics),
#              add_prefix_to_names("base2",res1$aics))
#   )
#   # keep the simple models and the top 2
#   l = keep_main_results_for_model_list(l)
#   gc()
#   return(l)
# }

keep_main_results_for_model_list<-function(l,num=2){
  aics = sort(l$aics)
  selected_models = union(names(aics)[1:num],names(aics)[grepl("base_model",names(aics))])
  l = list(models = l$models[selected_models],aics=l$aics[selected_models])
  model_names = names(l$aics)
  new_l = list()
  for(nn in model_names){
    if(is.na(nn)){next}
    curr_m = l$models[[nn]]
    coeffs = cbind(curr_m$beta,curr_m$se,curr_m$zval,curr_m$pval,curr_m$ci.lb,curr_m$ci.ub)
    colnames(coeffs) = c("beta","se","zval","pval","lb","ub")
    new_l[[nn]] = list(
      aic_c = l$aics[nn],
      coeffs = coeffs,
      mod_p = try_get_field(curr_m,"QMp"),
      het_p = try_get_field(curr_m,"QEp"),
      het_q = try_get_field(curr_m,"QE"),
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
print("Gene tables were loaded, the meta-analyses names are:")
print(names(meta_reg_datasets))
meta_reg_datasets = meta_reg_datasets[inds]
print("restricting the analysis to the specified meta analysis type, dataset is:")
print(names(meta_reg_datasets)[1])
nn = names(meta_reg_datasets)[1]
start = max(0,start)
end = min(end,length(meta_reg_datasets[[nn]]))
print(paste("running on genes from",start,"to",end))
curr_dataset = meta_reg_datasets[[nn]][start:end]
meta_reg_datasets = NULL
gc()
print("Moderators to be tested for each gene:")
curr_mods = meta_reg_to_mods[[nn]]
print(curr_mods)
analysis_res = list()
for(gg in names(curr_dataset)){
  print(paste("analyzing gene:",gg))
  analysis_res[[gg]] = gdata_metaanalysis(curr_dataset[[gg]],mod_names=curr_mods)
}
out_file = paste("meta_analysis_results",nn,start,end,".RData",sep="_")
out_file = gsub(",","_",out_file)
save(analysis_res,file=out_file)




