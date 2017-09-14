setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library(metafor)
library(org.Hs.eg.db)
source('repos/motrpac/helper_functions.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata

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

simplify_time_in_gdata<-function(gdata,func=simplify_time_acute){
  gdata$time = func(gdata$time)
  gdata$time = ordered(gdata$time)
  return(gdata)
}
acute_gene_tables_simpletime = lapply(acute_gene_tables,simplify_time_in_gdata)
acute_gene_tables_raw_simpletime = lapply(acute_gene_tables_raw,simplify_time_in_gdata)
longterm_gene_tables_simpletime = lapply(longterm_gene_tables,simplify_time_in_gdata,func=simplify_time_longterm)
longterm_gene_tables_raw_simpletime = lapply(longterm_gene_tables_raw,simplify_time_in_gdata,func=simplify_time_longterm)

# Our algorithm:
# Input: datasets for time point t and gene g
# 1. Run meta-analysis for endurence and resistance (and both)
# 2. Run meta-analysis for controls
# 3. Get controls' intercept estimation beta_c
# 4. Get current gene's overall diff expression: beta_e = max(beta_0,beta_0 + beta_endurence)
# 5. SCORE 1: abs(beta_e-beta_c)
# 6. SCORE 2: Egger's test of the exercise meta-analysis
# 7. Return SCORE 1, SCORE 2, and the betas (and their significance) from the exercise meta-analysis

# Input assumption: gdata has a single time point in the $time field
time_point_metaanalysis<-function(gdata,remove_other=T){
  if(remove_other){
    gdata = gdata[gdata$training!="other",]
  }
  gdata = gdata[gdata$vi>0,]
  gdata_e = gdata[gdata$training!="control",]
  gdata_c = gdata[gdata$training=="control",]
  beta_c = NA
  if(nrow(gdata_c)>0){
    obj_c = rma(yi,vi,weights = 1/vi,data=gdata_c)
    beta_c = obj_c$beta[[1]]
  }
  has_tr = length(unique(gdata_e$training))>1
  if(has_tr){
    success = F
    try({
      obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,mods = ~training,
                     control=list(iter.max=10000))
      success = T
    })
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,mods = ~training,
                       control=list(iter.max=10000,rel.tol=1e-8)) 
        success = T
      })
    }
    if(!success){
      obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,mods = ~training,
                     control=list(iter.max=10000,rel.tol=1e-7)) 
    }
  }
  else{
    success = F
    try({
      obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,
                     control=list(iter.max=10000))
      success = T
    })
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,
                       control=list(iter.max=10000,rel.tol=1e-8)) 
        success = T
      })
    }
    if(!success){
      obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,
                     control=list(iter.max=10000,rel.tol=1e-7)) 
    }
  }
  betas_e = obj_e$beta[,1]
  b_e_0 = betas_e[1]
  # SCORE 1
  beta_a = b_e_0
  if(has_tr){
    for(j in 2:length(betas_e)){
      beta_a = c(beta_a,b_e_0+betas_e[j])
    }
  }
  score1_data = c(beta_c,beta_a)
  names(score1_data)[1] = "beta_c"
  names(score1_data)[2:length(score1_data)] = paste("beta_a",2:length(score1_data),sep="_")
  # SCORE 2
  gdata_egger = gdata_e
  if(any(table(gdata$gse)>1)){
    gses = gdata_egger$gse
    gdata_egger = gdata_egger[,c("yi","vi")]
    gdata_egger = apply(gdata_egger,2,function(x,y)tapply(x,y,mean),y=gses)
  }
  egger_test_p = NA
  if(!is.null(dim(gdata_egger)) && nrow(gdata_egger)>2){
    egger_test = regtest(rma(yi,vi,weights = 1/vi,data=gdata_egger,control=list(maxiter=10000,stepadj=0.5)))
    egger_test_p = egger_test$pval
  }
  # Prepare outpus
  significance_e = c(obj_e$pval[1],obj_e$QMp)
  names(significance_e) = c("p_intercept","p_modifiers")
  names(betas_e) = paste("beta_e_",names(betas_e),sep="")
  out = c(betas_e,significance_e,score1_data,egger_test_p)
  names(out)[length(out)] = "egger_test_p"
  return(out)
}
perform_timepoint_metaanalyses<-function(gdata,time_simp_func=NULL,tissue="muscle"){
  gdata = gdata[grepl(tissue,gdata$tissue),]
  if(!is.null(time_simp_func)){
    gdata$time = time_simp_func(gdata$time)
  }
  out = c()
  for(tt in unique(gdata$time)){
    currg = gdata[gdata$time == tt,]
    v = time_point_metaanalysis(currg)
    names(v) = paste(paste("tp",tt,sep="_"),names(v),sep=';')
    out = c(out,v)
  }
  return(out)
}
perform_score_1_filter<-function(out,thr=0.5){
  beta_c = out[grepl("beta_c",names(out))]
  beta_as = out[grepl("beta_a",names(out))]
  s = abs(beta_as-beta_c)
  return(any(s>=thr))
}
perform_filter_1_on_tps<-function(out,thr=0.5){
  arrs = sapply(names(out),strsplit,split=';')
  tps = sapply(arrs,function(x)x[1])
  return(sapply(unique(tps),function(x,y,z)perform_score_1_filter(y[z==x]),y=out,z=tps))
}
get_tps_pvalues<-function(out){
  return(out[grepl(";p_",names(out))])
}
get_tps_egger_tests<-function(out){
  return(out[grepl("egger",names(out))])
}
get_tps_betas<-function(out){
  return(out[grepl("beta_e_",names(out))])
}

tp_meta_analysis_results = list()

tp_meta_analysis_results[["longterm,muscle"]] = list()
# for(g in setdiff(names(acute_gene_tables_raw_simpletime),names(tp_meta_analysis_results[["longterm,muscle"]]))){
#   tp_meta_analysis_results[["longterm,muscle"]][[g]] = perform_timepoint_metaanalyses(longterm_gene_tables_raw[[g]],
#                                                        time_simp_func=simplify_time_longterm_muscle,tissue="muscle")
# }
tp_meta_analysis_results[["acute,muscle"]] = t(sapply(acute_gene_tables_raw_simpletime,perform_timepoint_metaanalyses,tissue="muscle"))
tp_meta_analysis_results[["acute,blood"]] = t(sapply(acute_gene_tables_raw_simpletime,perform_timepoint_metaanalyses,tissue="blood"))
tp_meta_analysis_results[["longterm,muscle"]] = t(sapply(longterm_gene_tables_raw,perform_timepoint_metaanalyses,
         time_simp_func=simplify_time_longterm_muscle,tissue="muscle"))
tp_meta_analysis_results[["longterm,blood"]] = t(sapply(longterm_gene_tables_raw,perform_timepoint_metaanalyses,
         time_simp_func=simplify_time_longterm_blood,tissue="blood"))
save(tp_meta_analysis_results,file="tp_meta_analysis_results.RData")

# Analysis of the results
pval_matrices = lapply(tp_meta_analysis_results,function(x)t(apply(x,1,get_tps_pvalues)))
sapply(pval_matrices,colnames)
par(mfrow=c(1,2))
sapply(pval_matrices,hist)
all_ps = c(unlist(pval_matrices))
all_ps = all_ps[!is.na(all_ps)]
all_ps_thr = max(all_ps[p.adjust(all_ps,method='fdr')<0.1])
print(length(all_ps))
print(all_ps_thr)
pval_matrices_binary = lapply(pval_matrices,function(x,y)x<=y,y=all_ps_thr)
pval_selected_genes = lapply(pval_matrices_binary,function(x)rownames(x)[rowSums(x)>0])

control_filter_tests = lapply(tp_meta_analysis_results,function(x)t(apply(x,1,perform_filter_1_on_tps)))
non_all_false_filter1 = lapply(control_filter_tests,function(x)rownames(x)[rowSums(x)>0])

selected_genes_acute_muscle = unlist(entrez2symbol[intersect(non_all_false_filter1[[3]],pval_selected_genes[[1]])])
known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
                "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")
intersect(selected_genes_acute_muscle,known_genes)
selected_genes_acute_muscle[grepl("^COL",selected_genes_acute_muscle)]

# gdata = longterm_gene_tables_raw[["65121"]]
# #gdata = gdata[gdata$tissue=="muscle",]
# gdata$time = simplify_time_longterm_muscle(gdata$time)
# out = perform_timepoint_metaanalyses(gdata,tissue="muscle",time_simp_func=simplify_time_longterm_muscle)
# out = perform_timepoint_metaanalyses(gdata,tissue="blood",time_simp_func=simplify_time_longterm_blood)
# perform_filter_1_on_tps(out)
# get_tps_pvalues(out)
# get_tps_egger_tests(out)
# get_tps_betas(out)
# 
# gdata = acute_gene_tables_raw_simpletime[["10891"]]
# gdata = acute_gene_tables_raw_simpletime[["10161"]]
# gdata = gdata[gdata$tissue=="blood",]
# out = perform_timepoint_metaanalyses(gdata,tissue="muscle")
# out = perform_timepoint_metaanalyses(gdata,tissue="blood")
# out = perform_timepoint_metaanalyses(gdata,tissue="blood|muscle")
# perform_filter_1_on_tps(out)
# get_tps_pvalues(out)
# get_tps_egger_tests(out)
# get_tps_betas(out)

















