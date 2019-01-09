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
source('/Users/David/Desktop/repos/motrpac/metaanalysis/helper_functions.R')
source('/Users/David/Desktop/repos/motrpac/metaanalysis/helper_functions_meta_analaysis.R')
entrez2symbol = as.list(org.Hs.egSYMBOL)

# Get the datasets and their metadata
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
acute_datasets = cohort_data
acute_metadata = cohort_metadata
load("PADB_univariate_results_and_preprocessed_data_longterm.RData")
longterm_datasets = cohort_data
longterm_metadata = cohort_metadata
load("PADB_dataset_level_meta_analysis_data.RData")

#' Change the time field from numeric to ordered factor.
#' @param gdata A data frame with a column "time"
#' @return A new data frame.
simplify_time_in_gdata<-function(gdata,func=simplify_time_acute){
  gdata$time = func(gdata$time)
  gdata$time = ordered(gdata$time)
  return(gdata)
}
acute_gene_tables_simpletime = lapply(acute_gene_tables,simplify_time_in_gdata)
acute_gene_tables_raw_simpletime = lapply(acute_gene_tables_raw,simplify_time_in_gdata)

# Some simple tests
obj = acute_gene_tables_simpletime[["10891"]]
obj = obj[obj$tissue=="muscle",]
colnames(obj)
# Very simple models without specifying repeated measures
res = rma(yi,vi,data = obj,mods = ~1+training+time+as.numeric(avg_age)+prop_males,weights = 1/vi)
res0 = rma(yi,vi,data = obj,weights=1/vi)
# useful functions: funnel, forest, radial, confint, summary, rstandard, rstudent, residuals
# Add noise structure to random effects
# The nested structure is required to have different correlated
# random effects for time points within a cohort.
res_mv1 = rma.mv(yi,vi,data = obj,mods = ~1+training+time,
                 random = ~ V1|gse,struct="AR")
res_mv2 = rma.mv(yi,vi,data = obj,
                 mods = ~1+training+time+as.numeric(prop_males),
                 random = list(~ V1|gse),struct="AR")
plot(residuals(res_mv2),col=as.factor(obj$gse))
forest(res_mv2)
funnel(res_mv1)
funnel(res_mv2)

# A series of helper functions for implementing our algorithm
controls_meta_analysis<-function(gdata){
  gdata = gdata[!grepl("treatment",gdata$training),]
  gdata_c = gdata[grepl("untrained",gdata$training),]
  beta_c_est = NA; beta_c_lb = NA; beta_c_ub = NA
  if(nrow(gdata_c)>0){
    obj_c = rma(yi,vi,weights = 1/vi,data=gdata_c)
    beta_c_est = obj_c$beta[[1]]
    beta_c_lb = obj_c$ci.lb[1]
    beta_c_ub = obj_c$ci.ub[1]
  }
  return(c(beta_c_est,beta_c_lb,beta_c_ub))
}

#' Implementation of our single window analysis
#' @param gdata A data frame with a single time point in the $time field
#' @param remove_treatment A logical value indicating whether we want to remove all special treatment cohorts (e.g., drug+exercise).
#' @return A vector with the results of the meta-analysis.
time_window_metaanalysis<-function(gdata,remove_treatment=T){
  if(remove_treatment){gdata = gdata[!grepl("treatment",gdata$training),]}
  gdata$vi = pmax(0.00001,gdata$vi)
  gdata_e = gdata[!grepl("untrained",gdata$training),]
  gdata_c = gdata[grepl("untrained",gdata$training),]
  beta_c_est = NA; beta_c_lb = NA; beta_c_ub = NA
  if(nrow(gdata_c)>0){
    obj_c = rma(yi,vi,weights = 1/vi,data=gdata_c)
    beta_c_est = obj_c$beta[[1]]
    beta_c_lb = obj_c$ci.lb[1]
    beta_c_ub = obj_c$ci.ub[1]
  }
  has_tr = length(unique(gdata_e$training))>1
  if(has_tr){
    success = F
    gdata_e = tryCatch({
      obj_e = rma.mv(yi,vi,data = gdata_e,random= ~V1|gse,struct="AR",mods = ~training,
                     control=list(iter.max=10000))
      success = T
    },error = function (e){
          if(grepl("Ratio of largest to smallest sampling variance extremely large",e)){
            gdata_e$vi = pmax(gdata_e$vi,0.001)
          }
          return(gdata_e)}
    )
    gdata_e = gdata[!grepl("untrained",gdata$training),]
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,data = gdata_e,random= ~V1|gse,struct="AR",mods = ~training,
                       control=list(iter.max=10000,rel.tol=1e-8)) 
        success = T
      })
    }
    gdata_e = gdata[!grepl("untrained",gdata$training),]
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,data = gdata_e,random= ~V1|gse,struct="AR",mods = ~training,
                       control=list(iter.max=10000,rel.tol=1e-7))
      })
    }
  }
  else{
    success = F
    gdata_e = tryCatch({
      obj_e = rma.mv(yi,vi,data = gdata_e,random= ~V1|gse,struct="AR",
                     control=list(iter.max=10000))
      success = T
    },error = function (e){
      if(grepl("Ratio of largest to smallest sampling variance extremely large",e)){
        gdata_e$vi = pmax(gdata_e$vi,0.001)
      }
      return(gdata_e)}
    )
    gdata_e = gdata[!grepl("untrained",gdata$training),]
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,data = gdata_e,random= ~V1|gse,struct="AR",
                       control=list(iter.max=10000,rel.tol=1e-8)) 
        success = T
      })
    }
    if(!success){
      obj_e = rma.mv(yi,vi,data = gdata_e,random= ~V1|gse,struct="AR",
                     control=list(iter.max=10000,rel.tol=1e-7)) 
    }
  }
  # The tryCatch may had removed gdata_e, we recreate it here
  gdata_e = gdata[!grepl("untrained",gdata$training),]
  betas_e = obj_e$beta[,1]
  b_e_0 = betas_e[1]
  # SCORE 1
  beta_a = b_e_0
  if(has_tr){
    for(j in 2:length(betas_e)){
      beta_a = c(beta_a,b_e_0+betas_e[j])
    }
  }
  score1_data = c(beta_c_est,beta_c_lb,beta_c_ub,beta_a)
  names(score1_data)[1:3] = c("beta_c_est","beta_c_lb","beta_c_ub")
  names(score1_data)[4:length(score1_data)] = paste("beta_a",1:(length(score1_data)-3),sep="_")
  score1_ME_p = NA
  if(nrow(gdata_c)>0){
    is_ctrl = grepl("untrained",gdata$training)
    gdata = cbind(is_ctrl,gdata)
    try({
      score1_ME_p = rma.mv(yi,vi,mods = ~is_ctrl,random= ~ V1|gse,data=gdata,
                           control=list(iter.max=10000,rel.tol=1e-8))$pval[2]
    })
    if(is.na(score1_ME_p)){
      try({
        score1_ME_p = rma.mv(yi,vi,mods = ~is_ctrl,random= ~ V1|gse,data=gdata,
                             control=list(iter.max=10000,rel.tol=1e-7))$pval[2]
      })
    }
    if(is.na(score1_ME_p)){
      try({
        score1_ME_p = rma.mv(yi,vi,mods = ~is_ctrl,random= ~ V1|gse,data=gdata,
                             control=list(iter.max=10000,rel.tol=1e-7))$pval[2]
      })
    }
  }
  score1_data = c(score1_data,score1_ME_p)
  names(score1_data)[length(score1_data)] = "score1_ME_p"
  # SCORE 2
  gdata_egger = gdata_e
  if(any(table(gdata$gse)>1)){
    gses = gdata_egger$gse
    gdata_egger = gdata_egger[,c("yi","vi")]
    gdata_egger = apply(gdata_egger,2,function(x,y)tapply(x,y,mean),y=gses)
  }
  egger_test_p = NA
  if(!is.null(dim(gdata_egger)) && nrow(gdata_egger)>2){
    try({
      egger_test = regtest(rma(yi,vi,weights = 1/vi,data=gdata_egger,control=list(maxiter=10000,stepadj=0.5)))
      egger_test_p = egger_test$pval
    })
  }
  # Prepare outpus
  significance_e = c(obj_e$pval[1],obj_e$QMp)
  names(significance_e) = c("p_intercept","p_modifiers")
  names(betas_e) = paste("beta_e_",names(betas_e),sep="")
  out = c(betas_e,significance_e,score1_data,egger_test_p)
  names(out)[length(out)] = "egger_test_p"
  return(out)
}
perform_timepoint_metaanalyses<-function(gdata,time_simp_func=NULL,tissue="muscle",min_studies=3){
  gdata = gdata[grepl(tissue,gdata$tissue),]
  if(!is.null(time_simp_func)){
    gdata$time = time_simp_func(gdata$time)
  }
  out = c()
  for(tt in unique(gdata$time)){
    currg = gdata[gdata$time == tt,]
    if(length(unique(currg[currg$training!="untrained","gse"]))<min_studies){next}
    v = time_window_metaanalysis(currg)
    names(v) = paste(paste("tp",tt,sep="_"),names(v),sep=';')
    out = c(out,v)
  }
  return(out)
}
perform_score_1_filter<-function(out,thr=0.5){
  beta_c = out[grepl("beta_c_est",names(out))]
  beta_as = out[grepl("beta_a",names(out))]
  s = abs(beta_as-beta_c)
  return(any(s>=thr))
}
perform_score_1_filter_ci<-function(out){
  beta_c_lbs = out[grepl("beta_c_lb",names(out))]
  beta_c_ubs = out[grepl("beta_c_ub",names(out))]
  beta_as = out[grepl("beta_a",names(out))]
  s1 = beta_as>beta_c_ubs
  s2 = beta_as<beta_c_lbs
  return(any(s1|s2))
}
perform_filter_1_on_tps<-function(out,thr=0.5){
  arrs = sapply(names(out),strsplit,split=';')
  tps = sapply(arrs,function(x)x[1])
  return(sapply(unique(tps),function(x,y,z,thr)perform_score_1_filter(y[z==x],thr=thr),y=out,z=tps,thr=thr))
}
perform_filter_1_on_tps_ci<-function(out){
  arrs = sapply(names(out),strsplit,split=';')
  tps = sapply(arrs,function(x)x[1])
  return(sapply(unique(tps),function(x,y,z)perform_score_1_filter_ci(y[z==x]),y=out,z=tps))
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
get_beta_a_vec<-function(out){
  return(out[grepl("beta_a_",names(out))])
}
get_control_stat_test_p<-function(out){
  return(out[grepl("score1_ME_p",names(out))])
}

# Run the flow
tp_meta_analysis_results = list()
tp_meta_analysis_results[["acute,muscle"]] = t(sapply(acute_gene_tables_raw_simpletime,perform_timepoint_metaanalyses,tissue="muscle"))
tp_meta_analysis_results[["acute,blood"]] = t(sapply(acute_gene_tables_raw_simpletime,perform_timepoint_metaanalyses,tissue="blood"))
tp_meta_analysis_results[["longterm,muscle"]] = t(sapply(longterm_gene_tables_raw,perform_timepoint_metaanalyses,
         time_simp_func=simplify_time_longterm_muscle,tissue="muscle"))
tp_meta_analysis_results[["longterm,blood"]] = t(sapply(longterm_gene_tables_raw,perform_timepoint_metaanalyses,
         time_simp_func=simplify_time_longterm_blood,tissue="blood"))
# Select the large windows
tp_meta_analysis_results_all_windows = tp_meta_analysis_results
tp_meta_analysis_results[["acute,muscle"]] = tp_meta_analysis_results_all_windows[["acute,muscle"]][,grepl("tp_4;",
        colnames(tp_meta_analysis_results[["acute,muscle"]]))] 
tp_meta_analysis_results[["acute,blood"]] = tp_meta_analysis_results_all_windows[["acute,blood"]][,grepl("tp_1;|tp_24;",
        colnames(tp_meta_analysis_results[["acute,blood"]]))] 
sapply(tp_meta_analysis_results,dim)
sapply(tp_meta_analysis_results_all_windows,dim)
save(tp_meta_analysis_results,tp_meta_analysis_results_all_windows,file="tp_meta_analysis_results.RData")

get_tps_from_names<-function(x){unique(sapply(x,function(x)strsplit(x,split=';')[[1]][1]))}
tp_is_gene_selected<-function(out,fdr_thr=0.001,effect_thr=0.5,perform_ctrl_test=T,control_effect_thr=0.5){
  pvals = get_tps_pvalues(out)
  pval_test = pvals<=fdr_thr
  if(all(!pval_test)){return(F)}
  beta_as = get_beta_a_vec(out)
  effects_test = abs(beta_as)>=effect_thr
  if(all(!effects_test)){return(F)}
  curr_tps = intersect(get_tps_from_names(names(which(pval_test))),get_tps_from_names(names(which(effects_test))))
  if(length(curr_tps)==0){return(F)}
  if(perform_ctrl_test){
    filter_1_res = perform_filter_1_on_tps(out,control_effect_thr)
    if(all(is.na(filter_1_res))){return(F)}
    filter_1_res = abs(filter_1_res[!is.na(filter_1_res)])
    curr_tps = intersect(curr_tps,get_tps_from_names(names(which(filter_1_res>0))))
  }
  return(length(curr_tps)>0)
}
get_selected_tps_pvals_and_effects<-function(out,fdr_thr=0.001,effect_thr=0.5,perform_ctrl_test=T,control_effect_thr=0.5){
  pvals = get_tps_pvalues(out)
  effects = get_tps_betas(out)
  c_effects = out[grepl("beta_c_est",names(out))]
  tps = get_tps_from_names(names(out))
  meta_analysis_summary = c()
  for(tp in tps){
    effs = effects[grepl(names(effects),pattern=tp)]
    interc = effs[grepl(names(effs),pattern="intrc")]
    effs[!grepl(names(effs),pattern="intrc")] = effs[!grepl(names(effs),pattern="intrc")] + interc
    names(effs)[!grepl(names(effs),pattern="intrc")] = gsub( names(effs)[!grepl(names(effs),pattern="intrc")],pattern="beta_e_training",replace='')
    names(effs)[grepl(names(effs),pattern="intrc")] = "base"
    mean_control = c_effects[grepl(names(c_effects),pattern=tp)]
    v = c(effs,mean_control)
    names(v)[length(v)] = paste(tp,";untrained_fc",sep='')
    meta_analysis_summary = c(meta_analysis_summary,v)
  }
  return(meta_analysis_summary)
}

# Analysis of the results: p-values
pval_matrices = lapply(tp_meta_analysis_results,function(x)t(apply(x,1,get_tps_pvalues)))
sapply(pval_matrices,colnames)
par(mfrow=c(2,2))
for(nn in names(pval_matrices)){
  hist(c(pval_matrices[[nn]]),main=nn,xlab="P-value")
}
par(mfrow=c(4,4))
tmp = sapply(pval_matrices,function(x)apply(x,2,hist))
all_ps = c(unlist(pval_matrices))
all_ps = all_ps[!is.na(all_ps)]
par(mfrow=c(1,1));hist(all_ps,main="All p-values", xlab="P-value")
all_ps_thr = max(all_ps[p.adjust(all_ps,method="fdr")<=0.05])
print(length(all_ps));print(all_ps_thr)
pval_matrices_binary = lapply(pval_matrices,function(x,y)x<=y,y=all_ps_thr)
sapply(pval_matrices_binary,table)
m = pval_matrices_binary[[1]][,1:2]
colSums(m)
table(m[,1],m[,2])

# # vs. controls significance - 
# # NUMBER OF UNTRAINED COHORTS SEEMS TO BE TOO LOW
# ctrl_test_ps = lapply(tp_meta_analysis_results,function(x)t(apply(x,1,get_control_stat_test_p)))
# for(i in 1:length(ctrl_test_ps)){
#   if(nrow(ctrl_test_ps[[i]])==1){ctrl_test_ps[[i]] = t(ctrl_test_ps[[i]])}
# }
# all_ps = c(unlist(ctrl_test_ps))
# all_ps = all_ps[!is.na(all_ps)]
# par(mfrow=c(1,1));hist(all_ps,main="All p-values", xlab="P-value")
# all_ps_thr = max(all_ps[p.adjust(all_ps)<=0.1])
# print(length(all_ps));print(all_ps_thr)
# sapply(ctrl_test_ps,function(x)colSums(x<=all_ps_thr,na.rm = T))
# sapply(ctrl_test_ps,function(x)x["10891",])
# sapply(ctrl_test_ps,function(x)x["7422",])
# sapply(ctrl_test_ps,function(x)x["3791",])
# sapply(ctrl_test_ps,function(x)x["4521",])
# sapply(ctrl_test_ps,function(x)x["4619",])
# gdata = longterm_gene_tables[["4619"]]
# gdata = gdata[gdata$tissue=="muscle",]
# tp_meta_analysis_results[[3]]["4619",]

# gene selection based on all tests
# In the analyses below 0.001 threshold for the p-value is similar to using 0.1 BY FDR
gene_selection_all_tests = list()
gene_selection_all_tests[["acute,muscle"]]  = apply(tp_meta_analysis_results[["acute,muscle"]],1,tp_is_gene_selected,
      fdr_thr=0.001,effect_thr=0.5,perform_ctrl_test=T,control_effect_thr=0.5)
gene_selection_all_tests[["acute,blood"]]  = apply(tp_meta_analysis_results[["acute,blood"]],1,tp_is_gene_selected,
      fdr_thr=0.001,effect_thr=0.5,perform_ctrl_test=F)
gene_selection_all_tests[["longterm,muscle"]]  = apply(tp_meta_analysis_results[["longterm,muscle"]],1,tp_is_gene_selected,
      fdr_thr=0.001,effect_thr=0.25,perform_ctrl_test=T,control_effect_thr=0.25)
gene_selection_all_tests[["longterm,blood"]]  = apply(tp_meta_analysis_results[["longterm,blood"]],1,tp_is_gene_selected,
      fdr_thr=0.001,effect_thr=0.25,perform_ctrl_test=F)
sapply(gene_selection_all_tests,sum)
selected_genes_all_tests = sapply(gene_selection_all_tests,function(x)names(which(x)))
sapply(selected_genes_all_tests,length)
selected_genes_all_tests = selected_genes_all_tests[sapply(selected_genes_all_tests,length)>0]

save(tp_meta_analysis_results,tp_meta_analysis_results_all_windows,gene_selection_all_tests,
     selected_genes_all_tests,file="tp_meta_analysis_results.RData")

# Wrtie the results into a file
get_gene_summary_table<-function(l,e2s=entrez2symbol){
  all_genes = unique(unlist(l))
  mat = matrix(F,nrow=length(all_genes),ncol=length(l))
  colnames(mat) = names(l)
  rownames(mat) = all_genes
  for(nn in names(l)){
    mat[l[[nn]],nn] = T
  }
  sum_col = apply(mat,1,function(x,y)paste(y[x],collapse=" AND "),y=colnames(mat))
  table(sum_col)
  mat = cbind(sum_col,mat)
  colnames(mat)[1] = "DiscoveredIn"
  mat = cbind(e2s[all_genes],mat)
  colnames(mat)[1] = "Symbol"
  mat = cbind(rownames(mat),mat)
  colnames(mat)[1] = "Entrez"
  rownames(mat) = NULL
  return(mat)
}
get_gene_summary_table(selected_genes_all_tests)
write.table(get_gene_summary_table(selected_genes_all_tests),file="metaanalysis_genes.txt",
            sep="\t",col.names = T,row.names = F,quote=F)

library('Vennerable')
V = Venn(selected_genes_all_tests)
plot(V,doWeights=F)

selected_genes_all_tests_names = lapply(selected_genes_all_tests,function(x,y)unlist(y[x]),y=entrez2symbol)
known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
                "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")
sapply(selected_genes_all_tests_names,intersect,y=known_genes)
selected_genes_all_tests_names[[1]][grepl("^COL",selected_genes_all_tests_names[[1]])]
selected_genes_all_tests_names[[3]][grepl("^COL",selected_genes_all_tests_names[[3]])]

selected_genes_all_tests_topgo = run_topgo_enrichment_fisher(selected_genes_all_tests,
                                union(names(acute_gene_tables),names(longterm_gene_tables)))
table(extract_top_go_results(selected_genes_all_tests_topgo)[,1])
get_most_sig_enrichments_by_groups(extract_top_go_results(selected_genes_all_tests_topgo,0.1),2)

save(tp_meta_analysis_results,tp_meta_analysis_results_all_windows,gene_selection_all_tests,
     selected_genes_all_tests_topgo,selected_genes_all_tests_names,known_genes,
     selected_genes_all_tests,file="tp_meta_analysis_results.RData")


#######################################################
#######################################################
# Look at publication bias
#######################################################
#######################################################
egger_test_results = lapply(tp_meta_analysis_results,function(x)apply(x,1,get_tps_egger_tests))
par(mfrow=c(2,2))
all_egger_test_ps = unlist(sapply(egger_test_results[1:2],c))
all_egger_test_ps = all_egger_test_ps[!is.na(all_egger_test_ps)]
print(length(all_egger_test_ps))
hist(all_egger_test_ps,main="Egger test, acute studies",xlab="P-value")
egger_lfdrs = get_lfdrs(all_egger_test_ps)

all_egger_test_ps = unlist(sapply(egger_test_results[3:4],c))
print(length(all_egger_test_ps))
all_egger_test_ps = all_egger_test_ps[!is.na(all_egger_test_ps)]
hist(all_egger_test_ps,main="Egger test, long-term studies",xlab="P-value")
egger_lfdrs = get_lfdrs(all_egger_test_ps)

#######
par(mfrow=c(2,2))
for(nn in names(egger_test_results)){
  if(!is.null(dim(egger_test_results[[nn]]))){
    selected_genes_test_ps = c(egger_test_results[[nn]][,pvals_filter1_selected_genes[[nn]]])
  }
  else{
    selected_genes_test_ps = c(egger_test_results[[nn]][pvals_filter1_selected_genes[[nn]]])
  }
  selected_genes_test_ps = selected_genes_test_ps[!is.na(selected_genes_test_ps)]
  all_ps = c(egger_test_results[[nn]])
  all_ps = all_ps[!is.na(all_ps)]
  print(wilcox.test(qnorm(selected_genes_test_ps),qnorm(all_ps)));abline(0,1)
  print(table(selected_genes_test_ps<0.01)/length(selected_genes_test_ps))
  print(table(all_ps < 0.01)/length(all_ps))
  print(length(all_egger_test_ps))
}

save(tp_meta_analysis_results,tp_meta_analysis_results_all_windows,gene_selection_all_tests,
     selected_genes_all_tests_topgo,selected_genes_all_tests_names,known_genes,egger_test_results,
     selected_genes_all_tests,file="tp_meta_analysis_results.RData")

# Dec 2017: reformat the main results and extract the summary matrices for the selected genes
sapply(tp_meta_analysis_results,colnames)
tp_meta_analysis_results_statistics = list()
for(nn in names(selected_genes_all_tests)){
  mm = tp_meta_analysis_results[[nn]][selected_genes_all_tests[[nn]],]
  rownames(mm) = entrez2symbol[rownames(mm)]
  tp_meta_analysis_results_statistics[[nn]] = t(apply(mm,1,get_selected_tps_pvals_and_effects))
}
sapply(tp_meta_analysis_results_statistics,colnames)
save(tp_meta_analysis_results_statistics,file = "tp_meta_analysis_results_statistics.RData")

# # Some tests
# gdata = acute_gene_tables_raw_simpletime[["10891"]]
# gdata$vi = pmax(gdata$vi,0.001)
# gdata = gdata[gdata$time==4&gdata$tissue=="muscle"&gdata$training!="untrained",]
# res = rma.mv(yi,vi,data=gdata,random=~1|gse,mods = ~ training)$pval
# print(res)
# res = rma.mv(yi,vi,data=gdata,random=~V1/gse,mods = ~ training)$pval
# print(res)
# gdata = longterm_gene_tables_raw[[1]]
# gdata$vi = pmax(gdata$vi,0.001)
# gdata$time = simplify_time_longterm_muscle(gdata$time)
# gdata = gdata[gdata$tissue=="muscle"&gdata$training!="untrained",]
# gdata = gdata[gdata$time<=100&gdata$tissue=="muscle"&gdata$training!="untrained",]
# table(gdata$gse,gdata$time)
# res = rma.mv(yi,vi,data=gdata,random=~1|gse,mods = ~ training)$pval
# print(res)
# res = rma.mv(yi,vi,data=gdata,random=~V1/gse,mods = ~ training)$pval
# print(res)
# # compare samples
# samp = sample(1:10000)[1:300]
# scores1=c();scores2=c()
# for(nn in samp){
#   gdata = acute_gene_tables_raw[[nn]]
#   gdata$vi = pmax(gdata$vi,0.0001)
#   gdata$time = simplify_time_acute(gdata$time)
#   gdata$time = 100
#   gdata = gdata[gdata$time==100&gdata$tissue=="muscle"&gdata$training!="untrained",]
#   try({
#     res1 = rma.mv(yi,vi,data=gdata,random=~1|gse,mods = ~ training)$pval
#     res2 = rma.mv(yi,vi,data=gdata,random=~V1/gse,struct="HCS",mods = ~ training)$pval
#     scores1=rbind(scores1,res1)
#     scores2=rbind(scores2,res2)
#   })
# }
# plot(c(scores1),c(scores2));abline(0,1)
# hist(scores1);hist(scores2)
# 
# # Look at the number of gses per expected analysis group
# gdata = longterm_gene_tables_raw[[1]]
# length(unique(gdata[gdata$tissue=="muscle"&gdata$training!="untrained","gse"]))
# length(unique(gdata[gdata$tissue=="blood"&gdata$training!="untrained","gse"]))
# gdata = acute_gene_tables_raw[[1]]
# length(unique(gdata[gdata$tissue=="muscle"&gdata$training!="untrained","gse"]))
# length(unique(gdata[gdata$tissue=="blood"&gdata$training!="untrained","gse"]))
# get_num_gses_per_simplified_time_func<-function(gdata,tissue="muscle",func = function(x)x){
#   gdata$time = func(gdata$time)
#   gdata = gdata[gdata$tissue==tissue,]
#   gdata = gdata[!grepl("untrained",gdata$training),]
#   gdata = gdata[!grepl("treatment",gdata$training),]
#   return(table(gdata$gse,gdata$time))
# }
# gdata = longterm_gene_tables_raw[[1]]
# gdata = gdata[!grepl("treatment",gdata$training),]
# get_num_gses_per_simplified_time_func(gdata,"muscle")
# get_num_gses_per_simplified_time_func(gdata,"blood")
# gdata = acute_gene_tables_raw_simpletime[[1]]
# gdata = acute_gene_tables_raw[[1]]
# gdata = gdata[!grepl("treatment",gdata$training),]
# get_num_gses_per_simplified_time_func(gdata,"muscle")
# get_num_gses_per_simplified_time_func(gdata,"blood")
# gdata = acute_gene_tables_raw[[1]]
# simplify_time_acute_simple <-function(tt){
#   tt[tt<10] = 4; tt[tt>=10]=24
#   return(tt)
# }
# get_num_gses_per_simplified_time_func(gdata,"muscle",simplify_time_acute_simple)
# get_num_gses_per_simplified_time_func(gdata,"blood",simplify_time_acute_simple)

# tp_meta_analysis_results[["acute,muscle"]] = list()
# for(g in setdiff(names(acute_gene_tables_raw_simpletime),names(tp_meta_analysis_results[["acute,muscle"]]))){
#   tp_meta_analysis_results[["acute,muscle"]][[g]] = perform_timepoint_metaanalyses(acute_gene_tables_raw_simpletime[[g]],tissue="muscle")
# }
# tp_meta_analysis_results[["longterm,muscle"]] = list()
# for(g in setdiff(names(longterm_gene_tables_raw_simpletime),names(tp_meta_analysis_results[["longterm,muscle"]]))){
#   tp_meta_analysis_results[["longterm,muscle"]][[g]] = perform_timepoint_metaanalyses(longterm_gene_tables_raw[[g]],
#                                                        time_simp_func=simplify_time_longterm_muscle,tissue="muscle")
# }

# # Manual examinations
# pval_matrices[[1]]["5166",]
# pass_test_matrices[[1]]["5166",]
# control_filter_tests[[1]]["5166",]
# tp_meta_analysis_results[[1]]["7422",]
# gdata = acute_gene_tables_raw[["7422"]]
# gdata = longterm_gene_tables_raw[["5166"]]
# gdata = acute_gene_tables_raw_simpletime[["10891"]]
# gdata = acute_gene_tables_raw[["7139"]]
# rma.mv(yi,vi,mods = ~ tissue:factor(time) + training, random = list(~1|V1,~1|gse),data=gdata,btt=2:4)
# rma.mv(yi,vi,mods = ~ training + factor(time), 
#        random = list(~1|V1,~1|gse),data=gdata,subset = tissue=="muscle")
# gdata1 = acute_gene_tables_raw_simpletime[["10891"]]
# gdata2 = longterm_gene_tables[["10891"]]
# gdata = acute_gene_tables_raw_simpletime[["7070"]]
# gdata = longterm_gene_tables_raw[["7070"]]
# gdata = gdata[gdata$tissue=="muscle",]
# gdata$time = simplify_time_longterm_muscle(gdata$time)
# gdata$time = simplify_time_acute(gdata$time)
# out = perform_timepoint_metaanalyses(gdata,tissue="muscle",time_simp_func=simplify_time_longterm_muscle)
# out = perform_timepoint_metaanalyses(gdata,tissue="muscle")
# out = perform_timepoint_metaanalyses(gdata,tissue="blood",time_simp_func=simplify_time_longterm_blood)
# perform_filter_1_on_tps(out,0.5)
# get_tps_pvalues(out) < all_ps_thr
# get_tps_egger_tests(out)
# get_tps_betas(out)
# # gdata = acute_gene_tables_raw_simpletime[["10891"]]
# # gdata = acute_gene_tables_raw_simpletime[["10161"]]
# # gdata = gdata[gdata$tissue=="blood",]
# # out = perform_timepoint_metaanalyses(gdata,tissue="muscle")
# # out = perform_timepoint_metaanalyses(gdata,tissue="blood")
# # out = perform_timepoint_metaanalyses(gdata,tissue="blood|muscle")
# # perform_filter_1_on_tps(out)
# # get_tps_pvalues(out)
# # get_tps_egger_tests(out)
# # get_tps_betas(out)
# # check adding controls
# is_ctrl = grepl("control",gdata$training)
# gdata = cbind(is_ctrl,gdata)
# rma.mv(yi,vi,mods = ~is_ctrl,random= ~ V1|gse,data=gdata)
# gdata = acute_gene_tables_raw_simpletime[["10891"]]
# gdata = gdata[gdata$time==4,]
# gdata = gdata[gdata$training!="control",]
# get_subset_forest_plot(gdata,"muscle",main="PGC1, acute, muscle, 4h")








