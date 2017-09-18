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

# ####### Look at datasets with extremely low variances ######
# variances = list();ps = list()
# for(nn1 in names(acute_datasets)){
#   l = acute_datasets[[nn1]]$time2ttest_stats
#   variances[[nn1]] = sapply(l,function(x)x[,"vi"])
#   ps[[nn1]] = sapply(l,function(x)x[,"p"])
# }
# sapply(variances,dim)
# variances = variances[sapply(variances,function(x)!is.null(dim(x)))]
# ps = ps[names(variances)]
# low_vars_prop = sapply(variances,function(x){x=x*x;table(c(x)<0.001)["TRUE"]/length(x)})
# low_vars_prop[is.na(low_vars_prop)]=0
# names(low_vars_prop) = names(variances)
# cohort2size = sapply(acute_metadata,function(x)length(x$gsms))
# cohort2size = cohort2size[names(variances)]
# cohort2num_times = table(acute_gene_tables_raw[[1]]$V1)[names(variances)]
# cor(low_vars_prop,cohort2size,method='spearman')
# cor(low_vars_prop,cohort2size/cohort2num_times,method='spearman')
# plot(low_vars_prop,cohort2size/cohort2num_times)
# ####################
simplify_time_in_gdata<-function(gdata,func=simplify_time_acute){
  gdata$time = func(gdata$time)
  gdata$time = ordered(gdata$time)
  return(gdata)
}
acute_gene_tables_simpletime = lapply(acute_gene_tables,simplify_time_in_gdata)
acute_gene_tables_raw_simpletime = lapply(acute_gene_tables_raw,simplify_time_in_gdata)
longterm_gene_tables_simpletime = lapply(longterm_gene_tables,simplify_time_in_gdata,func=simplify_time_longterm)
longterm_gene_tables_raw_simpletime = lapply(longterm_gene_tables_raw,simplify_time_in_gdata,func=simplify_time_longterm)

# For display items
gd = acute_gene_tables_raw[[1]]
tb = table(gd$time,gd$tissue)
gd$time = simplify_time_acute(gd$time)
tb = table(gd$time,gd$tissue)
gd = longterm_gene_tables_raw[[1]]
tb = table(gd$time,gd$tissue)

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
      obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,mods = ~training,
                     control=list(iter.max=10000))
      success = T
    },error = function (e){
          if(grepl("Ratio of largest to smallest sampling variance extremely large",e)){
            gdata_e$vi = pmax(gdata_e$vi,0.001)
          }
          return(gdata_e)}
    )
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,mods = ~training,
                       control=list(iter.max=10000,rel.tol=1e-8)) 
        success = T
      })
    }
    if(!success){
      try({
        obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,mods = ~training,
                       control=list(iter.max=10000,rel.tol=1e-7))
      })
    }
  }
  else{
    success = F
    gdata_e = tryCatch({
      obj_e = rma.mv(yi,vi,1/vi,data = gdata_e,random= ~1|gse,
                     control=list(iter.max=10000))
      success = T
    },error = function (e){
      if(grepl("Ratio of largest to smallest sampling variance extremely large",e)){
        gdata_e$vi = pmax(gdata_e$vi,0.001)
      }
      return(gdata_e)}
    )
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
  # The tryCatch may had removed gdata_e, we recreate it here
  gdata_e = gdata[gdata$training!="control",]
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
  names(score1_data)[4:length(score1_data)] = paste("beta_a",4:length(score1_data),sep="_")
  score1_ME_p = NA
  if(nrow(gdata_c)>0){
    is_ctrl = grepl("control",gdata$training)
    gdata = cbind(is_ctrl,gdata)
    score1_ME_p = rma.mv(yi,vi,1/vi,mods = ~is_ctrl,random= ~ 1|gse,data=gdata)$pval[2]
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

tp_meta_analysis_results = list()
# tp_meta_analysis_results[["longterm,muscle"]] = list()
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
par(mfrow=c(2,2))
for(nn in names(pval_matrices)){
  hist(c(pval_matrices[[nn]]),main=nn,xlab="P-value")
}
all_ps = c(unlist(pval_matrices))
all_ps = all_ps[!is.na(all_ps)]
par(mfrow=c(1,1));hist(all_ps,main="All p-values", xlab="P-value")
all_ps_thr = max(all_ps[p.adjust(all_ps,method='fdr')<=0.1])
print(length(all_ps));print(all_ps_thr)
pval_matrices_binary = lapply(pval_matrices,function(x,y)x<=y,y=all_ps_thr)

# Look at the most significant p
ppp = 1e-100
ggg = names(which(apply(pval_matrices[[2]]<=ppp,1,any)))[1]
gdata = acute_gene_tables_raw_simpletime[[ggg]]
tp_meta_analysis_results[[2]][ggg,]
gdata = gdata[gdata$time==4,]
gdata = gdata[gdata$tissue=="blood",]
gdata$vi = pmax(gdata$vi,0.01)
gdata = gdata[gdata$training!="control",]
get_subset_forest_plot(gdata,"blood",main="PGC1, acute, muscle, 4h")
rma(yi,vi,1/vi,data=gdata)$pval

tp_vecs = sapply(pval_matrices_binary,function(x)sapply(colnames(x),function(x)strsplit(x,split=';')[[1]][1]))
for(i in 1:length(pval_matrices_binary)){
  tp_v = tp_vecs[[i]]
  m = t(apply(pval_matrices_binary[[i]],1,function(x,y)tapply(x,y,any),y=tp_v))
  if(nrow(m)==1){m=t(m)}
  pval_matrices_binary[[i]] = m
}
pval_selected_genes = lapply(pval_matrices_binary,function(x)rownames(x)[rowSums(x)>0])
sapply(pval_selected_genes,length)

control_filter_tests = lapply(tp_meta_analysis_results,function(x)t(apply(x,1,perform_filter_1_on_tps,0.5)))
for(i in 1:length(control_filter_tests)){
  if(nrow(control_filter_tests[[i]])==1){control_filter_tests[[i]] = t(control_filter_tests[[i]])}
}
filter1_passing_genes = lapply(control_filter_tests,function(x)rownames(x)[rowSums(x)>0])
sapply(filter1_passing_genes,length)

pass_test_matrices = list()
for(i in 1:length(pval_matrices_binary)){
  m1 = pval_matrices_binary[[i]]
  if(nrow(m1)==1){m1=t(m1)}
  m2 = control_filter_tests[[i]]
  if(nrow(m2)==1){m2=t(m2)}
  if(ncol(m2)>1){m2 = m2[,colnames(m1)]}
  m3 = m1 & m2
  if(nrow(m3)==1){m3 = t(m3)}
  pass_test_matrices[[i]] = m3
}
names(pass_test_matrices) = names(pval_matrices_binary)
pvals_filter1_selected_genes = lapply(pass_test_matrices,function(x)rownames(x)[rowSums(x)>0])
pvals_filter1_selected_genes$`longterm,blood` = pval_selected_genes$`longterm,blood`

# For blood longterm, select the top intercept genes
longterm_blood_gene_data = 
  tp_meta_analysis_results[["longterm,blood"]][pvals_filter1_selected_genes$`longterm,blood`,]
longterm_blood_effects_beta_a = t(apply(longterm_blood_gene_data,1,get_beta_a_vec))
longterm_blood_effects_beta_a = apply(longterm_blood_effects_beta_a,1,function(x)max(abs(x)))
pvals_filter1_selected_genes$`longterm,blood` = names(sort(longterm_blood_effects_beta_a,decreasing=T)[1:100])

pvals_filter1_selected_genes_names = lapply(pvals_filter1_selected_genes,function(x,y)unlist(y[x]),y=entrez2symbol)
pval_selected_genes_names = lapply(pval_selected_genes,function(x,y)unlist(y[x]),y=entrez2symbol)
known_genes = c("PPARGC1A","COX1","NDUFA","PDK4","VEGFA","KDR","THY1","MYL4",
                "MYH1","COL1A1","ACTC1","TNNT2","GADD45G","MMP9","NR4A1")
sapply(pvals_filter1_selected_genes_names,intersect,y=known_genes)
sapply(pval_selected_genes_names,intersect,y=known_genes)
pvals_filter1_selected_genes_names[[1]][grepl("^COL",pvals_filter1_selected_genes_names[[1]])]

# compare sizes
sapply(pvals_filter1_selected_genes_names,length)
sapply(pval_selected_genes,length)

# Look at publication bias
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
  print(table(all_ps<0.01)/length(all_ps))
  print(length(all_egger_test_ps))
  # hist(all_egger_test_ps,main="Egger test, acute studies",xlab="P-value")
  # zz = qnorm(all_egger_test_ps)
  # lf = locfdr(zz)
}

# # Manual examinations
# pval_matrices[[1]]["5166",]
# pass_test_matrices[[1]]["5166",]
# control_filter_tests[[1]]["5166",]
# tp_meta_analysis_results[[1]]["7422",]
# gdata = acute_gene_tables_raw[["7422"]]
# gdata = longterm_gene_tables_raw[["5166"]]
# gdata = acute_gene_tables_raw_simpletime[["10891"]]
# gdata = acute_gene_tables_raw[["7139"]]
# rma.mv(yi,vi,1/vi,mods = ~ tissue:factor(time) + training, random = list(~1|V1,~1|gse),data=gdata,btt=2:4)
# rma.mv(yi,vi,1/vi,mods = ~ training + factor(time), 
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
# rma.mv(yi,vi,1/vi,mods = ~is_ctrl,random= ~ 1|gse,data=gdata)
# gdata = acute_gene_tables_raw_simpletime[["10891"]]
# gdata = gdata[gdata$time==4,]
# gdata = gdata[gdata$training!="control",]
# get_subset_forest_plot(gdata,"muscle",main="PGC1, acute, muscle, 4h")

# Large-scale interpretation
# Load all functions from dataset_level_meta_analysis.R
# TODO: move all these to a new functions file

pval_selected_gene_enrichment = run_topgo_enrichment_fisher(pval_selected_genes,
                                union(names(acute_gene_tables),names(longterm_gene_tables)))
extract_top_go_results(pval_selected_gene_enrichment)
get_most_sig_enrichments_by_groups(extract_top_go_results(pval_selected_gene_enrichment))
pvals_filter1_selected_gene_enrichments = run_topgo_enrichment_fisher(pvals_filter1_selected_genes,
                                          union(names(acute_gene_tables),names(longterm_gene_tables)))
res = extract_top_go_results(pvals_filter1_selected_gene_enrichments)
table(res$setname)
get_most_sig_enrichments_by_groups(extract_top_go_results(pvals_filter1_selected_gene_enrichments,0.1,300),1)

# Gene patterns for the analysis
weighted_avg_matrices=list()
weighted_avg_matrices[["acute"]] = t(sapply(acute_gene_tables_raw,get_gene_weighted_avg_pattern))
weighted_avg_matrices[["longterm"]] = t(sapply(longterm_gene_tables_raw,get_gene_weighted_avg_pattern))
weighted_avg_matrices = lapply(weighted_avg_matrices,reorder_weighted_avg_matrix)
# sapply(weighted_avg_matrices,function(x)table(is.na(x)))
# rownames(weighted_avg_matrices[["longterm"]])[rowSums(is.na(weighted_avg_matrices[["longterm"]]))>0]

library(corrplot)
library(gplots)
run_corrmat_clustering_of_a_gene_set<-function(genes,m,exclude_cols_regex = "fat",
                                               min_homogn=0.4,cor_method="spearman",min_clust_size=5){
  m = m[genes,]
  m = m[,!grepl(exclude_cols_regex,colnames(m))]
  mc = cor(t(m),method = "spearman")
  m_clust = cluster_genes_by_homogeneity(mc,min_homogn,kmeans)
  selected_clusters = names(which(table(m_clust)>=min_clust_size))
  m_clust = m_clust[is.element(m_clust,set=selected_clusters)]
  return(list(clustering=m_clust,homogeneities=cluster_homogeneities(mc,m_clust),mc=mc))
}

# cluster our gene sets by their patterns
dir.create("timepoint_meta_analysis_gene_clustering")
pvals_filter1_selected_gene_clusters = list()
for(i in 1:length(pvals_filter1_selected_genes)){
  nn = names(pvals_filter1_selected_genes)[i]
  m = weighted_avg_matrices$acute
  if(grepl("longterm",nn)){m = weighted_avg_matrices$longterm}
  m1 = m[,grepl("muscle",colnames(m))]
  if(grepl("blood",nn)){m1 = m[,grepl("blood",colnames(m))]}
  pvals_filter1_selected_gene_clusters[[nn]] = run_corrmat_clustering_of_a_gene_set(pvals_filter1_selected_genes[[i]],m1)
}
lapply(pvals_filter1_selected_gene_clusters,function(x)table(x[[1]]))
pvals_filter1_clustered_gene_sets=list()
for(i in 1:length(pvals_filter1_selected_genes)){
  nn = names(pvals_filter1_selected_genes)[i]
  m_clust = pvals_filter1_selected_gene_clusters[[nn]][[1]]
  for(j in unique(m_clust)){
    pvals_filter1_clustered_gene_sets[[paste(nn,j,sep="_")]] = names(which(m_clust==j))
  }
}

pvals_filter1_clustered_gene_sets_topgo_res = run_topgo_enrichment_fisher(pvals_filter1_clustered_gene_sets,
                  union(names(acute_gene_tables),names(longterm_gene_tables)))
extract_top_go_results(pvals_filter1_clustered_gene_sets_topgo_res)
get_most_sig_enrichments_by_groups(extract_top_go_results(pvals_filter1_clustered_gene_sets_topgo_res),2)

# cluster plots
names(pvals_filter1_clustered_gene_sets)
m_1 = weighted_avg_matrices$acute[pvals_filter1_clustered_gene_sets[[5]],]
m_1 = weighted_avg_matrices$longterm[pvals_filter1_clustered_gene_sets[[11]],]
m_1[is.na(m_1)|is.nan(m_1)]=0
plot_gene_pattern(apply(m_1,2,mean),errs = apply(m_1,2,sd),tosmooth = T,mfrow=c(2,1),y_lim_add = 0.5,y_lim_min = 1)










