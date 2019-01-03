
########## !!!!!!!!!! #########
# TODO:
# A major change was performed
# Need to revise this entire code
###############################

# This script loads the PADB and performs 
# meta-analysis of acute studies
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery')
source('helper_functions.R')
GEO_destdir = "GEO"

###############################################
###############################################
############### Set constants #################
###############################################
###############################################

# Analysis constants
# set constants
# Set constants for the analysis
OMIC_TYPE = "mRNA" # mRNA, methylation
EXERCISE_TYPE = "both" # endurance, resistance, both, or other (yoga)
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
OUT_FILE = "PADB_univariate_results_and_preprocessed_data_acute.RData"
ANALYSIS_OUT_FILE = "PADB_statistical_analysis_results_acute.RData"

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
######## Load data and merge studies ##########
###############################################
###############################################
load("PADB_univariate_results_and_preprocessed_data_acute.RData")
load("PADB_sample_metadata_acute.RData")

# Get gene intersection and subject data
gene_sets = lapply(dataset2preprocessed_data,function(x)rownames(x$gene_data))
inters = gene_sets[[1]]
for(l in gene_sets){inters = intersect(inters,l)}
print(length(inters))

# map dataset to numbers
dataset2number = as.character(1:length(dataset2preprocessed_data))
names(dataset2number) = names(dataset2preprocessed_data)

# Get data for the gene fold change analysis
subject_gene_fc_matrix_original = c()
sample_gene_matrix = c()
subject_col2info = list()
subject_info_names = c("sex","age","dataset","subject","tissue","time","training","platform")
#for(nn in subject_info_names){subject_col2info[[nn]]=""}
sample2dataset_number = c()
for(dataset in names(dataset2preprocessed_data)){
  curr_m = dataset2preprocessed_data[[dataset]]$gene_fold_changes
  curr_m2 = dataset2preprocessed_data[[dataset]]$gene_data
  if(is.null(curr_m)){next}
  cols = colnames(curr_m)
  curr_subjects = sapply(cols,function(x)strsplit(x,split="_(\\.|\\d)+$",perl=T)[[1]][1])
  curr_subj_names_with_time = paste(dataset2number[dataset],cols,sep=";")
  names(curr_subjects) = curr_subj_names_with_time
  colnames(curr_m) = curr_subj_names_with_time
  curr_gsms = colnames(dataset2preprocessed_data[[dataset]]$probe_data)
  curr_sex = sapply(curr_subjects,get_subject_info_from_gsms,gsms=curr_gsms,metadata=metadata,
                    sample2subject=sample2subject,colname="Gender")
  curr_times = sapply(cols,get_time_from_subj_names)
  curr_age = sapply(curr_subjects,get_subject_info_from_gsms,gsms=curr_gsms,metadata=metadata,
                    sample2subject=sample2subject,colname="Age")
  curr_dataset_v = rep(dataset,length(cols))
  curr_tissue = sapply(curr_subjects,get_subject_info_from_gsms,gsms=curr_gsms,metadata=metadata,
                    sample2subject=sample2subject,colname="Tissue")
  curr_platform = sapply(curr_subjects,get_subject_info_from_gsms,gsms=curr_gsms,metadata=metadata,
                       sample2subject=sample2subject,colname="platform_id")
  names(curr_dataset_v) = curr_subj_names_with_time
  names(curr_age) = curr_subj_names_with_time
  names(curr_times) = curr_subj_names_with_time
  names(curr_sex) = curr_subj_names_with_time
  names(curr_tissue) = curr_subj_names_with_time
  names(curr_platform) = curr_subj_names_with_time
  
  subject_gene_fc_matrix_original = cbind(subject_gene_fc_matrix_original,curr_m[inters,])
  sample_gene_matrix = cbind(sample_gene_matrix,curr_m2[inters,])
  
  subject_col2info[["time"]] = c(subject_col2info[["time"]],curr_times)
  subject_col2info[["sex"]] = c(subject_col2info[["sex"]],curr_sex)
  subject_col2info[["age"]] = c(subject_col2info[["age"]],curr_age)
  subject_col2info[["dataset"]] = c(subject_col2info[["dataset"]],curr_dataset_v)
  curr_subjects = paste(dataset2number[dataset],curr_subjects,sep=";")
  names(curr_subjects) = curr_subj_names_with_time
  subject_col2info[["subject"]] = c(subject_col2info[["subject"]],curr_subjects)
  subject_col2info[["tissue"]] = c(subject_col2info[["tissue"]],curr_tissue)
  subject_col2info[["platform"]] = c(subject_col2info[["platform"]],curr_platform)
  
  curr_sample2dataset = rep(dataset2number[dataset],ncol(curr_m2))
  names(curr_sample2dataset) = colnames(curr_m2)
  sample2dataset_number = c(sample2dataset_number,curr_sample2dataset)
}
# boxplot(subject_gene_fc_matrix[,sample(1:400)[1:20]],las=2)
# boxplot(sample_gene_matrix[,sample(1:400)[1:20]],las=2)
sapply(subject_col2info,table)
library(preprocessCore)
sample_gene_matrix_quantile = normalize.quantiles.robust(sample_gene_matrix,
      remove.extreme="variance",use.log2=F)
rownames(sample_gene_matrix_quantile) = rownames(sample_gene_matrix)
colnames(sample_gene_matrix_quantile) = colnames(sample_gene_matrix)

# Simplify information
subject_col2info[["tissue"]] = simplify_tissue_info(subject_col2info[["tissue"]])
subject_col2info[["sex"]] = simplify_sex_info(subject_col2info[["sex"]])

# Additional data for the analyses
tissues = unique(subject_col2info[["tissue"]])
sample2subj_with_dataset_number = paste(sample2dataset_number,sample2subject[names(sample2dataset_number)],sep=';')
names(sample2subj_with_dataset_number) = names(sample2dataset_number)
#names(subject_col2info$subject)[which(subject_col2info$subject %in% setdiff(subject_col2info$subject,sample2subj_with_dataset_number))]
# Use the quantile normalized sample matrix to get fold changes
samps = colnames(sample_gene_matrix)
subject_gene_fc_matrix = get_fold_changes_vs_baseline(sample_gene_matrix_quantile,
                                                      sample2subj_with_dataset_number[samps],sample2time[samps])

# sanity checks
all(colnames(subject_gene_fc_matrix) == names(subject_col2info$time)) 
all(colnames(subject_gene_fc_matrix) %in% colnames(subject_gene_fc_matrix_original))
all(dim(subject_gene_fc_matrix)==dim(subject_gene_fc_matrix_original))
for(j in 1:ncol(subject_gene_fc_matrix)){
  print(cor(subject_gene_fc_matrix[,j],subject_gene_fc_matrix_original[,j]))
}

# subject to training 
for (sss in names(subject_col2info$subject)){
  curr_subj = subject_col2info$subject[sss]
  curr_gsms = names(sample2subj_with_dataset_number[sample2subj_with_dataset_number==curr_subj])
  curr_tr = sample2training_type[curr_gsms]
  subject_col2info$training[sss] = curr_tr[1]
}
sapply(subject_col2info,table)

# merge time points for the analysis
subject_cols = colnames(subject_gene_fc_matrix)
subject_col2gse = sapply(subject_col2info$dataset[subject_cols],function(x)strsplit(x,split=';')[[1]][1])
time2num_datasets = table(subject_col2info$time[subject_cols],subject_col2gse)
tp2datasets =  apply(time2num_datasets,1,function(x)names(which(x>0)))
sapply(tp2datasets,length)

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
####### Classification-based analysis #########
###############################################
###############################################

# Classification analysis
library(e1071);library(pROC)
# y labels: sex, tissue, training, time: early vs. later
# Use the FC matrices and the subjects
classification_ys = list();classification_xs = list()
d = subject_col2gse[subject_cols]
# each config has a name, an index to a y vector, and an index to an x matrix
classification_ys[["tissue"]] = as.factor(subject_col2info$tissue[subject_cols])
classification_ys[["sex"]] = as.factor(subject_col2info$sex[subject_cols])
y = subject_col2info$training[subject_cols]
y = y[y!="yoga"]
classification_ys[["training"]] = as.factor(y)
y = subject_col2info$time[subject_cols]
y[y<=1]=0;y[y>1]=1
classification_ys[["time_0_vs_later"]] = as.factor(y)
classification_xs[["fc_orig"]] = t(subject_gene_fc_matrix_original[,subject_cols])
classification_xs[["fc_quantile"]] = t(subject_gene_fc_matrix[,subject_cols])
# covariates for classification
x = cbind(subject_col2info$sex,subject_col2info$tissue,subject_col2info$training)
x_cov = data.frame(x[subject_cols,])
colnames(x_cov) = c("sex","tissue","training")

# Supervised analysis
# A configuration specifies the classification test.
# Each config has the following information in a list:
# The label name
# The features matrix name
# The sample set
# Logical: add covariates or not
# Classifier
# Classifier Arguments
svm_config = list(pred_args = list(probability=T),numFeatures=500,
                  classification_function=svm,probability=T,
                  kernel="linear")
rf_config = list(pred_args = list(type="vote"),numFeatures=500,
                 classification_function=randomForest)
main_configurations = list()
for(yname in names(classification_ys)){
  for(xname in names(classification_xs)){
    y = classification_ys[[yname]]
    y = y[!is.na(y)]
    sample_set = names(y)
    main_configurations[[paste("main",yname,xname,"rf",sep=",")]] = 
        list(yname=yname,xname=xname,include_covs=F,
             sample_set = sample_set,
             classification_args=rf_config)
  }
}
main_configurations_lso_results = run_lso_tests_on_a_configuration_set(main_configurations,
  classification_xs,classification_ys,d,x_cov)
main_configurations_lso_perf_scores = get_standard_classification_performance_scores_for_results_list(
  main_configurations_lso_results
)
main_configurations_lso_subj_acc_scores = get_subject_performance_scores_from_results_list(
  main_configurations_lso_results,subject_col2info$subject,x_cov=x_cov
)

# Additional configurations for time: split by tissue and training type
time_analysis_configs = list()
tissue_and_training = apply(x_cov[,c("tissue","training")],1,paste,collapse=';')
tissue_and_training_dataset_table = table(d[names(tissue_and_training)],tissue_and_training)
subset2num_datasets = colSums(tissue_and_training_dataset_table>0)
selected_new_lables = names(subset2num_datasets)[subset2num_datasets>2]
time_label_name = names(classification_ys)[grepl("time",names(classification_ys))][1]
for(sub_label in selected_new_lables){
  curr_samples = names(tissue_and_training)[tissue_and_training==sub_label]
  curr_table = table(classification_ys[[time_label_name]][curr_samples],d[curr_samples])
  if(any(rowSums(curr_table>0)<2)){next} # skip in lso cannot be run
  l = list(yname=time_label_name,xname=names(classification_xs)[1],include_covs=F,
           sample_set = curr_samples,
           classification_args=rf_config)
  time_analysis_configs[[paste(sub_label,"rf",sep=',')]] = l
}
time_analysis_configs_lso_results = run_lso_tests_on_a_configuration_set(time_analysis_configs,
  classification_xs,classification_ys,d,x_cov)
time_analysis_configs_lso_perf_scores = get_standard_classification_performance_scores_for_results_list(
  time_analysis_configs_lso_results
)
time_analysis_configs_lso_subj_acc_scores = get_subject_performance_scores_from_results_list(
  time_analysis_configs_lso_results,name2subj=subject_col2info$subject,x_cov=x_cov
)

# look at a selected model
x = classification_xs$fc_orig
y = classification_ys$tissue
m = featureSelectionClassifier(x,y,numFeatures = 500,ntree=1000)
rf_tissue = m$classifier
imp_tissue = rf$importance[,1]
y = classification_ys$time_0_vs_later
m = featureSelectionClassifier(x,y,numFeatures = 500,ntree=1000)
rf_time = m$classifier
imp_time = rf$importance[,1]

save(main_configurations,main_configurations_lso_perf_scores,main_configurations_lso_results,main_configurations_lso_subj_acc_scores,
     rf_time,rf_tissue,imp_time,imp_tissue,classification_ys,file="Acute_data_analysis_fc_data_classification_tests.RData")

# Use sample data
classification_ys = list();classification_xs = list()
sample2gse = as.character(metadata$GSE)
names(sample2gse) = metadata[,1]
samps = colnames(sample_gene_matrix_quantile)
d = sample2gse[samps]
# each config has a name, an index to a y vector, and an index to an x matrix
classification_ys[["tissue"]] = as.factor(sample2tissue[samps])
classification_ys[["sex"]] = as.factor(sample2sex[samps])
y = sample2training_type[samps]
y = y[y!="yoga"]
classification_ys[["training"]] = as.factor(y)
y = sample2time[samps]
pre_inds = y==-1
y[pre_inds]=0;y[!pre_inds]=1
classification_ys[["time_0_vs_later"]] = as.factor(y)
classification_xs[["quantile"]] = t(sample_gene_matrix_quantile[,samps])
# covariates for classification
x_cov = data.frame(cbind(sample2sex[samps],sample2tissue[samps],sample2training_type[samps]))
colnames(x_cov) = c("sex","tissue","training")

x_cov_stats = apply(x_cov,1,paste,collapse=';')
table(x_cov_stats)

main_configurations = list()
for(yname in names(classification_ys)){
  for(xname in names(classification_xs)){
    y = classification_ys[[yname]]
    y = y[!is.na(y)]
    sample_set = names(y)
    main_configurations[[paste("main",yname,xname,"rf",sep=",")]] = 
      list(yname=yname,xname=xname,include_covs=F,
           sample_set = sample_set,
           classification_args=rf_config)
  }
}
main_configurations_lso_results = run_lso_tests_on_a_configuration_set(main_configurations,
  classification_xs,classification_ys,d,x_cov)
main_configurations_lso_perf_scores = get_standard_classification_performance_scores_for_results_list(
  main_configurations_lso_results
)
main_configurations_lso_subj_acc_scores = get_subject_performance_scores_from_results_list(
  main_configurations_lso_results,sample2subj,x_cov=x_cov
)

# Additional configurations for time: split by tissue and training type
time_analysis_configs = list()
tissue_and_training = apply(x_cov[,c("tissue","training")],1,paste,collapse=';')
tissue_and_training_dataset_table = table(d[names(tissue_and_training)],tissue_and_training)
subset2num_datasets = colSums(tissue_and_training_dataset_table>0)
selected_new_lables = names(subset2num_datasets)[subset2num_datasets>2]
time_label_name = names(classification_ys)[grepl("time",names(classification_ys))][1]
for(sub_label in selected_new_lables){
  curr_samples = names(tissue_and_training)[tissue_and_training==sub_label]
  curr_table = table(classification_ys[[time_label_name]][curr_samples],d[curr_samples])
  if(any(rowSums(curr_table>0)<2)){next} # skip in lso cannot be run
  l = list(yname=time_label_name,xname=names(classification_xs)[1],include_covs=F,
           sample_set = curr_samples,
           classification_args=rf_config)
  time_analysis_configs[[paste(sub_label,"rf",sep=',')]] = l
}
time_analysis_configs_lso_results = run_lso_tests_on_a_configuration_set(time_analysis_configs,
                                                                         classification_xs,classification_ys,d,x_cov)
time_analysis_configs_lso_perf_scores = get_standard_classification_performance_scores_for_results_list(
  time_analysis_configs_lso_results
)
time_analysis_configs_lso_subj_acc_scores = get_subject_performance_scores_from_results_list(
  time_analysis_configs_lso_results,name2subj=sample2subj,x_cov=x_cov
)

# look at a selected model
x = classification_xs$quantile
y = classification_ys$tissue
m = featureSelectionClassifier(x,y,numFeatures = 500,ntree=1000)
rf_tissue = m$classifier
imp_tissue = rf$importance[,1]
y = classification_ys$time_0_vs_later
m = featureSelectionClassifier(x,y,numFeatures = 500,ntree=1000)
rf_time = m$classifier
imp_time = rf$importance[,1]

save(main_configurations,main_configurations_lso_perf_scores,main_configurations_lso_results,main_configurations_lso_subj_acc_scores,
     rf_time,rf_tissue,imp_time,imp_tissue,classification_ys,file="Acute_data_analysis_expression_intensities_data_classification_tests.RData")

# Display items
load("Acute_data_analysis_expression_intensities_data_classification_tests.RData")
load("Acute_data_analysis_fc_data_classification_tests.RData")
plot(sort(imp_time),ylab="Importance")
sapply(main_configurations_lso_perf_scores,function(x)x$aucs[,3])
sapply(main_configurations_lso_subj_acc_scores,function(x)x$acc)
selected_genes = names(imp_time)[imp_time>2]
sort(imp_time[selected_genes])
top_go_results = run_topgo_enrichment_fisher(selected_genes,rownames(subject_gene_fc_matrix))
extract_top_go_results(top_go_results)

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
##### Fold change replicabilty analysis #######
###############################################
###############################################

# look at the number of datasets per time point
# merge, we want >3 datasets in each cell if possible
time_vs_dataset = table(sample2time[samps],sample2dataset_number[samps])
rowSums(time_vs_dataset>0)
tt = sample2time[samps]
tt[tt==0.5 | tt==1] = 0.75
tt[tt==2.5] = 3
tt[tt==20] = 24
tt[tt==72 | tt==96] = 80
tt[tt==4 | tt==5] = 4.5
tt[tt==48 | tt==80] = 48
time_vs_dataset = table(tt,sample2dataset_number[samps])
rowSums(time_vs_dataset>0)
# Will be used in the mixed effect analysis also

rep_data = get_paired_ttest_matrices(dataset_ids[samps],tt,sample2subj_with_dataset_number,
                                     dataset2preprocessed_data,rownames(sample_gene_matrix))

# Screen
library(kernlab)
source('~/Desktop/screen/supplementary_data/submission_code/SCREEN_code_for_submission.R')
source('~/Desktop/screen/supplementary_data/submission_code/twogroups_methods_for_submission.R')
run_leadingeigen_clustering = function(x,cor_thr=0.2,toplot=F){
  x = x >= cor_thr
  mode(x) = 'numeric';diag(x)=0
  g = graph.adjacency(x,mode='undirected',weighted=T)
  if(toplot){plot(igraph::simplify(g))}
  return (cluster_infomap(g)$membership)
}
dataset_pvals[is.na(dataset_pvals)] = 0.5
screen_res = SCREEN(rep_data$dataset_pvals,ks=2:ncol(rep_data$dataset_pvals),nH=10000)
screen_ind_res = SCREEN_ind(rep_data$dataset_pvals,ks=2:ncol(rep_data$dataset_pvals))
colSums(screen_res<=0.2)
save(rep_data,screen_res,screen_ind_res,file="Acute_replicability_analysis.RData")

###############################################
###############################################
########## Mixed effects analysis #############
###############################################
###############################################

samps = colnames(sample_gene_matrix)
samps = samps[sample2training_type[samps]!="yoga"]
sample2subj = sample2subj_with_dataset_number[samps]
sample2tissue_slim = sample2tissue[samps]
sample2training = sample2training_type[samps]
sample2sex = sample2sex[samps]
table(sample2training)

d1 = data.frame(subjs=factor(sample2subj),timev=tt[samps],
                tissue = sample2tissue_slim,training=sample2training,
                sex=sample2sex,dataset=factor(sample2dataset_number[samps]))


# New analysis: corrected p-values
# Mixed effects with interactions
# Analyze a gene sample to get empirical p-values
gene_sample = sample(1:nrow(sample_gene_matrix))[1:500]
mixed_effect_pvals_emp_pvals = apply(sample_gene_matrix_quantile[gene_sample,samps],1,get_mixed_effect_model_time_empirical_p,
                                     frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                     frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                     d1=d1)
mixed_effect_pvals_apprx_pvals = apply(sample_gene_matrix_quantile[gene_sample,samps],1,get_mixed_effect_model_time_apprx_p,
                                       frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       d1=d1)
hist(mixed_effect_pvals_apprx_pvals);hist(mixed_effect_pvals_emp_pvals)
# Analyze all genes using the standard statistics
mixed_effect_pvals_apprx_statistic = apply(sample_gene_matrix_quantile[,samps],1,get_mixed_effect_model_time_apprx_p,
                                           frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                           frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                           d1=d1,statistic="Chisq")
mixed_effect_pvals_apprx_bic_diff = apply(sample_gene_matrix_quantile[,samps],1,get_mixed_effect_model_time_apprx_stat_diff,
                                          frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                          frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                          d1=d1,statistic="BIC")

original_df = get_stats_df_chisq(mixed_effect_pvals_apprx_statistic[names(mixed_effect_pvals_emp_pvals)],mixed_effect_pvals_apprx_pvals,1:200)
corrected_df = get_stats_df_chisq3(mixed_effect_pvals_apprx_statistic[names(mixed_effect_pvals_emp_pvals)],mixed_effect_pvals_emp_pvals,1:200)
corrected_ps = pchisq(mixed_effect_pvals_apprx_statistic,corrected_df,lower.tail = F)
par(mfrow=c(1,2))
plot(mixed_effect_pvals_apprx_pvals,corrected_ps[names(mixed_effect_pvals_emp_pvals)],ylab="corrected p-values");abline(0,1)
plot(mixed_effect_pvals_emp_pvals,corrected_ps[names(mixed_effect_pvals_emp_pvals)],ylab="corrected p-values");abline(0,1)
par(mfrow=c(1,1))
plot(mixed_effect_pvals_apprx_bic_diff,corrected_ps,ylab="corrected p-values")

save(mixed_effect_pvals_emp_pvals,mixed_effect_pvals_apprx_pvals,mixed_effect_pvals_apprx_statistic,mixed_effect_pvals_apprx_bic_diff,
     original_df, corrected_df,corrected_ps,file="Acute_data_analysis_mixed_effect_analysis.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

# # try lmlist analysis
# lms = lmList(x~timev+timev^2|subjs,data=d)
# lms_coef = coef(lms)
# subject2tissue = unique(cbind(subject_col2subject,simplify_tissue_info(subject_col2tissue)))
# rownames(subject2tissue) = subject2tissue[,1]
# subject2tissue = subject2tissue[,2]
# plot(lms_coef[,2],factor(subject2tissue[rownames(lms_coef)]))

# Cluster genes
all_genes = unique(unlist(set2gene))
study_scaled_matrix = c()
for(dd in unique(d$dataset)){
  curr_samples = rownames(d)[d$dataset==dd]
  currm = t(sample_gene_matrix_quantile[all_genes,curr_samples])
  currm = scale(currm)
  study_scaled_matrix = cbind(study_scaled_matrix,t(currm))
}
dim(study_scaled_matrix)
gene_corrs = cor(t(study_scaled_matrix))
hist(gene_corrs[lower.tri(gene_corrs)])

# Look at a specific gene or a pattern
gene = "4632"
x = sample_gene_matrix_quantile[gene,]
d1 = data.frame(subjs=factor(sample2subj),timev=sample2time[samps],
                tissue = sample2tissue_slim,training=sample2training,
                sex=sample2sex,dataset=factor(sample2dataset_number[samps]))
d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
d_for_design = data.frame(x=x,subjs=factor(sample2subj),timev=ordered(sample2time[samps]),
                          tissue = sample2tissue_slim,training=sample2training,
                          sex=sample2sex,dataset=factor(sample2dataset_number[samps]))
plot.design(d_for_design)

# d_model_matrix_fixed = model.matrix(~1+factor(timev) * factor(tissue) * factor(training) * factor(sex),d1)
# d_model_matrix_rand = model.matrix(~factor(dataset)+factor(subjs),d1)
# #d_model_matrix_fixed = d_model_matrix_fixed[,apply(d_model_matrix_fixed,2,sd)>0]
# d_model_matrix_fixed = d_model_matrix_fixed[,apply(d_model_matrix_fixed,2,sum)>10]
# heatmap.2(d_model_matrix_fixed,trace="none",scale="none")

lmer_est = run_mixed_effect_model_with_drop1_tests(x,d1,F)
lmer_all_interactions1 =  
  lmer(x ~  ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
lmer_all_interactions2 =  
  lmer(x ~  factor(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
lmer_0 =  lmer(x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
lmer_0_2 =  lmer(x ~ factor(tissue) * factor(training) * factor(sex) + (timev|dataset/subjs) ,REML=F,data=d)
anova(lmer_all_interactions1,lmer_0)
anova(lmer_est$lmer_obj_0,lmer_0)
anova(lmer_all_interactions1,lmer_est$lmer_obj)
summary(lmer_all_interactions2)
# Plots of the patterns by groups
# grp1 - the main biological classification
# grp2 - subjects
library(ggplot2)    
grp1 = apply(d[,c("tissue","training","sex")],1,paste,collapse=";")
grp2 = d$dataset
# for gene expression standardize by datasets, no need for the fold changes
x_stand = d$x
for(g2 in unique(grp2)){
  inds = g2==grp2
  x_stand[inds] = scale(x_stand[inds])
}
time_stats = c()
par(mfrow=c(3,3))
for(g1 in sort(unique(grp1))){
  inds = g1==grp1
  currt = d$timev[inds]
  currx = x_stand[inds]
  boxplot(currx~currt,main=g1)
}
ps = list()
for(g1 in sort(unique(grp1))){
  inds = g1==grp1
  currt = d$timev[inds]
  names(currt) = rownames(d)[inds]
  currx = x_stand[inds]
  currg = d$subjs[inds]
  ps[[g1]] = get_profile_longi_line_plot(currx,currt,currg) + 
    geom_line() +ggtitle(g1)
  ord = order(currt)
  currm = tapply(currx[ord],INDEX=as.character(currg[ord]),FUN=function(x)x)
}
multiplot(plotlist=ps,cols=3)

# Plots of the fold changes by groups
# grp1 - the main biological classification
# grp2 - subjects
subject_cols = colnames(subject_gene_fc_matrix)
all(is.element(subject_cols,set=names(subject_col2age)))
curr_fc_pattern = subject_gene_fc_matrix[gene,]
grp1 = cbind(subject_col2info$tissue[subject_cols],
             subject_col2info$training[subject_cols],
             subject_col2info$sex[subject_cols])
grp1 = grp1[!apply(is.na(grp1),1,any),]
grp1 = apply(grp1,1,paste,collapse=';')
grp2 = subject_col2dataset
time_vec = subject_col2info$time[subject_cols]
# for gene expression standardize by datasets, no need for the fold changes
x_stand = subject_gene_fc_matrix[gene,]
time_stats = c()
par(mfrow=c(3,3))
for(g1 in sort(unique(grp1))){
  inds = g1==grp1
  currt = time_vec[inds]
  currx = x_stand[inds]
  boxplot(currx~currt,main=g1,ylim=c(-2,2))
}

ps = list()
for(g1 in sort(unique(grp1))){
  inds = g1==grp1
  currt = time_vec[inds]
  currx = x_stand[inds]
  currg = subject_col2info$subject[names(currt)]
  ps[[g1]] = get_profile_longi_line_plot(currx,currt,currg) + 
    geom_line() +ggtitle(g1)
  ord = order(currt)
  currm = tapply(currx[ord],INDEX=as.character(currg[ord]),FUN=function(x)x)
}
multiplot(plotlist=ps,cols=3)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

anova(lmer_all_interactions1,lmer_est$lmer_obj)
anova(lmer_all_interactions1,lmer_all_interactions2)
par(mfrow=c(1,3))
plot(lmer_est$lmer_obj_0)
plot(lmer_est$lmer_obj)
plot(lmer_all_interactions)
plot(lmer_all_interactions2)
summary(lmer_est$lmer_obj_0)
summary(lmer_est$lmer_obj)
summary(lmer_all_interactions2)
summary(lmer_all_interactions1)
anova(lmer_all_interactions2)
# # residuals by dataset
plot(lmer_est$lmer_obj,form = resid(., type = "pearson") ~ fitted(.) | dataset, abline = 0 )
# longitudinal plots
d1 = subset(d,d$training!="untrained" & d$tissue=="muscle" & d$sex=="male")
d1$timev = as.numeric(as.factor(d1$timev))
d2 = subset(d,d$training=="untrained")
d2$timev = as.numeric(as.factor(d2$timev))
get_profile_longi_line_plot(d2$x,d2$timev,d2$subjs,d2$tissue) + geom_line()
get_profile_longi_line_plot(d1$x,d1$timev,d1$subjs,d1$training) + geom_line()
# # interaction between time and tissue
interaction.plot(d$timev,d$tissue,d$x)
interaction.plot(d$timev,d$sex,d$x)
interaction.plot(d$tissue,d$sex,d$x)
# #####################

# (3) Replicability analysis
# TODO: add this later













