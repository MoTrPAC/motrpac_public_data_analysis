# This script loads the PADB and performs 
# meta-analysis of acute studies
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery');library(corrplot)
source('helper_functions.R')
GEO_destdir = "GEO"

###############################################
###############################################
############# Load the data ###################
###############################################
###############################################

metadata_file = 'GEO_sample_metadata.xlsx'
# Sheet 1 has the acute samples metadata
metadata = read.xlsx2(file=metadata_file,sheetIndex=2)
metadata = metadata[as.character(metadata[,1])!="",]
#gsm_duplications = names(which(table(metadata[,1])>1)) # no duplications
rownames(metadata) = metadata[,1]
print(dim(metadata))

# exclude samples without time info - happens due to some acute/longterm mixed
# datasets such as GSE28392
sample2time = as.numeric(as.character(metadata$Long.term..Standardized.Time..days.))
metadata = metadata[!is.na(sample2time) & sample2time!="",]
print(dim(metadata))
# exclude samples without subject ids - we should contact the authors
sample2subject = as.character(metadata[,"Subject.id"])
names(sample2subject) = metadata[,1]
samples_without_subject = names(which(is.na(sample2subject)|sample2subject==""))
# get these problematic datasets and remove these samples from the 
# metadata table
GSEs_without_subjects = as.character(metadata[samples_without_subject,"GSE"])
table(GSEs_without_subjects)
metadata = metadata[!is.element(rownames(metadata),set=samples_without_subject),]
print(dim(metadata))

# Get sample information
# The short description of the training program
sample2training_type = simplify_training_type(metadata)
table(sample2training_type)
# Sanity checks
table(is.na(sample2training_type)|sample2training_type=="")
table(is.na(metadata$GSE)|metadata$GSE=="")
# Study ids are primarily based on pubmed data
study_ids = as.character(metadata$pmid)
study_ids[study_ids==""] = as.character(metadata[study_ids=="","GSE"])
names(study_ids) = metadata[,1]
# Subject names as given in the original GEO datasets
sample2subject = as.character(metadata[,"Subject.id"])
names(sample2subject) = metadata[,1]
# Dataset id -  a level below study id. It separates samples by:
# GSE id, tissue, platform, training, study subgroup
# The dataset ids contain the above information, separated by ';'
dataset_ids = paste(metadata$GSE,metadata$Tissue,metadata$platform_id,sample2training_type,metadata$Study.subgroup,sep=';')
names(dataset_ids) = metadata[,1]
sort(table(dataset_ids))
# Time series by sample
sample2time = as.numeric(as.character(metadata$Long.term..Standardized.Time..days.))
names(sample2time) = metadata[,1]
table(sample2time)
# Other important data
sample2tissue = simplify_tissue_info(tolower(as.character(metadata$Tissue)))
names(sample2tissue) = metadata[,1]
sample2age = tolower(as.character(metadata$Age))
sample2age = gsub(sample2age,pattern = "age: ",replace="")
names(sample2age) = metadata[,1]
names(sample2tissue) = metadata[,1]
sample2sex = simplify_sex_info(tolower(as.character(metadata$Gender)))
names(sample2sex) = metadata[,1]
sample2replicate_info = as.character(metadata[,"Replicate.info"])
names(sample2replicate_info) = metadata[,1]

# Check for irregularities in subject ids
# Printed tables - do not have the exact same values
# Datasets with some zeroes are okay - these are datasets with
# some missing time points per subject. Not a problem.
for(dataset in unique(dataset_ids)){
  curr_samples = names(which(dataset_ids==dataset))
  curr_samples_rep_info = sample2replicate_info[curr_samples]
  # Keep samples that are not repeats and have times
  curr_samples = curr_samples[curr_samples_rep_info == ""]
  curr_times = sample2time[curr_samples]
  curr_samples = curr_samples[curr_times !=""]
  curr_samples = curr_samples[sample2subject[curr_samples]!=""]
  curr_subjects = sample2subject[curr_samples]
  curr_times = sample2time[curr_samples]
  curr_table = table(curr_subjects,curr_times)
  if(length(curr_table)==0 || is.null(dim(curr_table))){next}
  # Exclude for this test, subjects without the baseline time point
  curr_table = curr_table[curr_table[,1]>0,]
  if(length(curr_table)==0 || is.null(dim(curr_table))){next}
  if(ncol(curr_table)>1){
    curr_table = curr_table[curr_table[,2]>0,]
  }
  if(all(c(curr_table)==c(curr_table)[1])){next}
  print(dataset)
  print(curr_table)
}

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############### Set constants #################
###############################################
###############################################

OMIC_TYPE = "mRNA" # mRNA, methylation
EXERCISE_TYPE = "both" # endurance, resistance, both, or other (yoga)
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
OUT_FILE = "PADB_univariate_results_and_preprocessed_data_longterm.RData"
ANALYSIS_OUT_FILE = "PADB_statistical_analysis_results_longterm.RData"

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
############### Preprocessing #################
###############################################
###############################################
# load the gene expression database
load('gpl_mappings_to_entrez.RData')
load("PA_database_profiles.RData")
load(OUT_FILE)

# Select the relevant samples
# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){curr_rows = (metadata[,3] == "RNA" | metadata[,3] == "SRA")}
print(table(curr_rows))
analysis_samples = rownames(metadata)[curr_rows]

# We start with a simple fold-change exploratory analysis
# The loop below gets the gene and probe matrices for our datasets
dataset2preprocessed_data = list()
for (dataset in unique(dataset_ids[analysis_samples])){
  if(is.element(dataset,set=names(dataset2preprocessed_data))){next}
  assign("last.warning", NULL, envir = baseenv())
  dataset_samples = names(which(dataset_ids==dataset))
  # cut by samples that have time points and are not repeats
  dataset_samples = dataset_samples[!is.na(sample2time[dataset_samples])]
  dataset_samples = dataset_samples[sample2replicate_info[dataset_samples] == ""]
  dataset_samples = dataset_samples[sample2subject[dataset_samples]!=""]
  if(length(dataset_samples)==0){next}
  curr_platforms = as.character(metadata[dataset_samples,"platform_id"])
  platform = curr_platforms[1]
  if(metadata[dataset_samples[1],"Type"]=="SRA"){platform = "recount"}
  
  # for convinience we order the samples by the time
  dataset_samples = dataset_samples[order(sample2time[dataset_samples])]
  
  # Order: fRMA, RMA, GSE object, GSM profiles
  frma_mat = get_data_matrix_from_matrix_lists(CEL_frma_profiles,dataset_samples)
  rma_mat = get_data_matrix_from_matrix_lists(CEL_rma_profiles,dataset_samples)
  gse_mat = get_data_matrix_from_matrix_lists(gse_matrices,dataset_samples)
  print(dim(gse_mat))
  data_matrix = frma_mat
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Not enough samples in fRMA object, use RMA:")
    print(dataset)
    data_matrix = rma_mat
  }
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Not enough samples in RMA object, use GSEs:")
    print(dataset)
    data_matrix = gse_mat
  }
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Not enough samples in GSE object:")
    print(dataset)
    print("Skipping for now, address later")
    next
  }
  na_rows = apply(is.na(data_matrix),1,any) | apply(is.nan(data_matrix),1,any)
  print(table(na_rows))
  data_matrix = data_matrix[!na_rows,]
  if(is.null(data_matrix)|| is.null(dim(data_matrix)) || length(data_matrix)<=1 || ncol(data_matrix)<length(dataset_samples)){
    print("Data matrix has too many rows with NA or NaN values:")
    print(dataset)
    print("Skipping for now, address later")
    next
  }
  if(nrow(data_matrix) < MIN_NROW){
    print("Data matrix has less rows than the minimal number specified:")
    print(dataset)
    print("Skipping")
    next
  }
  
  # Analyze the data matrix
  # 1. Check if the data was logged and correct if needed
  if(max(data_matrix,na.rm=T)>50){
    data_matrix = log(1+data_matrix,base=2)
    data_matrix[is.na(data_matrix)] = min(data_matrix,na.rm=T)
  }
  
  # 2. get current features
  if(length(warnings())>0){print("Got some warnings, Please debug:")}
  
  print("Adding a new dataset to the preprocessed data containers:")
  print(paste("****************",dataset))
  
  # For fold change calculations exclude subjects without a baseline or subjects with baseline only
  curr_times = sample2time[dataset_samples]
  curr_table = table(sample2subject[dataset_samples],curr_times)
  dataset_samples_for_fchange = dataset_samples
  if(is.null(dim(curr_table))||sd(curr_times,na.rm=T)==0){
    dataset_samples_for_fchange = c()
  }
  else{
    subjects_to_remove = rownames(curr_table)[curr_table[,1]==0 | rowSums(curr_table)<=1]
    if(length(subjects_to_remove)>0){
      curr_subjects = sample2subject[dataset_samples]
      samples_to_exclude = dataset_samples[is.element(curr_subjects,set=subjects_to_remove)]
      dataset_samples_for_fchange = setdiff(dataset_samples_for_fchange,samples_to_exclude)
    }
  }

  # Get the data matrices for further analysis
  # Transform to genes (entrez or symbols)
  genes_data_matrix_obj = transform_matrix_into_genes(data_matrix,gpl = platform,gpl_mappings_entrez2probes)
  genes_data_matrix = genes_data_matrix_obj$entrez_mat
  # exclude genes with NAs
  if(sum(genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0)>10000){
    genes_data_matrix = genes_data_matrix[genes_data_matrix_obj$entrez_mat_na_stats[["row NA counts"]]==0,]
  }
  # get the fold change matrices
  gene_fold_changes = NULL; probe_fold_changes = NULL
  if(length(dataset_samples_for_fchange)>0){
    curr_subjects = sample2subject[dataset_samples_for_fchange]
    probe_fold_changes = get_fold_changes_vs_baseline(data_matrix[,dataset_samples_for_fchange],curr_subjects,curr_times[dataset_samples_for_fchange])
    gene_fold_changes = get_fold_changes_vs_baseline(genes_data_matrix[,dataset_samples_for_fchange],curr_subjects,curr_times[dataset_samples_for_fchange])
    table(rowSums(is.na(gene_fold_changes) | is.nan(gene_fold_changes)| is.infinite(gene_fold_changes)))
    table(rowSums(is.na(data_matrix) | is.nan(data_matrix)| is.infinite(data_matrix)))
    corrplot(cor(probe_fold_changes),order="hclust")
  }
  
  # add the results to the data containers
  dataset2preprocessed_data[[dataset]] = list()
  dataset2preprocessed_data[[dataset]][["probe_data"]] = data_matrix
  dataset2preprocessed_data[[dataset]][["gene_data"]] = genes_data_matrix
  dataset2preprocessed_data[[dataset]][["probe_fold_changes"]] = probe_fold_changes
  dataset2preprocessed_data[[dataset]][["gene_fold_changes"]] = gene_fold_changes
  
  # release memory and save
  rm(data_matrix);rm(genes_data_matrix);rm(genes_data_matrix_obj);gc()
  save(dataset2preprocessed_data,file=OUT_FILE)
}

rm(CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes,gpl_mappings_to_entrez,gpl_tables)
gc()

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
load(OUT_FILE)
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
boxplot(sample_gene_matrix[,sample(1:100)[1:20]],las=2)
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
library(e1071);library(pROC)

# Long term analysis: use sample data
classification_ys = list();classification_xs = list()
sample2gse = as.character(metadata$GSE)
names(sample2gse) = metadata[,1]
samps = colnames(sample_gene_matrix_quantile)
d = sample2gse[samps]
# each config has a name, an index to a y vector, and an index to an x matrix
classification_ys[["tissue"]] = as.factor(sample2tissue[samps])
classification_ys[["sex"]] = as.factor(sample2sex[samps])
y = sample2training_type[samps]
classification_ys[["training"]] = as.factor(y)
y = sample2time[samps]
pre_inds = y==0
y[pre_inds]=0;y[!pre_inds]=1
classification_ys[["time_0_vs_later"]] = as.factor(y)
classification_xs[["quantile"]] = t(sample_gene_matrix_quantile[,samps])
# covariates for classification
x_cov = data.frame(cbind(sample2sex[samps],sample2tissue[samps],sample2training_type[samps]))
colnames(x_cov) = c("sex","tissue","training")
x_cov_stats = apply(x_cov,1,paste,collapse=';')
table(x_cov_stats)

# Supervised analysis
# A configuration specifies the classification test.
# Each config has the following information in a list:
# The label name
# The features matrix name
# The sample set
# Logical: add covariates or not
# Classifier
# Classifier Arguments
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
  main_configurations_lso_results,sample2subj_with_dataset_number,x_cov=x_cov
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
     rf_time,rf_tissue,imp_time,imp_tissue,classification_ys,file="Longterm_data_analysis_expression_intensities_data_classification_tests.RData")

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

dataset_pvals = c();dataset_stats = c();dataset_pvals_two_tail=c()
par(mfrow=c(3,4))
for (dataset in unique(dataset_ids)){
  dataset_samples = names(dataset_ids[dataset_ids==dataset])
  dataset_subjects = sample2subj_with_dataset_number[dataset_samples]
  # all NAs == we do not have the data
  if(all(is.na(dataset_subjects))){next}
  tab = table(dataset_subjects)
  paired_subjects = names(tab[tab==2])
  if(length(paired_subjects)==0){next}
  # select subjects with paired samples
  dataset_samples = dataset_samples[is.element(dataset_subjects,set=paired_subjects)]
  dataset_subjects = sample2subj_with_dataset_number[dataset_samples]
  dataset_times = sample2time[dataset_samples]
  # reorder  
  ord = order(dataset_times,dataset_subjects,decreasing=F)
  dataset_samples = dataset_samples[ord]
  dataset_times = dataset_times[ord]
  dataset_subjects = dataset_subjects[ord]
  print(all(dataset_subjects[dataset_times==0] == dataset_subjects[dataset_times>0]))
  # run t-test
  mat = dataset2preprocessed_data[[dataset]]$gene_data[rownames(subject_gene_fc_matrix),dataset_samples]
  ttests = apply(mat,1,run_paired_test,subjs = dataset_subjects,timev=dataset_times,alternative="less")
  ttest_stats = sapply(ttests,function(x)x$statistic)
  ttest_pvals = sapply(ttests,function(x)x$p.value)
  hist(ttest_pvals,main=dataset)
  dataset_pvals = cbind(dataset_pvals,ttest_pvals)
  colnames(dataset_pvals)[ncol(dataset_pvals)] = dataset
  dataset_stats = cbind(dataset_stats,ttest_stats)
  colnames(dataset_stats)[ncol(dataset_stats)] = dataset
  dataset_pvals_two_tail = cbind(dataset_pvals_two_tail,2*pmin(1-ttest_pvals,ttest_pvals))
  colnames(dataset_pvals_two_tail)[ncol(dataset_pvals_two_tail)] = dataset
}

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
screen_res = SCREEN(dataset_pvals)
screen_ind_res = SCREEN_ind(dataset_pvals)

save(dataset_pvals,dataset_stats,dataset_pvals_two_tail,screen_res,screen_ind_res,file="Longterm_replicability_analysis.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
########## Mixed effects analysis #############
###############################################
###############################################

samps = colnames(sample_gene_matrix)
sample2subj = sample2subj_with_dataset_number[samps]
sample2tissue_slim = sample2tissue[samps]
sample2training = sample2training_type[samps]
sample2sex = sample2sex[samps]
table(sample2training)

# look at the number of datasets per time point
# merge, we want >2 datasets in each cell if possible
time_vs_dataset = table(sample2time[samps],sample2dataset_number[samps])
rowSums(time_vs_dataset>0)
tt = sample2time[samps]
tt[tt==70.1 | tt==70.2] = 80
tt[tt==84] = 80
tt[tt>100 & tt<200] = 150
time_vs_dataset = table(tt,sample2dataset_number[samps])
rowSums(time_vs_dataset>0)

d1 = data.frame(subjs=factor(sample2subj),timev=tt,
                    tissue = sample2tissue_slim,training=sample2training,
                    sex=sample2sex,dataset=factor(sample2dataset_number[samps]))


# New analysis: corrected p-values
# Mixed effects with interactions
# Analyze a gene sample to get empirical p-values
gene_sample = sample(1:nrow(sample_gene_matrix))[1:500]
mixed_effect_pvals_emp_pvals = apply(sample_gene_matrix_quantile[gene_sample,],1,get_mixed_effect_model_time_empirical_p,
                                     frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                     frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                     d1=d1)
mixed_effect_pvals_apprx_pvals = apply(sample_gene_matrix_quantile[gene_sample,],1,get_mixed_effect_model_time_apprx_p,
                                       frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       d1=d1)
hist(mixed_effect_pvals_apprx_pvals);hist(mixed_effect_pvals_emp_pvals)
# Empricial vs. approx analysis
plot(mixed_effect_pvals_emp_pvals,mixed_effect_pvals_apprx_pvals);abline(0,1)

# Analyze all genes using the standard statistics
mixed_effect_pvals_apprx_statistic = apply(sample_gene_matrix_quantile,1,get_mixed_effect_model_time_apprx_p,
                                           frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                           frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                           d1=d1,statistic="Chisq")
mixed_effect_pvals_apprx_bic_diff = apply(sample_gene_matrix_quantile,1,get_mixed_effect_model_time_apprx_stat_diff,
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
     original_df, corrected_df,corrected_ps,file="Longterm_data_analysis_mixed_effect_analysis.RData")

###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
######## Results and display items ############
###############################################
###############################################

# Load all the results
load("Longterm_data_analysis_expression_intensities_data_classification_tests.RData")
load("Longterm_data_analysis_mixed_effect_analysis.RData")
load("Longterm_replicability_analysis.RData")

# Display items from the classification analysis
sapply(main_configurations_lso_perf_scores,function(x)x$aucs[,3])
sapply(main_configurations_lso_subj_acc_scores,function(x)x$acc)
selected_genes = names(imp_time)[imp_time>2]
sort(imp_time[selected_genes])
top_go_results = run_topgo_enrichment_fisher(selected_genes,rownames(subject_gene_fc_matrix))
extract_top_go_results(top_go_results)

# Analyze the mixed effects results
table(mixed_effect_pvals_apprx_bic_diff>0)
names(which(mixed_effect_pvals_apprx_bic_diff>0))
corrected_ps_bonf = p.adjust(corrected_ps,method='fdr')
mixed_effects_genes = names(which(corrected_ps_bonf<0.01))
length(mixed_effects_genes)
mixed_effects_genes_topgo = run_topgo_enrichment_fisher(mixed_effects_genes,names(mixed_effect_pvals_apprx_statistic))
extract_top_go_results(mixed_effects_genes_topgo,0.05,500)

pathways = reactomePathways(names(mixed_effect_pvals_apprx_statistic))
#system.time({gsea_res = fgsea(pathways,mixed_effect_pvals_apprx_statistic,nperm=1000,minSize = 5)})
system.time({gsea_res = fgsea_wrapper(pathways,mixed_effect_pvals_apprx_statistic,nperm=20000,minSize = 5)})
qvals = gsea_res$padj; ES = gsea_res$ES
selected_pathways = gsea_res$pathway[qvals<0.01 & ES > 0]
pathway_name = "Mitochondrial translation"
pathway_name = "Hemostasis"
plotEnrichment(pathways[[pathway_name]],mixed_effect_pvals_apprx_statistic)
quantile(mixed_effect_pvals_apprx_statistic[pathways[[pathway_name]]])
quantile(mixed_effect_pvals_apprx_statistic)
qqplot(mixed_effect_pvals_apprx_statistic,mixed_effect_pvals_apprx_statistic[pathways[[pathway_name]]]);abline(0,1)

# Replicability analysis results
# Simple correction, e.g., Shukla et al.
ps = c(dataset_pvals_two_tail)
qs = p.adjust(ps,method = 'BY')
thr = max(ps[qs<0.1],na.rm = T)
bin_ps = dataset_pvals_two_tail<=thr
table(rowSums(bin_ps))
pathways = reactomePathways(rownames(bin_ps))
fgs = fgsea_wrapper(pathways,rowSums(bin_ps))
rep_genes = names(which(rowSums(bin_ps)>=4))
rep_topgo = run_topgo_enrichment_fisher(rep_genes,rownames(bin_ps))
extract_top_go_results(rep_topgo)

colSums(screen_ind_res < 0.2)
rep_genes = rownames(dataset_pvals)[screen_res[,1]<0.2]
rep_topgo = run_topgo_enrichment_fisher(rep_genes,rownames(bin_ps))
extract_top_go_results(rep_topgo)

# Look at the selected gene set
library(gplots)
dataset2simple_name = c()
for(nn in colnames(dataset_stats)){
  samp = names(which(dataset_ids==nn))[1]
  dataset2simple_name[nn] = paste(sample2tissue[samp],sample2training_type[samp])
}
rownames(dataset_stats) = gsub(rownames(dataset_stats),pattern = "\\.t$",replace="",perl=T)
hist(dataset_stats)
xx = dataset_stats[rep_genes,]
xx[xx>6]=6;xx[xx< (-6)] = -6
xx[xx<2 & xx > -2]=0
colnames(xx) = dataset2simple_name
heatmap.2(t(xx),scale="none",trace="none",mar=c(4,10),
          col=colorRampPalette(c("red","white","blue"))(256))

study_two_groups = get_study_marginal_estimation(dataset_pvals)
study_correlations = get_study_pair_corr_matrix(study_two_groups$f_1_mat,study_two_groups$f_0_mat)
colnames(study_correlations) = dataset2simple_name
rownames(study_correlations) = dataset2simple_name
corrplot(study_correlations,order='hclust',bg='gray')

# Simple clustering
xx = dataset_stats[rep_genes,]
colnames(xx) = dataset2simple_name
xx[is.na(xx)|is.nan(xx)]=0
xx[xx>6]=6;xx[xx< (-6)] = -6
xx[xx<2 & xx > -2]=0
xx = xx[!apply(xx==0,1,all),]
dim(xx)
ks = 2:20
kmeans_objs = list()
for(k in ks){kmeans_objs[[as.character(k)]] = kmeans(xx,centers = k)}
totwss = sapply(kmeans_objs,function(x)x$tot.withinss)
plot(totwss)
clust = kmeans_objs[["2"]]$cluster
table(clust)
j = 2
curr_cluster_samples = names(which(clust==j))
heatmap.2(t(xx[curr_cluster_samples,]),scale="none",trace="none",mar=c(4,10),
          col=colorRampPalette(c("red","white","blue"))(256))
sol_as_list = sapply(unique(clust),function(x,y)names(which(y==x)),y=clust)
names(sol_as_list) = paste("Cluster",1:length(sol_as_list))
clustering_enrichment = run_topgo_enrichment_fisher(sol_as_list,rownames(dataset_pvals),5)
extract_top_go_results(clustering_enrichment)


###############################################
###############################################
#################### End ######################
###############################################
###############################################

###############################################
###############################################
######### Additional tests and plots ##########
###############################################
###############################################

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

