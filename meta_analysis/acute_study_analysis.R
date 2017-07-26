# This script loads the PADB and performs 
# meta-analysis of acute studies
setwd('/Users/David/Desktop/MoTrPAC/PA_database')
library('xlsx');library('GEOquery')
source('helper_functions.R')
GEO_destdir = "GEO"

metadata_file = 'GEO_sample_metadata.xlsx'
# Sheet 1 has the acute samples metadata
metadata = read.xlsx2(file=metadata_file,sheetIndex=1)
metadata = metadata[as.character(metadata[,1])!="",]
gsm_duplications = names(which(table(metadata[,1])>1))
# Solve duplications in the GSMs (happened in one dataset with several GEO ids)
to_keep = rep(T,nrow(metadata))
for(curr_gsm in gsm_duplications){
  curr_ids = which(metadata[,1]==curr_gsm)
  to_keep[curr_ids[-1]]=F
}
metadata = metadata[to_keep,]
rownames(metadata) = metadata[,1]

# Get sample information
# The short description of the training program
sample2training_type = simplify_training_type(metadata)
table(sample2training_type)
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
sort(table(dataset_ids))[1:5]
# Time series by sample
sample2time = as.numeric(as.character(metadata$Acute..Standardized.Time..hours...1.is.baseline.))
names(sample2time) = metadata[,1]
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
  if(length(curr_table)==0){next}
  # Exclude for this test, subjects without the baseline time point
  curr_table = curr_table[curr_table[,1]>0,]
  if(all(c(curr_table)==c(curr_table)[1])){next}
  print(dataset)
  print(curr_table)
  #break
}

# Analysis constants
# set constants
# Set constants for the analysis
OMIC_TYPE = "mRNA" # mRNA, methylation
EXERCISE_TYPE = "both" # endurance, resistance, both, or other (yoga)
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
OUT_FILE = "PADB_univariate_results_and_preprocessed_data_acute.RData"
ANALYSIS_OUT_FILE = "PADB_statistical_analysis_results_acute.RData"

########################################################################
######################## Preprocessing #################################
########################################################################

# load the gene expression database
load('gpl_mappings_to_entrez.RData')
load("PA_database_profiles.RData")

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
    data_matrix = log(data_matrix,base=2)
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

########################################################################
######################## Preprocessing end #############################
########################################################################
##### Load without running the loop above
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

#################################
#################################
#################################
# (0) Sanity checks, unsupervised analysis
#################################
#################################
#################################

# Use fold change data
# kmeans
subject_col_clustering = kmeans(scale(t(subject_gene_fc_matrix)),10)
reverse_mapping_list(as.list(subject_col_clustering$cluster))
table(subject_col_clustering$cluster,subject_col2info$sex[subject_cols])
table(subject_col_clustering$cluster,subject_col2info$dataset[subject_cols])
table(subject_col_clustering$cluster,subject_col2gse[subject_cols])
table(subject_col_clustering$cluster)

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

##################################################################
##################################################################
##################################################################
############# Functions for analysis of configurations ###########
##################################################################

run_lso_tests_on_a_configuration_set<-function(configurations,classification_xs,classification_ys,d,x_cov){
  lso_results = list()
  for(cname in names(configurations)){
    config = configurations[[cname]]
    y_names = config$sample_set
    y = classification_ys[[config$yname]][y_names]
    x = classification_xs[[config$xname]][y_names,]
    features_inds_to_keep=c()
    if(config$include_covs=="TRUE"){
      curr_covs = x_cov[,setdiff(colnames(x_cov),config$yname)]
      curr_covs = model.matrix(~.,curr_covs)[,-1]
      new_y_names = intersect(rownames(x),rownames(curr_covs))
      x = x[new_y_names,];y=y[new_y_names]
      x = cbind(curr_covs[new_y_names,],x)
      y_names = new_y_names
      features_inds_to_keep = 1:ncol(curr_covs)
    }
   curr_lso_test_args = c(list(y,x,d[y_names],func=featureSelectionClassifier,
                              features_inds_to_keep=features_inds_to_keep),config$classification_args)
    lso = do.call(leave_study_out2,curr_lso_test_args)
    lso_results[[cname]] = lso
  }
  return(lso_results)
}

get_standard_classification_performance_scores_for_results_list<-function(lso_results){
  performance_scores = list()
  for(cname in names(lso_results)){
    lso = lso_results[[cname]]
    lso_preds = lapply(lso,function(x)x$preds)
    if(! "matrix" %in% class(lso_preds[[1]])){
     lso_preds = lapply(lso,function(x)attr(x$preds,"probabilities"))
   }
   class_names = colnames(lso_preds[[1]])
    m = c()
    for(mm in lso_preds){m = rbind(m,mm[,class_names])}
    colnames(m) = class_names
    lso_preds = m
    aucs = c()
    for(class_name in class_names){
      lso_y = unlist(sapply(lso,function(x)x$real))
      class_inds = lso_y==class_name
      lso_y = rep(0,length(lso_y))
      lso_y[class_inds]=1
      lso_p = lso_preds[,class_name]
      roc = calcAupr(lso_p,as.numeric(lso_y),roc=T)
      proc = as.numeric(pROC::auc(lso_y, lso_p))
      aupr = calcAupr(lso_p,as.numeric(lso_y),roc=F)
      aucs = rbind(aucs,c(roc=roc,aupr=aupr,pROC=proc))
    }
    rownames(aucs) = class_names
    acc = lso_y == (lso_p>=0.5)
    acc = sum(acc)/length(acc)
    performance_scores[[cname]] = list(aucs=aucs,acc=acc)
  }
  return(performance_scores)
}

get_subject_performance_scores_from_results_list<-function(lso_results,name2subj=sample2subj,paired_tests_analysis=T,x_cov){
  subject_performance_scores = list()
  for(cname in names(lso_results)){
    subject_performance_scores[[cname]]=list()
    lso = lso_results[[cname]]
    lso_preds = lapply(lso,function(x)x$preds)
    if(! "matrix" %in% class(lso_preds[[1]])){
      lso_preds = lapply(lso,function(x)attr(x$preds,"probabilities"))
    }
    class_names = colnames(lso_preds[[1]])
    m = c()
    for(mm in lso_preds){m = rbind(m,mm[,class_names])}
    colnames(m) = class_names
    lso_preds = m
    #curr_subjects = subject_col2info$subject[rownames(lso_preds)]
    #curr_subjects = sample2subj[rownames(lso_preds)]
    curr_subjects = name2subj[rownames(lso_preds)]
    lso_y = unlist(sapply(lso,function(x)x$real))
    class1 = colnames(lso_preds)[1]
    if(paired_tests_analysis){
      get_profile_longi_line_plot(lso_preds[,class1],lso_y,curr_subjects) + 
        geom_line() +ggtitle(g1)
     xx1 = c();xx2=c()
     for(subj in unique(curr_subjects)){
       subj_inds = curr_subjects==subj
       if(sum(subj_inds)!=2){next} 
       curr_y = lso_y[subj_inds]
       if(length(unique(curr_y))!=2){next}
       curr_p = lso_preds[subj_inds,class1]
       xx1 = c(xx1,curr_p[curr_y==class1])
       xx2 = c(xx2,curr_p[curr_y!=class1])
     }
     curr_p = NA
     if(length(xx1)>0){curr_p = wilcox.test(xx1,xx2,paired=T)$p.value}
     d1 = data.frame(subjs=factor(curr_subjects),timev=lso_y,
                     tissue = x_cov[rownames(lso_preds),"tissue"],
                     training=x_cov[rownames(lso_preds),"training"],
                     sex=x_cov[rownames(lso_preds),"sex"])
     mixed_effects_model_emp_pval = NA
     # Should not work when splitting the database by one of the covariates
     # TODO: fixed later
     try({ 
      mixed_effects_model_emp_pval = get_mixed_effect_model_time_empirical_p(x=lso_preds[,class1],
          frm0 = x ~ factor(tissue) + factor(training) + factor(sex) + (1|subjs),
          frm1 = x ~ ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|subjs),
          d1=d1,reps=100,min_reps = 100)
     })
     subject_performance_scores[[cname]][["paired_tests_analysis"]] = 
       list(x1=xx1,x2=xx2,binary_wilcox_p=curr_p,mixed_effects_emp_p = mixed_effects_model_emp_pval)
    }
    subj2acc = c()
    for(subj in curr_subjects){
      subj_inds = curr_subjects==subj
      if(sum(subj_inds)<2){next}
      curr_y = lso_y[subj_inds]
      if(!grepl("time",cname) || length(unique(curr_y))==1){
        curr_y = lso_y[subj_inds] == class1
        curr_p = mean(lso_preds[subj_inds,class1])>0.5
        subj2acc[subj] = sum(curr_p==curr_y[1])
      }else{
        curr_p = as.numeric(lso_preds[subj_inds,class1]>0.5)
        subj2acc[subj] = mean(curr_p==curr_y)
      }
    }
    subject_performance_scores[[cname]][["acc"]] = mean(subj2acc,na.rm=T)
    print(paste(cname,mean(subj2acc,na.rm=T)))
  }
  return(subject_performance_scores)
}

##################################################################
##################################################################
##################################################################

#################################
# (1) Simple fold change analysis
#################################
subject_cols = colnames(subject_gene_fc_matrix)
grps = cbind(subject_col2info$tissue[subject_cols],subject_col2info$training[subject_cols],
             subject_col2info$sex[subject_cols],subject_col2info$time[subject_cols])
grps = apply(grps,1,paste,collapse=';')
table(grps)
grp2 = subject_col2info$dataset[subject_cols]

time_vec = subject_col2info$time[subject_cols]
fold_change_sliced_matrices = list()
for(g1 in sort(unique(grps))){
  fold_change_sliced_matrices[[g1]]=list()
  inds = g1==grps
  currt = time_vec[inds]
  currx = subject_gene_fc_matrix_original[,inds]
  fold_change_sliced_matrices[[g1]][["time"]]=unique(currt)
  fold_change_sliced_matrices[[g1]][["mat"]] = currx
}
sort(sapply(fold_change_sliced_matrices,function(x)ncol(x$mat)))
# Analyze a sliced fold change matrix
group2gene_set = list()
geo_mean_thr = 0.5
thr = 1
library(psych)

for(g1 in unique(grps)){
  mat = fold_change_sliced_matrices[[g1]]$mat
  mat_datasets = subject_col2info$dataset[colnames(mat)]
  by_dataset_sum_fc = apply(mat>=thr,1,function(x,y)tapply(x,INDEX=y,FUN=mean),y=mat_datasets)
  #by_dataset_sum_fc = apply(mat,1,function(x,y)tapply(x,INDEX=y,FUN=run_simple_ttest),y=mat_datasets)
  if(is.null(dim(by_dataset_sum_fc))){
    gms=by_dataset_sum_fc
  }
  else{
    gms = apply(by_dataset_sum_fc,2,geometric.mean)
  }
  plot(sort(gms,decreasing = T))
  group2gene_set[[paste(g1,"up",sep=';')]] = names(which(gms>=geo_mean_thr))
  by_dataset_sum_fc = apply(mat<=-thr,1,function(x,y)tapply(x,INDEX=y,FUN=mean),y=mat_datasets)
  if(is.null(dim(by_dataset_sum_fc))){
    gms=by_dataset_sum_fc
  }
  else{
    gms = apply(by_dataset_sum_fc,2,geometric.mean)
  }
  plot(sort(gms,decreasing = T))
  group2gene_set[[paste(g1,"down",sep=';')]] = names(which(gms>=geo_mean_thr))
}
sapply(group2gene_set,length)
topgo_results = run_topgo_enrichment_fisher(group2gene_set,rownames(subject_gene_fc_matrix))

run_simple_ttest<-function(fcs,ret="p.value"){
  if(length(fcs)<2){return(1)}
  return(t.test(fcs)[[ret]])
}

group2paired_ttest_res = list()
for(g1 in unique(grps)){
  mat = fold_change_sliced_matrices[[g1]]$mat
  mat_datasets = subject_col2info$dataset[colnames(mat)]
  dataset_pvals = apply(mat,1,function(x,y)tapply(x,INDEX=y,FUN=run_simple_ttest),y=mat_datasets)
  dataset_stats = apply(mat,1,function(x,y)tapply(x,INDEX=y,FUN=run_simple_ttest,ret="statistic"),y=mat_datasets)
  group2paired_ttest_res[[g1]]=list()
  group2paired_ttest_res[[g1]][["p"]] = dataset_pvals
  group2paired_ttest_res[[g1]][["t"]] = dataset_stats #pretty stupid to run this twice :(
}
# run queries to get gene sets and enrichments
group2num_datasets = sapply(group2paired_ttest_res,function(x)nrow(x$p))
group2num_datasets[sapply(group2num_datasets,is.null)]=1
group2num_datasets = unlist(group2num_datasets)
min_datasets = 2
selected_groups = names(which(group2num_datasets>=min_datasets))
all_ps = unlist(sapply(group2paired_ttest_res[selected_groups],function(x)x$p))
length(all_ps);hist(all_ps)
corrected_ps = p.adjust(all_ps,method='fdr')
hist(corrected_ps)
table(corrected_ps<=0.1)
p_thr = max(all_ps[corrected_ps<=0.1])
group2rep_genes = list()
rep_thr = 0.5
for(g1 in selected_groups){
  currp = group2paired_ttest_res[[g1]]$p
  currt = group2paired_ttest_res[[g1]]$t
  if(is.null(dim(currp))){
    currp = matrix(currp,nrow=1,dimnames = list(g1,names(currp)))
    currt = matrix(currt,nrow=1,dimnames = list(g1,names(currt)))
  }
  print(dim(currp))
  currp_rep = apply(currp <= p_thr,2,mean)
  print(table(currp_rep))
  curr_genes = colnames(currp)[currp_rep>=rep_thr]
  curr_genes_dir = c()
  for(g in curr_genes){
    t_s = currt[,g]
    p_inds = currp[,g]<=p_thr
    t_s = t_s[p_inds]
    curr_genes_dir[g]="mixed"
    if(all(t_s>0)){curr_genes_dir[g]="up"}
    if(all(t_s<0)){curr_genes_dir[g]="down"}
  }
  group2rep_genes[[g1]] = curr_genes_dir
}
sapply(group2rep_genes,table)
up_reg_genes = lapply(group2rep_genes,function(x)names(which(x=="up")))
down_reg_genes = lapply(group2rep_genes,function(x)names(which(x=="down")))
mixed_reg_genes = lapply(group2rep_genes,function(x)names(which(x=="mixed")))
up_reg_genes_go_enrichment = run_topgo_enrichment_fisher(up_reg_genes,rownames(subject_gene_fc_matrix))
extract_top_go_results(up_reg_genes_go_enrichment,0.1)
library(org.Hs.eg.db)
xx <- as.list(org.Hs.egSYMBOL)
up_reg_genes_symb = lapply(up_reg_genes,function(x,y)sort(unlist(y[x])),y=xx)
down_reg_genes_symb = lapply(down_reg_genes,function(x,y)sort(unlist(y[x])),y=xx)
mixed_reg_genes_symb = lapply(mixed_reg_genes,function(x,y)sort(unlist(y[x])),y=xx)

# # Get broader platform information, may be useful later as a random effect
# subject_col2platform_manufacurer = subject_col2platform
# subject_col2platform_technology = subject_col2platform
# for(pl in unique(subject_col2platform)){
#   pl_obj = getGEO(pl,destdir=GEO_destdir)
#   pl_meta = Meta(pl_obj)
#   subject_col2platform_manufacurer[subject_col2platform_manufacurer==pl] = pl_meta$manufacturer
#   subject_col2platform_technology[subject_col2platform_technology==pl] = pl_meta$technology
# }
# table(subject_col2platform_technology)
# table(subject_col2platform_manufacurer)

# # PCA plots
# gene_fold_change_pca = prcomp(t(subject_gene_fc_matrix_quantile),center=T,retx=T)
# gene_fold_change_pca = prcomp(t(subject_gene_fc_matrix),center=T,retx=T)
# plot(gene_fold_change_pca)
# PC1 = gene_fold_change_pca$x[,3]
# PC2 = gene_fold_change_pca$x[,2]
# plot(PC1,PC2,col=as.factor(subject_col2tissue_slim),lwd=2)
# plot(PC1,PC2,col=as.factor(subject_col2platform_manufacurer),lwd=2)
# plot(PC1,PC2,col=as.factor(subject_col2platform_technology),lwd=2)

# The previous step was about preprocessing to get the data matrices.
# We now present a series of analyses:
# (1) Exploratory analysis of the fold changes
# 1.1 Select genes with high fc values
# # 1.2 simple correlation with time
# gene_fc_time_corrs1 = cor(t(subject_gene_fc_matrix),subject_col2time)[,1]
# gene_fc_time_corrs2 = cor(t(subject_gene_fc_matrix_quantile),subject_col2time)[,1]
# samps = colnames(sample_gene_matrix)
# gene_time_corrs1 = cor(t(sample_gene_matrix),sample2time[samps])[,1]
# gene_time_corrs2 = cor(t(sample_gene_matrix_quantile),sample2time[samps])[,1]
# plot(gene_time_corrs1,gene_time_corrs2)
# plot(gene_fc_time_corrs1,gene_fc_time_corrs2)
# plot(gene_time_corrs2,gene_fc_time_corrs2)
# sort(abs(gene_time_corrs2),decreasing=T)[1:10]
# sort(abs(gene_fc_time_corrs2),decreasing=T)[1:10]
# set_nonfc = names(sort(abs(gene_time_corrs2),decreasing=T)[1:50])
# set_fc = names(sort(abs(gene_fc_time_corrs2),decreasing=T)[1:50])
# intersect(set_fc,set_nonfc)
# set_union = union(set_nonfc,set_fc)
# set_union_fc_cors = gene_fc_time_corrs2[set_union]
# set_union_nonfc_cors = gene_time_corrs2[set_union]
# out_matrix = cbind(set_union_nonfc_cors,set_union_fc_cors)
# out_matrix = cbind(is.element(set_union,set=set_nonfc),is.element(set_union,set=set_fc),out_matrix)
# rownames(out_matrix) = set_union
# colnames(out_matrix) = c("is_non_fc","is_fc","nonfc_cor","fc_cor")
# write.table(out_matrix,file="Simple_time_correlation_analysis_selected_genes.txt",sep="\t",quote=F)
# set_nonfc_go_enrichment = run_topgo_enrichment_fisher(set_nonfc,names(gene_time_corrs2))
# extract_top_go_results(set_nonfc_go_enrichment,0.1,500)
# set_fc_go_enrichment = run_topgo_enrichment_fisher(set_fc,names(gene_time_corrs2))
# extract_top_go_results(set_fc_go_enrichment,0.1,500)

################################
# (2) Mixed-effect meta-analysis
################################

# Functions
get_dummy_subject_randomized_vector<-function(times,subjects){
  dummy_times = times
  for(su in unique(subjects)){
    inds = which(subjects==su)
    if(length(inds)<=1){next}
    curr_t = sample(times[inds])
    dummy_times[inds] = curr_t
  }
  return(dummy_times)
}
run_mixed_effect_model1<-function(x,d1,return_pvals=T){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  lmer_obj_0 = lmer(x ~  factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj = lmer(x ~  ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj_int1 = list()
  lmer_obj_int1[["tissue"]] =  lmer(x ~  ordered(timev) * factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj_int1[["training"]] =  lmer(x ~  ordered(timev) * factor(training) + factor(tissue) + factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj_int1[["sex"]] =  lmer(x ~  ordered(timev) * factor(sex) + factor(tissue) + factor(training) + (1|dataset/subjs) ,REML=F,data=d)
  a1 = get_pairwise_anova_aic_bic_p(lmer_obj,lmer_obj_0)
  int_a = sapply(lmer_obj_int1,get_pairwise_anova_aic_bic_p,obj2=lmer_obj)
  int_a_0 = sapply(lmer_obj_int1,get_pairwise_anova_aic_bic_p,obj2=lmer_obj_0)
  if(return_pvals){
    pvals = unlist(c(a1[3],int_a[3,],int_a_0[3,]))
    return(pvals)
  }
  return(list(lmer_obj=lmer_obj,lmer_obj_0=lmer_obj_0,interaction_lmers = lmer_obj_int1,pvals=pvals))
}
run_mixed_effect_model2<-function(x,d1,return_pvals=T,empirical_pval=T,...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  lmer_obj_0 = lmer(x ~  factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  lmer_obj = lmer(x ~  ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs) ,REML=F,data=d)
  a1 = get_pairwise_anova_aic_bic_p(lmer_obj,lmer_obj_0)
  if(empirical_pval){
    rand_scores = get_lmer_model_empirical_null(
      x ~  ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
      d,RMEL=F)
    emp_p = (1+sum(rand_scores[,2]<=a1[2]))/(1+nrow(rand_scores))
    return(emp_p)
  }
  if(return_pvals){return(as.numeric(a1))}
  return(list(lmer_obj=lmer_obj,lmer_obj_0=lmer_obj_0))
}
get_mixed_effect_model_time_empirical_p<-function(x,frm1,frm0,d1,reps=1000,min_reps=200,statistic="Chisq",...){
  min_reps = min(min_reps,reps)
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  rand_scores = c()
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  anova_real = as.numeric(anova(m1,m0)[2,statistic])
  d_copy=d
  for(j in 1:reps){
    d_copy$timev = get_dummy_subject_randomized_vector(d$timev,d$subjs)
    rand_m1 = lmer(frm1,d_copy,REML=F,control=lmerControl(check.rankX =  "silent.drop.cols"))
    rand_m0 = lmer(frm0,d_copy,REML=F,control=lmerControl(check.rankX =  "silent.drop.cols"))
    anova_rand = as.numeric(anova(rand_m1,rand_m0)[2,statistic])
    rand_scores[j] = anova_rand
    if(j>=min_reps){
      emp_p = (1+sum(rand_scores>=anova_real))/(j+1)
      if(emp_p>0.1){break}
      print(emp_p)
    }
  }
  return(emp_p)
}
get_pairwise_anova_aic_bic_p<-function(obj1,obj2){
  return(anova(obj1,obj2)[2,c(2,3,8)])
}
get_mixed_effect_model_time_apprx_p<-function(x,frm1,frm0,d1,statistic="Pr(>Chisq)",...){
  d = cbind(data.frame(x=x,row.names=rownames(d1)),d1)
  m1 = lmer(frm1,d,REML=F)
  m0 = lmer(frm0,d,REML=F)
  return(as.numeric(anova(m1,m0)[2,statistic]))
}

# Analysis starts here
samps = colnames(sample_gene_matrix)
sample2subj = sample2subj_with_dataset_number[samps]
sample2tissue_slim = sample2tissue[samps]
sample2training = sample2training_type[samps]
sample2sex = sample2sex[samps]
table(sample2training)
d1 = data.frame(subjs=factor(sample2subj),timev=sample2time[samps],
                    tissue = sample2tissue_slim,training=sample2training,
                    sex=sample2sex,dataset=factor(sample2dataset_number[samps]))
# look at the number of datasets per time point
time_vs_dataset = table(d1$timev,sample2dataset_number[samps])
rowSums(time_vs_dataset>0)

# Old analysis before p-value correction
mixed_effect_pvals = apply(sample_gene_matrix_quantile,1,run_mixed_effect_model1,d1=d1)
mixed_effect_pvals2 = apply(sample_gene_matrix_quantile,1,run_mixed_effect_model2,d1=d1)
save(mixed_effect_pvals,mixed_effect_pvals2,file='mixed_effect_pvals.RData')

# Save the data to the gene-based mm analysis
save(sample_gene_matrix_quantile,d1,file="acute_data_analysis_data_for_mixed_effects_analysis.RData")

# New analysis: corrected p-values
# Simple fixed effects
mixed_effect_pvals_emp_pvals = apply(sample_gene_matrix_quantile[1:50,],1,get_mixed_effect_model_time_empirical_p,
                                     frm0 = x ~ factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs),
                                     frm1 = x ~ ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs),
                                     d1=d1)
mixed_effect_pvals_apprx_pvals = apply(sample_gene_matrix_quantile[1:50,],1,get_mixed_effect_model_time_apprx_p,
                                       frm0 = x ~ factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs),
                                       frm1 = x ~ ordered(timev) + factor(tissue) + factor(training) + factor(sex) + (1|dataset/subjs),
                                       d1=d1)
# Fixed effects with interactions
gene_sample = sample(1:nrow(sample_gene_matrix))[1:500]
mixed_effect_pvals_emp_pvals = apply(sample_gene_matrix_quantile[gene_sample,],1,get_mixed_effect_model_time_empirical_p,
                                     frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                     frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                     d1=d1)
mixed_effect_pvals_apprx_pvals = apply(sample_gene_matrix_quantile[gene_sample,],1,get_mixed_effect_model_time_apprx_p,
                                       frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       d1=d1)
mixed_effect_pvals_apprx_statistic = apply(sample_gene_matrix_quantile[gene_sample,],1,get_mixed_effect_model_time_apprx_p,
                                       frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                       d1=d1,statistic="Chisq")
get_stats_df_chisq<-function(stats,ps,range,use_log=T,func=abs){
  mse_minus_log = c()
  if(use_log){log_ps = -log(ps,base=10)}
  for(j in range[-1]){
    curr_ps = pchisq(stats,df=j,lower.tail = F)
    curr_ps[curr_ps<min(ps)] = min(ps)
    mse_minus_log[as.character(j)] = mean(func(ps-curr_ps))
    if(use_log){
      curr_ps = -log(curr_ps,base=10)
      mse_minus_log[as.character(j)] = mean(func(log_ps-curr_ps))
    }
    plot(curr_ps,log_ps,main=j);abline(0,1)
  }
  return(range[which(mse_minus_log==min(mse_minus_log))[1]])
}
get_stats_df_chisq2<-function(stats,ps,range){
  comparison_ps = c()
  for(j in range[-1]){
    curr_ps = pchisq(stats,df=j,lower.tail = F)
    curr_ps[curr_ps<min(ps)] = min(ps)
    p = wilcox.test(ps,curr_ps,alternative = "less",paired=T)$p.value
    print(p)
    comparison_ps[j]=p
    if (p < 1e-50){break}
  }
  plot(comparison_ps)
  return(j)
}
get_stats_df_chisq3<-function(stats,ps,range){
  comparison_ps = c()
  for(j in range[-1]){
    curr_ps = pchisq(stats,df=j,lower.tail = F)
    curr_ps[curr_ps<min(ps)] = min(ps)
    if(all(curr_ps>=ps)){break}
  }
  return(j)
}
hist(pchisq(mixed_effect_pvals_apprx_statistic,50,lower.tail = F))
# Empricial vs. approx analysis
plot(mixed_effect_pvals_emp_pvals,mixed_effect_pvals_apprx_pvals);abline(0,1)
original_df = get_stats_df_chisq(mixed_effect_pvals_apprx_statistic[names(mixed_effect_pvals_emp_pvals)],mixed_effect_pvals_apprx_pvals,1:200)
corrected_df = get_stats_df_chisq3(mixed_effect_pvals_apprx_statistic[names(mixed_effect_pvals_emp_pvals)],mixed_effect_pvals_emp_pvals,1:200)

mixed_effect_pvals_apprx_statistic = apply(sample_gene_matrix_quantile,1,get_mixed_effect_model_time_apprx_p,
                                           frm0 = x ~ factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                           frm1 = x ~ ordered(timev) * factor(tissue) * factor(training) * factor(sex) + (1|dataset/subjs),
                                           d1=d1,statistic="Chisq")
write.table(t(t(mixed_effect_pvals_apprx_statistic)),
            file="Acute_mixed_effect_pvals_apprx_statistic.txt",sep="\t",quote=F,col.names = F)


corrected_ps = pchisq(mixed_effect_pvals_apprx_statistic,corrected_df,lower.tail = F)
par(mfrow=c(1,2))
plot(mixed_effect_pvals_apprx_pvals,corrected_ps[names(mixed_effect_pvals_emp_pvals)],ylab="corrected p-values");abline(0,1)
plot(mixed_effect_pvals_emp_pvals,corrected_ps[names(mixed_effect_pvals_emp_pvals)],ylab="corrected p-values");abline(0,1)

corrected_ps_bonf = p.adjust(corrected_ps,method='fdr')
mixed_effects_genes = names(which(corrected_ps_bonf<0.01))
length(mixed_effects_genes)
mixed_effects_genes_topgo = run_topgo_enrichment_fisher(mixed_effects_genes,names(mixed_effect_pvals_apprx_statistic))
extract_top_go_results(mixed_effects_genes_topgo,0.05,500)

# GSEA
library(fgsea)
fgsea_wrapper <- function(pathways,scores,nperm=2000,run_nperm=1000,...){
  num_runs = nperm/run_nperm
  l = list()
  for(j in 1:num_runs){
    l[[j]] = fgsea(pathways,scores,nperm = run_nperm,...)
  }
  emp_pvals = sapply(l,function(x)x$pval)
  emp_pvals = emp_pvals*run_nperm
  min_to_add = min(emp_pvals)
  emp_pvals = emp_pvals-min_to_add
  new_pvals = rowSums(emp_pvals)+min_to_add
  new_pvals = new_pvals/nperm
  new_qvals = p.adjust(new_pvals,method='fdr')
  res = l[[1]]
  res[,"pval"] = new_pvals
  res[,"padj"] = new_qvals
  return(res)
}
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

save(mixed_effect_pvals_apprx_pvals,mixed_effect_pvals_apprx_statistic, corrected_ps, original_df, corrected_df, mixed_effects_genes_topgo,
     mixed_effect_pvals_emp_pvals,file="Acute_data_analysis_empirical_mixed_effect_pval_analysis.RData")

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













