# This script loads the raw profiles and the metadata
# creates standardized datasets and metadata and saves
# it in easy to use objects

###############################################
###############################################
############### Config workspace ##############
###############################################
###############################################
# Paths
WD = "/Users/David/Desktop/MoTrPAC/project_release_feb_2018/data"
SCRIPTS = "/Users/David/Desktop/MoTrPAC/project_release_feb_2018/rcode/"
metadata_file = 'GEO_sample_metadata.xlsx'
raw_data_output_obj = 'human_ge_profiles_db.RData'
raw_data_output_obj = '/Users/David/Desktop/MoTrPAC/PA_database/PA_database_profiles.RData'

# Analysis constants
OMIC_TYPE = "mRNA" # mRNA, methylation
MIN_NROW = 5000 # minimal number of rows in a datamatrix of a specific dataset
OUT_FILE_ACUTE = "human_ge_cohort_preprocessed_db_acute.RData"
OUT_FILE_LONGTERM = "human_ge_cohort_preprocessed_db_longterm.RData"

setwd(WD)
library('xlsx');library(corrplot)
source(paste(SCRIPTS,'ge_download_preprocessing_helper_functions.R',sep=''))
load(raw_data_output_obj)
load('gpl_mappings_to_entrez.RData')

# Comments about the acute metadata
# Time series:
#   -1 means pre-treatment
#   Time is measured in hours
# Training - columns "Training program during experiment type"
#   The type is generally endurance or resistance
#   Untrained == controls
#   There can also be: training with treatment (e.g., LPS)
#   Column "Study subgroup" contains some additional raw information from the study description
# We exclude samples without time info - happens due to some acute/longterm mixed
# datasets such as GSE28392. We also exclude samples without subject ids, and samples 
# whose intervention includes a treatment (very few).
# Study subgroup represetns different training types or treatments within a study - 
# i.e., different interventions within a study.

# Comments about the long-term metadata
# Time series:
#   0 means pre-treatment
#   Time is measured in days
# Training - columns "Training program during experiment type"
#   The type is generally endurance or resistance
#   Untrained == controls
#   There can also be: training with treatment (e.g., LPS), or both endurance and resistance
#   Column "Study subgroup" contains some additional raw information from the study description
# We exclude samples without time info - happens due to some acute/longterm mixed
# datasets such as GSE28392. We also exclude samples without subject ids, and samples 
# whose intervention includes a treatment (very few).
# Study subgroup represetns different training types or treatments within a study - 
# i.e., different interventions within a study.

###############################################
###############################################
########## Preprocessing functions ############
###############################################
###############################################
clean_raw_metadata<-function(metadata){
  metadata = metadata[as.character(metadata[,1])!="",]
  # Solve duplications in the GSMs (happened in one dataset with several GEO ids), only relevant for acute
  gsm_duplications = names(which(table(metadata[,1])>1))
  to_keep = rep(T,nrow(metadata))
  for(curr_gsm in gsm_duplications){
    curr_ids = which(metadata[,1]==curr_gsm)
    if(length(curr_ids)<2){next}
    to_keep[curr_ids[-1]]=F
  }
  metadata = metadata[to_keep,];rownames(metadata) = metadata[,1]
  # Exclude samples without time info - happens due to some acute/longterm mixed datasets such as GSE28392
  stand_time_col = colnames(metadata)[grepl("standard",colnames(metadata),ignore.case = T)&grepl("time",colnames(metadata),ignore.case = T)]
  sample2time = as.numeric(as.character(metadata[,stand_time_col]))
  metadata = metadata[!is.na(sample2time) & sample2time!="",]
  # exclude samples without subject ids - we should contact the authors
  sample2subject = as.character(metadata[,"Subject.id"])
  names(sample2subject) = metadata[,1]
  samples_without_subject = names(which(is.na(sample2subject)|sample2subject==""))
  # get these problematic datasets and remove these samples from the  metadata table
  GSEs_without_subjects = as.character(metadata[samples_without_subject,"GSE"])
  metadata = metadata[!is.element(rownames(metadata),set=samples_without_subject),]  
}
simplify_training_type<-function(metadata){
  raw_training_data = as.character(metadata$Training.program.during.experiment.type)
  dataset_subgroup = as.character(metadata$Study.subgroup)
  sample_is_endurance = grepl("endur",raw_training_data,ignore.case = T)
  sample_is_resistance = grepl("resis",raw_training_data,ignore.case = T)| 
    grepl("strength",raw_training_data,ignore.case = T)
  sample_has_treatment = grepl("treatment",raw_training_data,ignore.case = T)
  # The short description of the training program
  sample2training_type = rep("other",nrow(metadata))
  sample2training_type[sample_is_endurance] = "endurance"
  sample2training_type[sample_is_resistance] = "resistance"
  sample2training_type[sample_is_endurance & sample_is_resistance] = "both"
  sample2training_type[grepl("untrained",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("no training",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("no exercise",raw_training_data,ignore.case = T)] = "untrained"
  sample2training_type[grepl("yoga",raw_training_data,ignore.case = T)] = "yoga"
  sample2training_type[sample_has_treatment] = paste(sample2training_type[sample_has_treatment],"treatment",sep="_")
  sample2training_type[raw_training_data==""] = ""
  names(sample2training_type) = metadata[,1]
  return(sample2training_type)
}
simplify_tissue_info<-function(subject_col2tissue){
  subject_col2tissue_slim = subject_col2tissue
  subject_col2tissue_slim[
    grepl("muscl",subject_col2tissue,ignore.case = T) |
      grepl("vastus",subject_col2tissue,ignore.case = T) |
      grepl("bicep",subject_col2tissue,ignore.case = T)
    ] = "muscle"
  subject_col2tissue_slim[
    grepl("adipo",subject_col2tissue,ignore.case = T) |
      grepl("fat",subject_col2tissue,ignore.case = T)
    ] = "fat"
  subject_col2tissue_slim[
    grepl("blood",subject_col2tissue,ignore.case = T) |
      grepl("cytes",subject_col2tissue,ignore.case = T) |
      grepl("pbmc",subject_col2tissue,ignore.case = T) |
      grepl("phil",subject_col2tissue,ignore.case = T)
    ] = "blood"
  return(subject_col2tissue_slim)
}
simplify_sex_info<-function(sexinfo){
  newv = as.character(sexinfo)
  names(newv) = names(sexinfo)
  newv[grepl(": f",newv)] = "female"
  newv[newv!="" & newv!="female"] = "male"
  newv[newv==""] = NA
  return(newv)
}
get_simplified_sample_information<-function(metadata){
  # The short description of the training program
  sample2training_type = simplify_training_type(metadata)
  # Sanity checks: should be all false
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
  # Time series by sample
  stand_time_col = colnames(metadata)[grepl("standard",colnames(metadata),ignore.case = T)&grepl("time",colnames(metadata),ignore.case = T)]
  sample2time = as.numeric(as.character(metadata[,stand_time_col]))
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
  sample_metadata = list(age=sample2age,sex=sample2sex,replicates=sample2replicate_info,
                         time=sample2time,tissue=sample2tissue,subject=sample2subject,training=sample2training_type,dataset = dataset_ids)
  return(sample_metadata)
}
preprocess_expression_data<-function(metadata,analysis_samples,sample_metadata,
                                     CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes){
  dataset_ids = sample_metadata$dataset
  sample2time = sample_metadata$time
  sample2replicate_info = sample_metadata$replicates
  sample2subject = sample_metadata$subject
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
      data_matrix = log(data_matrix+1,base=2)
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
      gene_fold_changes = get_fold_changes_vs_baseline(genes_data_matrix[,dataset_samples_for_fchange],
                                                       curr_subjects,curr_times[dataset_samples_for_fchange])
    }
    
    # add the results to the data containers
    dataset2preprocessed_data[[dataset]] = list()
    dataset2preprocessed_data[[dataset]][["probe_data"]] = data_matrix
    dataset2preprocessed_data[[dataset]][["gene_data"]] = genes_data_matrix
    dataset2preprocessed_data[[dataset]][["gene_fold_changes"]] = gene_fold_changes
    
    # release memory and save
    rm(data_matrix);rm(genes_data_matrix);rm(genes_data_matrix_obj);gc()
  }
  return(dataset2preprocessed_data)
}
# We currently take the first repeat when subjects have more than
# a single time series data
get_fold_changes_vs_baseline<-function(x,subjs,timev,baseline=NULL,func = function(a,b){a-b},metadata=NULL){
  if(is.null(baseline)){baseline = sort(unique(timev))[1]}
  print(baseline)
  newx = c()
  for(subj in unique(subjs)){
    inds = subjs==subj
    if(sum(inds)<=1){next}
    inds_t = timev[inds]
    currx = x[,inds]
    #print(subj);print(timev[inds])
    base_t_ind = which(inds_t==baseline)
    if(ncol(currx)>2){
      currx_diff = apply(currx[,-base_t_ind],2,func,b=currx[,base_t_ind])
    }
    else{
      currx_diff = as.matrix(func(currx[,-base_t_ind],currx[,base_t_ind]),ncol=1)
    }
    colnames(currx_diff) = paste(subj,inds_t[-base_t_ind],sep="_")
    newx = cbind(newx,currx_diff)
  }
  rownames(newx) = rownames(x)
  return(newx)
}
## Simple functions for cohort-level analyses
# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_yi_vi <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  d = x2-x1
  n = length(d)
  sdd = sd(d)/sqrt(length(d))
  return(c(yi=mean(d),vi=sdd^2))
}
# # test
# x1 = rnorm(10);x2=rnorm(10);x=c(x1,x2)
# s2t = c(rep(0,10),rep(1,10))
# pp =get_paired_ttest_yi_vi(x,s2t,0,1)
# pp
# tt = t.test(x1,x2,paired = T)
# tt$estimate/tt$statistic
# tt$estimate

get_ttest_yi_vi_per_dataset<-function(mat,metadata){
  dataset_times = metadata$time[colnames(mat)]
  if(any(is.na(dataset_times))){
    mat = mat[,!is.na(dataset_times)]
    dataset_times = metadata$time[colnames(mat)]
  }
  min_time = min(dataset_times)
  other_times = setdiff(unique(dataset_times),min_time)
  times2effects = list()
  for(other_time in other_times){
    curr_mat = mat[,dataset_times==min_time | dataset_times==other_time]
    curr_subjects = metadata$subject[colnames(curr_mat)]
    subjects_to_keep = names(which(table(curr_subjects)==2))
    curr_mat = curr_mat[,is.element(curr_subjects,set = subjects_to_keep)]
    ord = order(metadata$subject[colnames(curr_mat)],metadata$time[colnames(curr_mat)])
    curr_mat = curr_mat[,ord]
    curr_times = metadata$time[colnames(curr_mat)]
    paired_test_data = t(apply(curr_mat,1,get_paired_ttest_yi_vi,sample2time=curr_times,t1=min_time,t2=other_time))
    times2effects[[as.character(other_time)]] = paired_test_data
  }
  return(times2effects)
}
get_ttest_pval_per_dataset<-function(mat,metadata){
  dataset_times = metadata$time[colnames(mat)]
  if(any(is.na(dataset_times))){
    mat = mat[,!is.na(dataset_times)]
    dataset_times = metadata$time[colnames(mat)]
  }
  min_time = min(dataset_times)
  other_times = setdiff(unique(dataset_times),min_time)
  times2pvals = c()
  for(other_time in other_times){
    curr_mat = mat[,dataset_times==min_time | dataset_times==other_time]
    curr_subjects = metadata$subject[colnames(curr_mat)]
    subjects_to_keep = names(which(table(curr_subjects)==2))
    curr_mat = curr_mat[,is.element(curr_subjects,set = subjects_to_keep)]
    ord = order(metadata$subject[colnames(curr_mat)],metadata$time[colnames(curr_mat)])
    curr_mat = curr_mat[,ord]
    curr_times = metadata$time[colnames(curr_mat)]
    paired_test_data = apply(curr_mat,1,get_paired_ttest_pval,sample2time=curr_times,t1=min_time,t2=other_time)
    times2pvals = cbind(times2pvals,paired_test_data)
    colnames(times2pvals)[ncol(times2pvals)] = as.character(other_time)
  }
  return(times2pvals)
}
# some simple functions to analyze a time point in a dataset's matrix
get_paired_ttest_pval <-function(x,sample2time,t1,t2){
  x1 = x[sample2time==t1]
  x2 = x[sample2time==t2]
  return(t.test(x1,x2,paired = T)$p.value)
}
##
# Gene selection and display items
get_matrix_p_adjust<-function(x,q=0.1,...){
  v = c(x)
  v = v[!is.na(x)]
  vq = p.adjust(v,...)
  thr = max(v[vq<=q])
  return(thr)
}
###############################################
###############################################
######### Acute data preprocessing ############
###############################################
###############################################

# Sheet 1 has the acute samples metadata
metadata = clean_raw_metadata(read.xlsx2(file=metadata_file,sheetIndex=1))
sample_metadata = get_simplified_sample_information(metadata)

# Select the relevant samples
# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){curr_rows = (metadata[,3] == "RNA" | metadata[,3] == "SRA")}
print(table(curr_rows))
analysis_samples = rownames(metadata)[curr_rows]

dataset2preprocessed_data = preprocess_expression_data(metadata,analysis_samples,sample_metadata,
    CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes)
  
# Encode the dataset ids and create their metadata info
cohort_ids = paste("GE_A_",1:length(dataset2preprocessed_data),sep='')
cohort_info = sapply(names(dataset2preprocessed_data),function(x)strsplit(x,split=';')[[1]])
cohort_metadata = list()
for(j in 1:length(dataset2preprocessed_data)){
  c_id = cohort_ids[j]
  d_id = names(dataset2preprocessed_data)[j]
  gsms = colnames(dataset2preprocessed_data[[d_id]]$gene_data)
  c_info = cohort_info[[j]]
  gse = c_info[1]
  full_tissue = c_info[2]
  gpl = c_info[3]
  training = c_info[4]
  additional_info = NA
  if(length(c_info)>4){additional_info = c_info[5]}
  tissue = simplify_tissue_info(full_tissue)
  cohort_metadata[[c_id]] = list(gsms=gsms,tissue=tissue,training=training,gse=gse,gpl=gpl,
                                 full_tissue=full_tissue,additional_info=additional_info)
  curr_times = sample_metadata$time[gsms]
  if(length(unique(curr_times))<2){next}
  gene_fchanges = get_fold_changes_vs_baseline(dataset2preprocessed_data[[j]]$gene_data,sample_metadata$subject[gsms],curr_times)
  dataset2preprocessed_data[[j]][["gene_fchanges"]] = gene_fchanges
}
names(dataset2preprocessed_data) = cohort_ids
cohort_data = dataset2preprocessed_data
sample2time=sample_metadata$time;sample2sex=sample_metadata$sex;sample2age=sample_metadata$age
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_ACUTE)

### Add fold changes and t-tests #######
# compute ttest p-values, yi's and vi's
for(j in 1:length(cohort_data)){
  if(length(cohort_data[[j]])<3){next}
  res1 = get_ttest_yi_vi_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  res2 = get_ttest_pval_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  for(nn in names(res1)){
    res1[[nn]] = cbind(res1[[nn]],res2[,nn])
    colnames(res1[[nn]]) = c("yi","vi","p")
  }
  cohort_data[[j]][["time2ttest_stats"]] = res1
}
table(sapply(cohort_data,function(x)is.null(x$gene_fold_changes)))
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_ACUTE)


###############################################
###############################################
####### longterm data preprocessing ###########
###############################################
###############################################

# Sheet 1 has the acute samples metadata
metadata = clean_raw_metadata(read.xlsx2(file=metadata_file,sheetIndex=2))
sample_metadata = get_simplified_sample_information(metadata)

# Select the relevant samples
# Get the current metadata
curr_rows = rep(F,nrow(metadata))
names(curr_rows) = rownames(metadata)
# Of the selected rows use the selected omics and exercise types
if(OMIC_TYPE == "mRNA"){curr_rows = (metadata[,3] == "RNA" | metadata[,3] == "SRA")}
print(table(curr_rows))
analysis_samples = rownames(metadata)[curr_rows]

dataset2preprocessed_data = preprocess_expression_data(metadata,analysis_samples,sample_metadata,
      CEL_frma_profiles,CEL_rma_profiles,gse_matrices,gpl_mappings_entrez2probes)

# Encode the dataset ids and create their metadata info
cohort_ids = paste("GE_L_",1:length(dataset2preprocessed_data),sep='')
cohort_info = sapply(names(dataset2preprocessed_data),function(x)strsplit(x,split=';')[[1]])
cohort_metadata = list()
for(j in 1:length(dataset2preprocessed_data)){
  c_id = cohort_ids[j]
  d_id = names(dataset2preprocessed_data)[j]
  gsms = colnames(dataset2preprocessed_data[[d_id]]$gene_data)
  c_info = cohort_info[[j]]
  gse = c_info[1]
  full_tissue = c_info[2]
  gpl = c_info[3]
  training = c_info[4]
  additional_info = NA
  if(length(c_info)>4){additional_info = c_info[5]}
  tissue = simplify_tissue_info(full_tissue)
  cohort_metadata[[c_id]] = list(gsms=gsms,tissue=tissue,training=training,gse=gse,gpl=gpl,
                                 full_tissue=full_tissue,additional_info=additional_info)
  curr_times = sample_metadata$time[gsms]
  if(length(unique(curr_times))<2){next}
  gene_fchanges = get_fold_changes_vs_baseline(dataset2preprocessed_data[[j]]$gene_data,sample_metadata$subject[gsms],curr_times)
  dataset2preprocessed_data[[j]][["gene_fchanges"]] = gene_fchanges
}
names(dataset2preprocessed_data) = cohort_ids
cohort_data = dataset2preprocessed_data
sample2time=sample_metadata$time;sample2sex=sample_metadata$sex;sample2age=sample_metadata$age
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_LONGTERM)

### Add fold changes and t-tests #######
# compute ttest p-values, yi's and vi's
for(j in 1:length(cohort_data)){
  if(length(cohort_data[[j]])<3){next}
  res1 = get_ttest_yi_vi_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  res2 = get_ttest_pval_per_dataset(cohort_data[[j]]$gene_data,sample_metadata)
  for(nn in names(res1)){
    res1[[nn]] = cbind(res1[[nn]],res2[,nn])
    colnames(res1[[nn]]) = c("yi","vi","p")
  }
  cohort_data[[j]][["time2ttest_stats"]] = res1
}
table(sapply(cohort_data,function(x)is.null(x$gene_fold_changes)))
save(sample_metadata,cohort_data,cohort_metadata,sample2time,sample2sex,sample2age,file = OUT_FILE_LONGTERM)





